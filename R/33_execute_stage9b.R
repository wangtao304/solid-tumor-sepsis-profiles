#!/usr/bin/env Rscript

# Stage 9B: mortality association analysis using the frozen Stage 9A profiles.
# stay_id is used only transiently for linkage QC and is never exported.

options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(mice)
  library(survival)
  library(rstpm2)
})

mice_path <- "data/processed/inflammatory_metabolic_mice_solid_m10_maxit20_seed20260721.rds"
lpa_path <- "data/processed/inflammatory_metabolic_lpa_solid_2to6_mi10_seed20260722.rds"
sql_path <- "sql/30_profile_survival_input.sql"
out30 <- "outputs/tables/30_profile_survival_results.csv"
out31 <- "outputs/tables/31_profile_rmst.csv"
out32 <- "outputs/tables/32_profile_complete_case_sensitivity.csv"
out33 <- "outputs/tables/33_profile_outcome_linkage_qc.csv"
fig_path <- "outputs/figures/figure_profile_km.pdf"
dir.create(dirname(out30), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(fig_path), recursive = TRUE, showWarnings = FALSE)

mi <- readRDS(mice_path)
lpa <- readRDS(lpa_path)
stopifnot(inherits(mi, "mids"), mi$m == 10L, nrow(mi$data) == 4380L)
vars <- c("wbc", "platelet_count", "lactate", "glucose", "creatinine", "bicarbonate")
stopifnot(identical(lpa$configuration$variables, vars),
          identical(lpa$outcome_variables_used, FALSE),
          identical(lpa$survival_information_used, FALSE))

# Reconstruct the previously frozen modal consensus; no LPA is fitted here.
assignments <- matrix(NA_integer_, 4380L, 10L)
posterior <- matrix(0, 4380L, 2L)
for (i in 1:10) {
  key <- sprintf("imp%02d_k2", i)
  rec <- lpa$fit_records[[key]]
  alignment <- lpa$alignments[[key]]
  stopifnot(length(rec$model$classification) == 4380L, length(alignment) == 2L)
  assignments[, i] <- match(rec$model$classification, alignment)
  posterior <- posterior + rec$model$posterior[, alignment, drop = FALSE] / 10
}
v1 <- rowSums(assignments == 1L)
v2 <- rowSums(assignments == 2L)
profile <- ifelse(v1 > v2, 1L, ifelse(v2 > v1, 2L,
  max.col(posterior, ties.method = "first")))
stopifnot(!anyNA(profile), identical(tabulate(profile, 2L), c(3197L, 1183L)), sum(v1 == v2) == 26L)

# Read-only database query; returned identifiers remain in memory only.
psql <- Sys.getenv("MIMIC_PSQL_BIN", "")
pgpass <- Sys.getenv("MIMIC_PGPASS_FILE", "")
db_host <- Sys.getenv("MIMIC_DB_HOST", "")
db_port <- Sys.getenv("MIMIC_DB_PORT", "")
db_user <- Sys.getenv("MIMIC_DB_USER", "")
db_name <- Sys.getenv("MIMIC_DB_NAME", "")
stopifnot(
  file.exists(psql), file.exists(pgpass), file.exists(sql_path),
  all(nzchar(c(db_host, db_port, db_user, db_name)))
)
Sys.setenv(PGPASSFILE = pgpass)
args <- c("-h", db_host,
          "-p", db_port,
          "-U", db_user,
          "-d", db_name,
          "-X", "-w", "-v", "ON_ERROR_STOP=1", "-q", "--csv", "-f", sql_path)
txt <- system2(psql, args, stdout = TRUE, stderr = "")
if (!is.null(attr(txt, "status")) && attr(txt, "status") != 0L) stop("Read-only linkage query failed")
db <- read.csv(text = paste(txt, collapse = "\n"), check.names = FALSE, na.strings = "")
rm(txt)

as_bool <- function(x) tolower(as.character(x)) %in% c("t", "true", "1", "yes")
num_eq <- function(x, y, tol = 1e-8) {
  x <- as.numeric(x); y <- as.numeric(y)
  identical(is.na(x), is.na(y)) && all(abs(x[!is.na(x)] - y[!is.na(y)]) <= tol)
}
local <- mi$data
link_checks <- c(
  age = num_eq(db$age, local$age),
  sex = identical(as.character(db$sex), as.character(local$sex)),
  cancer_phenotype = identical(as.character(db$cancer_phenotype), as.character(local$cancer_phenotype)),
  charlson = num_eq(db$charlson_comorbidity_index, local$charlson_comorbidity_index),
  sofa = num_eq(db$sofa, local$sofa),
  sapsii = num_eq(db$sapsii, local$sapsii),
  complete_case = identical(as_bool(db$complete_case_primary), complete.cases(local[vars]))
)
qc <- data.frame(
  check = c("cohort_row_count", "unique_icu_stay", "duplicate_record_count",
            "profile_assignment_complete", "mortality_28d_complete",
            "follow_up_28d_complete", "covariate_row_order_linkage"),
  observed = c(nrow(db), length(unique(db$stay_id)), sum(duplicated(db$stay_id)),
               sum(!is.na(profile)), sum(!is.na(db$event_28d)), sum(!is.na(db$time_28d_days)),
               sum(link_checks)),
  expected = c(4380L, 4380L, 0L, 4380L, 4380L, 4380L, length(link_checks)),
  status = "PASS",
  notes = c("One row per eligible ICU stay", "stay_id used only in memory",
            "No duplicate linkage records", "Frozen two-profile consensus assignment",
            "Binary 28-day mortality complete", "28-day event/censoring time complete",
            paste(names(link_checks), collapse = "; ")),
  stringsAsFactors = FALSE
)
qc$status[qc$observed != qc$expected] <- "FAIL"
if (any(qc$status != "PASS") || nrow(db) != 4380L || anyDuplicated(db$stay_id)) stop("Profile-outcome linkage QC failed")

db$event_28d <- as.integer(db$event_28d)
db$time_28d_days <- as.numeric(db$time_28d_days)
stopifnot(all(db$event_28d %in% 0:1), sum(db$event_28d) == 1537L,
          all(is.finite(db$time_28d_days)), all(db$time_28d_days > 0), all(db$time_28d_days <= 28))
d <- data.frame(
  profile = factor(profile, 1:2, c("Profile 1", "Profile 2")),
  age10 = (as.numeric(local$age) - mean(as.numeric(local$age))) / 10,
  sex = factor(local$sex),
  charlson = as.numeric(local$charlson_comorbidity_index) - mean(as.numeric(local$charlson_comorbidity_index)),
  cancer_phenotype = factor(local$cancer_phenotype),
  sofa = as.numeric(local$sofa) - mean(as.numeric(local$sofa)),
  sapsii = as.numeric(local$sapsii) - mean(as.numeric(local$sapsii)),
  time = db$time_28d_days,
  event = db$event_28d,
  complete_case = as_bool(db$complete_case_primary)
)
rm(db)
stopifnot(!anyNA(d), sum(d$complete_case) == 2649L)

spec <- data.frame(
  model = c("Model 1", "Model 2", "Model 3", "Model 4", "Sensitivity: SAPS-II instead of SOFA"),
  rhs = c("profile", "profile + age10 + sex",
          "profile + age10 + sex + charlson + cancer_phenotype",
          "profile + age10 + sex + charlson + cancer_phenotype + sofa",
          "profile + age10 + sex + charlson + cancer_phenotype + sapsii"),
  adjustment = c("Profile only", "+ age + sex", "+ Charlson + cancer phenotype",
                 "+ SOFA", "SAPS-II substituted for SOFA"),
  stringsAsFactors = FALSE
)
stopifnot(!any(vapply(vars, function(x) any(grepl(x, spec$rhs, fixed = TRUE)), logical(1))))

extract_effect <- function(fit, term, method, model_row, data, ph_profile = NA_real_, ph_global = NA_real_) {
  b <- coef(fit)[term]
  se <- sqrt(vcov(fit)[term, term])
  data.frame(
    outcome = "28-day mortality", analysis_method = method,
    model = spec$model[model_row], adjustment = spec$adjustment[model_row],
    n = nrow(data), events = sum(data$event), reference = "Profile 1",
    exposure = "Profile 2", effect_measure = "Hazard ratio",
    estimate = exp(b), confidence_interval_lower = exp(b - qnorm(.975) * se),
    confidence_interval_upper = exp(b + qnorm(.975) * se),
    p_value = 2 * pnorm(-abs(b / se)), ph_test_profile_p = ph_profile,
    ph_test_global_p = ph_global, aic = AIC(fit), stringsAsFactors = FALSE
  )
}

fit_set <- function(data) {
  rows <- list()
  for (i in seq_len(nrow(spec))) {
    f <- as.formula(paste("Surv(time, event) ~", spec$rhs[i]))
    fp <- stpm2(f, data = data, df = 4L, scale = "hazard",
      control = list(maxit = 10000L, use.gr = TRUE,
                     reltol.search = 1e-8, reltol.final = 1e-10, reltol.outer = 1e-8))
    term <- grep("^profileProfile 2$", names(coef(fp)), value = TRUE)
    stopifnot(length(term) == 1L, fp@details$convergence == 0L)
    rows[[length(rows) + 1L]] <- extract_effect(fp, term, "Flexible parametric survival model (PH; 4-df baseline spline)", i, data)
    cx <- coxph(f, data = data, ties = "efron", x = TRUE)
    z <- cox.zph(cx, transform = "km", terms = TRUE, singledf = TRUE)
    term_c <- grep("^profileProfile 2$", names(coef(cx)), value = TRUE)
    profile_row <- grep("^profile$", rownames(z$table))
    rows[[length(rows) + 1L]] <- extract_effect(cx, term_c, "Cox proportional hazards sensitivity model", i, data,
      z$table[profile_row, "p"], z$table["GLOBAL", "p"])
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

main <- fit_set(d)
cc <- fit_set(d[d$complete_case, ])

km <- survfit(Surv(time, event) ~ profile, data = d, conf.type = "log-log")
rs <- summary(km, rmean = 28)$table
est <- as.numeric(rs[, "rmean"]); se <- as.numeric(rs[, "se(rmean)"])
rmst <- data.frame(
  outcome = "28-day mortality", horizon_days = 28,
  estimand = c("RMST: Profile 1", "RMST: Profile 2", "RMST difference: Profile 2 minus Profile 1"),
  estimate_days = c(est, est[2] - est[1]),
  standard_error_days = c(se, sqrt(sum(se^2))),
  n = c(as.integer(table(d$profile)), nrow(d)),
  events = c(tapply(d$event, d$profile, sum), sum(d$event)),
  method = "Unadjusted Kaplan-Meier RMST",
  stringsAsFactors = FALSE
)
rmst$confidence_interval_lower_days <- rmst$estimate_days - qnorm(.975) * rmst$standard_error_days
rmst$confidence_interval_upper_days <- rmst$estimate_days + qnorm(.975) * rmst$standard_error_days
rmst$p_value <- c(NA, NA, 2 * pnorm(-abs(rmst$estimate_days[3] / rmst$standard_error_days[3])))

cols <- c("#0072B2", "#D55E00")
pdf(fig_path, width = 8.8, height = 7.2, useDingbats = FALSE)
layout(matrix(c(1, 2), 2), heights = c(4.6, 1.6))
par(mar = c(1.5, 4.5, 3.8, 1), las = 1, xaxs = "i", yaxs = "i")
plot(km, col = cols, lwd = 2.4, conf.int = FALSE, xlim = c(0, 28), ylim = c(0, 1),
     xaxt = "n", xlab = "", ylab = "Survival probability",
     main = "28-day survival by inflammatory-metabolic profile")
axis(1, at = seq(0, 28, 7), labels = FALSE); grid(nx = NA, ny = 5, col = "#E5E5E5")
legend("bottomleft", c(sprintf("Profile 1 (n=%d)", sum(d$profile == "Profile 1")),
                       sprintf("Profile 2 (n=%d)", sum(d$profile == "Profile 2"))),
       col = cols, lwd = 2.4, bty = "n")
times <- seq(0, 28, 7); sm <- summary(km, times = times, extend = TRUE)
risk <- matrix(sm$n.risk, nrow = 2, byrow = TRUE)
par(mar = c(4.5, 4.5, .2, 1), xaxs = "i", yaxs = "i")
plot(NA, xlim = c(-5.0, 30.5), ylim = c(.3, 2.7), axes = FALSE,
     xlab = "Days since ICU admission", ylab = "")
text(-4.7, 2.4, "Number at risk", adj = 0, font = 2, cex = .82)
for (i in 1:2) {
  y <- 2.0 - (i - 1) * .65
  text(-4.7, y, paste("Profile", i), adj = 0, col = cols[i], cex = .80)
  text(times, y, risk[i, ], cex = .86)
}
axis(1, at = times)
dev.off()

round_cols <- c("estimate", "confidence_interval_lower", "confidence_interval_upper", "p_value",
                "ph_test_profile_p", "ph_test_global_p", "aic")
main[round_cols] <- lapply(main[round_cols], signif, 8)
cc[round_cols] <- lapply(cc[round_cols], signif, 8)
rmst[vapply(rmst, is.numeric, logical(1))] <- lapply(rmst[vapply(rmst, is.numeric, logical(1))], signif, 8)
write.csv(main, out30, row.names = FALSE, na = "")
write.csv(rmst, out31, row.names = FALSE, na = "")
write.csv(cc, out32, row.names = FALSE, na = "")
write.csv(qc, out33, row.names = FALSE, na = "")

# Final aggregate-only safeguards.
stopifnot(nrow(main) == 10L, nrow(cc) == 10L, nrow(rmst) == 3L, nrow(qc) == 7L,
          all(qc$status == "PASS"), all(main$n == 4380L), all(cc$n == 2649L),
          !any(c("subject_id", "stay_id", "profile_assignment") %in% names(main)),
          file.info(fig_path)$size > 1000)
cat(sprintf("STAGE9B_OK full_n=%d cc_n=%d deaths=%d\n", nrow(d), sum(d$complete_case), sum(d$event)))
