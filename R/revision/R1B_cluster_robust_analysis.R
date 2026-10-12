#!/usr/bin/env Rscript

# Reviewer Revision R1B: authorized subject-clustered robust-variance analysis.
# Identifiers remain in protected memory and only aggregate model comparisons
# are written. MICE, LPA, consensus labels, outcomes, and covariates are frozen.

options(stringsAsFactors = FALSE, warn = 1)

if (!identical(Sys.getenv("R1B_EXECUTE_CLUSTER_ROBUST"), "YES")) {
  cat("R1B cluster-robust execution gate is closed.\n")
  quit(save = "no", status = 0L)
}

suppressPackageStartupMessages({
  library(mice)
  library(survival)
  library(rstpm2)
})

if (packageVersion("rstpm2") != package_version("1.7.1")) {
  stop("This script was reviewed for rstpm2 1.7.1; another version requires re-review.")
}

mice_path <- "data/processed/inflammatory_metabolic_mice_solid_m10_maxit20_seed20260721.rds"
lpa_path <- "data/processed/inflammatory_metabolic_lpa_solid_2to6_mi10_seed20260722.rds"
sql_path <- "sql/41_first_eligible_stay_sensitivity_input.sql"
frozen_results_path <- "outputs/tables/30_profile_survival_results.csv"
output_path <- "outputs/revision/R1B_cluster_robust_results.csv"
diagnostic_path <- "outputs/revision/R1B_cluster_robust_diagnostics.csv"
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

required <- c(mice_path, lpa_path, sql_path, frozen_results_path)
if (!all(file.exists(required))) stop("A required frozen input is missing.")
if (file.exists(output_path) || file.exists(diagnostic_path)) {
  stop("Refusing to overwrite an existing R1B cluster-robust output.")
}

profile_variables <- c(
  "wbc", "platelet_count", "lactate", "glucose", "creatinine", "bicarbonate"
)
mi <- readRDS(mice_path)
lpa <- readRDS(lpa_path)
if (
  !inherits(mi, "mids") || mi$m != 10L || nrow(mi$data) != 4380L ||
  !identical(lpa$configuration$variables, profile_variables) ||
  !identical(lpa$outcome_variables_used, FALSE) ||
  !identical(lpa$survival_information_used, FALSE)
) stop("Frozen MICE/LPA safeguard check failed.")

# Reconstruct the frozen consensus assignment without fitting LPA.
assignment_votes <- matrix(NA_integer_, nrow = 4380L, ncol = 10L)
mean_posterior <- matrix(0, nrow = 4380L, ncol = 2L)
for (imputation in seq_len(10L)) {
  key <- sprintf("imp%02d_k2", imputation)
  record <- lpa$fit_records[[key]]
  alignment <- lpa$alignments[[key]]
  if (is.null(record) || is.null(alignment) || length(alignment) != 2L) {
    stop("Malformed frozen LPA record: ", key)
  }
  assignment_votes[, imputation] <- match(record$model$classification, alignment)
  mean_posterior <- mean_posterior +
    record$model$posterior[, alignment, drop = FALSE] / 10
}
vote_1 <- rowSums(assignment_votes == 1L)
vote_2 <- rowSums(assignment_votes == 2L)
profile <- ifelse(
  vote_1 > vote_2, 1L,
  ifelse(vote_2 > vote_1, 2L, max.col(mean_posterior, ties.method = "first"))
)
if (!identical(as.integer(tabulate(profile, nbins = 2L)), c(3197L, 1183L))) {
  stop("Frozen profile consensus could not be reproduced.")
}

# Read-only deterministic linkage. No identifier is written.
psql_bin <- Sys.getenv("MIMIC_PSQL_BIN", unset = "")
pgpass_file <- Sys.getenv(
  "MIMIC_PGPASS_FILE",
  unset = ""
)
db_host <- Sys.getenv("MIMIC_DB_HOST", "")
db_port <- Sys.getenv("MIMIC_DB_PORT", "")
db_user <- Sys.getenv("MIMIC_DB_USER", "")
db_name <- Sys.getenv("MIMIC_DB_NAME", "")
if (!file.exists(psql_bin) || !file.exists(pgpass_file) ||
    any(!nzchar(c(db_host, db_port, db_user, db_name)))) {
  stop("Set all MIMIC database environment variables to authorized local values.")
}
Sys.setenv(PGPASSFILE = pgpass_file)
psql_args <- c(
  "-h", db_host,
  "-p", db_port,
  "-U", db_user,
  "-d", db_name,
  "-X", "-w", "-v", "ON_ERROR_STOP=1", "-q", "--csv", "-f", sql_path
)
query_text <- system2(psql_bin, psql_args, stdout = TRUE, stderr = "")
if (!is.null(attr(query_text, "status")) && attr(query_text, "status") != 0L) {
  stop("Read-only R1B cluster linkage query failed.")
}
db <- utils::read.csv(
  text = paste(query_text, collapse = "\n"), check.names = FALSE, na.strings = ""
)
rm(query_text)

numeric_equal <- function(x, y, tolerance = 1e-8) {
  x <- as.numeric(x); y <- as.numeric(y)
  identical(is.na(x), is.na(y)) &&
    all(abs(x[!is.na(x)] - y[!is.na(y)]) <= tolerance)
}

local <- mi$data
link_checks <- c(
  age = numeric_equal(db$age, local$age),
  sex = identical(as.character(db$sex), as.character(local$sex)),
  cancer_phenotype = identical(
    as.character(db$cancer_phenotype), as.character(local$cancer_phenotype)
  ),
  charlson = numeric_equal(db$charlson_comorbidity_index, local$charlson_comorbidity_index),
  sofa = numeric_equal(db$sofa, local$sofa),
  sapsii = numeric_equal(db$sapsii, local$sapsii)
)

db$event_28d <- as.integer(db$event_28d)
db$time_28d_days <- as.numeric(db$time_28d_days)
if (
  nrow(db) != 4380L || length(unique(db$stay_id)) != 4380L ||
  length(unique(db$subject_id)) != 3772L || anyDuplicated(db$stay_id) ||
  anyNA(db[c("subject_id", "stay_id", "event_28d", "time_28d_days")]) ||
  !all(link_checks) || sum(db$event_28d) != 1537L ||
  any(!is.finite(db$time_28d_days)) || any(db$time_28d_days <= 0) ||
  any(db$time_28d_days > 28)
) stop("R1B cluster linkage/cohort gate failed; no model was fitted.")

d <- data.frame(
  subject_id = db$subject_id,
  profile = factor(profile, 1:2, c("Profile 1", "Profile 2")),
  age10 = (as.numeric(local$age) - mean(as.numeric(local$age))) / 10,
  sex = factor(local$sex),
  charlson = as.numeric(local$charlson_comorbidity_index) -
    mean(as.numeric(local$charlson_comorbidity_index)),
  cancer_phenotype = factor(local$cancer_phenotype),
  sofa = as.numeric(local$sofa) - mean(as.numeric(local$sofa)),
  sapsii = as.numeric(local$sapsii) - mean(as.numeric(local$sapsii)),
  time = db$time_28d_days,
  event = db$event_28d,
  check.names = FALSE
)
rm(db, local)
if (anyNA(d)) stop("Prepared analysis data contain missing values.")

spec <- data.frame(
  model = c(
    "Model 1", "Model 2", "Model 3", "Model 4",
    "Sensitivity: SAPS-II instead of SOFA"
  ),
  rhs = c(
    "profile",
    "profile + age10 + sex",
    "profile + age10 + sex + charlson + cancer_phenotype",
    "profile + age10 + sex + charlson + cancer_phenotype + sofa",
    "profile + age10 + sex + charlson + cancer_phenotype + sapsii"
  ),
  adjustment = c(
    "Profile only", "+ age + sex", "+ Charlson + cancer phenotype",
    "+ SOFA", "SAPS-II substituted for SOFA"
  ),
  stringsAsFactors = FALSE
)

control_settings <- list(
  maxit = 10000L, use.gr = TRUE,
  reltol.search = 1e-8, reltol.final = 1e-10, reltol.outer = 1e-8
)

extract_term <- function(fit) {
  term <- grep("^profileProfile 2$", names(stats::coef(fit)), value = TRUE)
  if (length(term) != 1L || fit@details$convergence != 0L) {
    stop("A model did not converge or lacked the frozen profile term.")
  }
  beta <- unname(stats::coef(fit)[term])
  variance <- unname(stats::vcov(fit)[term, term])
  if (!is.finite(beta) || !is.finite(variance) || variance <= 0) {
    stop("A model returned an invalid profile coefficient or variance.")
  }
  c(beta = beta, se = sqrt(variance))
}

results <- vector("list", nrow(spec))
diagnostics <- vector("list", nrow(spec))
z <- stats::qnorm(0.975)

for (i in seq_len(nrow(spec))) {
  formula_i <- stats::as.formula(paste("survival::Surv(time, event) ~", spec$rhs[i]))
  model_matrix <- stats::model.matrix(stats::delete.response(stats::terms(formula_i)), d)

  ordinary <- rstpm2::stpm2(
    formula_i, data = d, df = 4L, scale = "hazard", control = control_settings
  )
  clustered <- rstpm2::stpm2(
    formula_i, data = d, df = 4L, scale = "hazard", control = control_settings,
    robust = TRUE, cluster = "subject_id", frailty = FALSE
  )
  o <- extract_term(ordinary)
  r <- extract_term(clustered)
  # rstpm2 1.7.1 stores the optimized negative log-likelihood in @min;
  # its stpm2 class does not expose a working stats::logLik method here.
  ll_o <- -as.numeric(ordinary@min)
  ll_r <- -as.numeric(clustered@min)
  coefficient_difference <- unname(r["beta"] - o["beta"])
  loglik_difference <- ll_r - ll_o
  if (
    !is.finite(ll_o) || !is.finite(ll_r) ||
    abs(coefficient_difference) > 1e-8 || abs(loglik_difference) > 1e-8
  ) stop("Robust refit changed the coefficient or likelihood beyond numerical tolerance.")

  results[[i]] <- data.frame(
    model = spec$model[i],
    adjustment = spec$adjustment[i],
    N_stays = nrow(d),
    N_subjects = length(unique(d$subject_id)),
    death_records = sum(d$event),
    HR = exp(o["beta"]),
    original_SE_logHR = o["se"],
    cluster_robust_SE_logHR = r["se"],
    original_CI_lower = exp(o["beta"] - z * o["se"]),
    original_CI_upper = exp(o["beta"] + z * o["se"]),
    cluster_robust_CI_lower = exp(r["beta"] - z * r["se"]),
    cluster_robust_CI_upper = exp(r["beta"] + z * r["se"]),
    cluster_robust_P = 2 * stats::pnorm(-abs(r["beta"] / r["se"])),
    SE_ratio = r["se"] / o["se"],
    pooling_method = "Single frozen consensus assignment; no Rubin pooling",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  diagnostics[[i]] <- data.frame(
    model = spec$model[i],
    n_model_matrix_rows = nrow(model_matrix),
    n_model_matrix_columns = ncol(model_matrix),
    ordinary_convergence = ordinary@details$convergence,
    clustered_convergence = clustered@details$convergence,
    ordinary_log_likelihood = ll_o,
    clustered_log_likelihood = ll_r,
    log_likelihood_difference = loglik_difference,
    logHR_difference = coefficient_difference,
    robust = TRUE,
    cluster_unit = "subject_id",
    frailty = FALSE,
    random_effect = FALSE,
    weighting = "none",
    baseline_spline_df = 4L,
    scale = "hazard",
    stringsAsFactors = FALSE
  )
}

result <- do.call(rbind, results)
diagnostic <- do.call(rbind, diagnostics)
rownames(result) <- NULL
rownames(diagnostic) <- NULL

frozen <- utils::read.csv(frozen_results_path, check.names = FALSE)
frozen <- frozen[grepl("Flexible parametric", frozen$analysis_method), ]
frozen_match <- frozen[match(result$model, frozen$model), ]
if (
  nrow(frozen_match) != nrow(spec) || anyNA(frozen_match$estimate) ||
  any(abs(result$HR - frozen_match$estimate) > 1e-5)
) stop("Ordinary refit did not reproduce the frozen primary estimates; output withheld.")

result$frozen_primary_HR <- frozen_match$estimate
result$HR_difference_from_frozen <- result$HR - result$frozen_primary_HR
numeric_result <- vapply(result, is.numeric, logical(1L))
result[numeric_result] <- lapply(result[numeric_result], signif, digits = 10L)
numeric_diag <- vapply(diagnostic, is.numeric, logical(1L))
diagnostic[numeric_diag] <- lapply(diagnostic[numeric_diag], signif, digits = 12L)

for (x in list(result, diagnostic)) {
  if (any(grepl("subject_id|stay_id|hadm_id", names(x), ignore.case = TRUE))) {
    stop("Identifier field detected in aggregate output.")
  }
}
utils::write.csv(result, output_path, row.names = FALSE, na = "")
utils::write.csv(diagnostic, diagnostic_path, row.names = FALSE, na = "")
cat("R1B_CLUSTER_ROBUST_COMPLETE aggregate outputs written; no identifiers exported.\n")
