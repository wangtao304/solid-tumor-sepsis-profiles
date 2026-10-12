#!/usr/bin/env Rscript

# Stage 12L-2: first eligible ICU stay sensitivity analysis.
# Frozen MICE/LPA objects are read but never modified or refitted.
# Identifiers and timestamps remain in memory and are never exported.

options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(mice)
  library(survival)
  library(rstpm2)
})

mice_path <- "data/processed/inflammatory_metabolic_mice_solid_m10_maxit20_seed20260721.rds"
lpa_path <- "data/processed/inflammatory_metabolic_lpa_solid_2to6_mi10_seed20260722.rds"
sql_path <- "sql/41_first_eligible_stay_sensitivity_input.sql"
primary_results_path <- "outputs/tables/30_profile_survival_results.csv"
output_path <- "outputs/tables/41_first_eligible_stay_sensitivity.csv"

required_files <- c(mice_path, lpa_path, sql_path, primary_results_path)
if (!all(file.exists(required_files))) {
  stop("One or more required frozen Stage 12L-2 inputs are missing.")
}
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

mi <- readRDS(mice_path)
lpa <- readRDS(lpa_path)
profile_variables <- c(
  "wbc", "platelet_count", "lactate",
  "glucose", "creatinine", "bicarbonate"
)
if (
  !inherits(mi, "mids") ||
  mi$m != 10L ||
  nrow(mi$data) != 4380L ||
  !identical(lpa$configuration$variables, profile_variables) ||
  !identical(lpa$outcome_variables_used, FALSE) ||
  !identical(lpa$survival_information_used, FALSE)
) {
  stop("Frozen MICE/LPA inputs failed safeguard checks.")
}

# Reconstruct the frozen Stage 9A modal consensus without fitting LPA.
assignments <- matrix(NA_integer_, nrow = 4380L, ncol = 10L)
mean_posterior <- matrix(0, nrow = 4380L, ncol = 2L)
for (imputation in 1:10) {
  key <- sprintf("imp%02d_k2", imputation)
  record <- lpa$fit_records[[key]]
  alignment <- lpa$alignments[[key]]
  if (
    is.null(record) ||
    is.null(alignment) ||
    length(alignment) != 2L ||
    length(record$model$classification) != 4380L
  ) {
    stop("Frozen two-profile fit is missing or malformed for ", key, ".")
  }
  assignments[, imputation] <- match(record$model$classification, alignment)
  mean_posterior <- mean_posterior +
    record$model$posterior[, alignment, drop = FALSE] / 10
}
profile_1_votes <- rowSums(assignments == 1L)
profile_2_votes <- rowSums(assignments == 2L)
profile <- ifelse(
  profile_1_votes > profile_2_votes,
  1L,
  ifelse(
    profile_2_votes > profile_1_votes,
    2L,
    max.col(mean_posterior, ties.method = "first")
  )
)
if (
  anyNA(profile) ||
  !identical(as.integer(tabulate(profile, nbins = 2L)), c(3197L, 1183L)) ||
  sum(profile_1_votes == profile_2_votes) != 26L
) {
  stop("Frozen Stage 9A consensus assignment could not be reproduced.")
}

# Stream controlled patient/stay linkage data from PostgreSQL into memory only.
psql_bin <- Sys.getenv(
  "MIMIC_PSQL_BIN",
  unset = ""
)
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
query_text <- system2(psql_bin, args = psql_args, stdout = TRUE, stderr = "")
query_status <- attr(query_text, "status")
if (!is.null(query_status) && query_status != 0L) {
  stop("Read-only first-stay linkage query failed.")
}
db <- read.csv(
  text = paste(query_text, collapse = "\n"),
  check.names = FALSE,
  na.strings = ""
)
rm(query_text)

numeric_equal <- function(x, y, tolerance = 1e-8) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  identical(is.na(x), is.na(y)) &&
    all(abs(x[!is.na(x)] - y[!is.na(y)]) <= tolerance)
}
local <- mi$data
link_checks <- c(
  age = numeric_equal(db$age, local$age),
  sex = identical(as.character(db$sex), as.character(local$sex)),
  cancer_phenotype = identical(
    as.character(db$cancer_phenotype),
    as.character(local$cancer_phenotype)
  ),
  charlson = numeric_equal(
    db$charlson_comorbidity_index,
    local$charlson_comorbidity_index
  ),
  sofa = numeric_equal(db$sofa, local$sofa),
  sapsii = numeric_equal(db$sapsii, local$sapsii)
)

db$event_28d <- as.integer(db$event_28d)
db$time_28d_days <- as.numeric(db$time_28d_days)
db$icu_intime_parsed <- as.POSIXct(
  db$icu_intime,
  format = "%Y-%m-%d %H:%M:%S",
  tz = "UTC"
)

# Hard cohort and linkage gate. No model is fitted before these checks pass.
if (
  nrow(db) != 4380L ||
  length(unique(db$stay_id)) != 4380L ||
  anyDuplicated(db$stay_id) != 0L ||
  length(unique(db$subject_id)) != 3772L ||
  anyNA(db$subject_id) ||
  anyNA(db$stay_id) ||
  anyNA(db$icu_intime_parsed) ||
  anyNA(profile) ||
  anyNA(db$event_28d) ||
  anyNA(db$time_28d_days) ||
  !all(link_checks) ||
  !all(db$event_28d %in% 0:1) ||
  sum(db$event_28d) != 1537L ||
  !all(is.finite(db$time_28d_days)) ||
  !all(db$time_28d_days > 0 & db$time_28d_days <= 28)
) {
  stop("Full-cohort linkage/QC gate failed before first-stay selection.")
}

# Identify the first eligible stay using actual ICU admission time.
minimum_intime <- ave(
  as.numeric(db$icu_intime_parsed),
  db$subject_id,
  FUN = min
)
at_minimum <- as.numeric(db$icu_intime_parsed) == minimum_intime
earliest_counts <- table(db$subject_id[at_minimum])
subjects_with_tied_earliest_time <- sum(earliest_counts > 1L)
if (subjects_with_tied_earliest_time != 0L) {
  stop("At least one subject had tied earliest eligible ICU admission times.")
}

first_order <- order(db$subject_id, db$icu_intime_parsed, db$stay_id)
first_index <- first_order[!duplicated(db$subject_id[first_order])]
first_subjects <- db$subject_id[first_index]
if (
  length(first_index) != 3772L ||
  length(unique(first_subjects)) != 3772L ||
  anyDuplicated(first_subjects) != 0L ||
  anyDuplicated(db$stay_id[first_index]) != 0L ||
  !all(at_minimum[first_index])
) {
  stop("First-eligible-stay cohort QC failed; expected 3,772 unique subjects.")
}

# Preserve the primary analysis parameterization, including full-cohort centering.
full_data <- data.frame(
  profile = factor(profile, levels = 1:2, labels = c("Profile 1", "Profile 2")),
  age10 = (as.numeric(local$age) - mean(as.numeric(local$age))) / 10,
  sex = factor(local$sex),
  charlson = as.numeric(local$charlson_comorbidity_index) -
    mean(as.numeric(local$charlson_comorbidity_index)),
  cancer_phenotype = factor(
    local$cancer_phenotype,
    levels = c("metastatic solid", "nonmetastatic solid")
  ),
  sofa = as.numeric(local$sofa) - mean(as.numeric(local$sofa)),
  sapsii = as.numeric(local$sapsii) - mean(as.numeric(local$sapsii)),
  time = db$time_28d_days,
  event = db$event_28d,
  stringsAsFactors = FALSE
)
first_data <- full_data[first_index, , drop = FALSE]
if (
  nrow(first_data) != 3772L ||
  anyNA(first_data) ||
  sum(table(first_data$profile)) != 3772L
) {
  stop("First-stay model data failed completeness checks.")
}

specification <- data.frame(
  model = c("FS1", "FS2", "FS3", "FS4", "FS-SAPS"),
  rhs = c(
    "profile",
    "profile + age10 + sex",
    "profile + age10 + sex + charlson + cancer_phenotype",
    "profile + age10 + sex + charlson + cancer_phenotype + sofa",
    "profile + age10 + sex + charlson + cancer_phenotype + sapsii"
  ),
  adjustment = c(
    "Profile only",
    "+ age + sex",
    "+ Charlson + cancer phenotype",
    "+ SOFA",
    "SAPS-II substituted for SOFA"
  ),
  stringsAsFactors = FALSE
)
forbidden_adjustments <- c(
  "wbc", "platelet_count", "lactate",
  "glucose", "creatinine", "bicarbonate"
)
if (any(vapply(forbidden_adjustments, function(variable) {
  any(grepl(variable, specification$rhs, fixed = TRUE))
}, logical(1L)))) {
  stop("A profile-defining variable entered a sensitivity model.")
}

profile_counts <- as.integer(table(first_data$profile))
if (length(profile_counts) != 2L || any(profile_counts == 0L)) {
  stop("Both frozen profiles must be represented in the first-stay cohort.")
}

rows <- vector("list", nrow(specification))
for (model_index in seq_len(nrow(specification))) {
  model_formula <- as.formula(
    paste("survival::Surv(time, event) ~", specification$rhs[model_index])
  )
  fit <- stpm2(
    model_formula,
    data = first_data,
    df = 4L,
    scale = "hazard",
    control = list(
      maxit = 10000L,
      use.gr = TRUE,
      reltol.search = 1e-8,
      reltol.final = 1e-10,
      reltol.outer = 1e-8
    )
  )
  fit_details <- fit@details
  profile_term <- grep(
    "^profileProfile 2$",
    names(coef(fit)),
    value = TRUE
  )
  if (
    length(profile_term) != 1L ||
    !identical(as.integer(fit_details$convergence), 0L) ||
    !is.finite(fit_details$maxgrad) ||
    fit_details$maxgrad > 1e-2
  ) {
    stop("Flexible parametric first-stay model failed convergence checks: ",
         specification$model[model_index], ".")
  }
  beta <- coef(fit)[profile_term]
  standard_error <- sqrt(vcov(fit)[profile_term, profile_term])
  z_value <- beta / standard_error
  rows[[model_index]] <- data.frame(
    model = specification$model[model_index],
    N = nrow(first_data),
    deaths = sum(first_data$event),
    profile_1_n = profile_counts[1L],
    profile_2_n = profile_counts[2L],
    HR = exp(beta),
    lower_95_CI = exp(beta - qnorm(0.975) * standard_error),
    upper_95_CI = exp(beta + qnorm(0.975) * standard_error),
    p_value = 2 * pnorm(-abs(z_value)),
    adjustment = specification$adjustment[model_index],
    reference = "Profile 1: lower inflammatory-metabolic stress",
    exposure = "Profile 2: higher inflammatory-metabolic stress",
    baseline_spline_df = 4L,
    scale = "hazard",
    convergence_code = as.integer(fit_details$convergence),
    maximum_absolute_gradient = fit_details$maxgrad,
    aic = AIC(fit),
    stringsAsFactors = FALSE
  )
}
results <- do.call(rbind, rows)
row.names(results) <- NULL

# Reconcile the two adjusted comparisons against the immutable primary table.
primary <- read.csv(primary_results_path, check.names = FALSE)
primary_flexible <- primary[
  primary$analysis_method ==
    "Flexible parametric survival model (PH; 4-df baseline spline)",
  , drop = FALSE
]
primary_sofa <- primary_flexible$estimate[primary_flexible$model == "Model 4"]
primary_saps <- primary_flexible$estimate[
  primary_flexible$model == "Sensitivity: SAPS-II instead of SOFA"
]
if (
  length(primary_sofa) != 1L ||
  length(primary_saps) != 1L ||
  abs(primary_sofa - 1.4317824) > 1e-7 ||
  abs(primary_saps - 1.4286779) > 1e-7
) {
  stop("Frozen primary comparison estimates could not be reconciled.")
}

numeric_columns <- c(
  "HR", "lower_95_CI", "upper_95_CI", "p_value",
  "maximum_absolute_gradient", "aic"
)
results[numeric_columns] <- lapply(results[numeric_columns], signif, 8)
required_columns <- c(
  "model", "N", "deaths", "profile_1_n", "profile_2_n",
  "HR", "lower_95_CI", "upper_95_CI", "p_value"
)
prohibited_output_fields <- c(
  "subject_id", "stay_id", "hadm_id", "profile_assignment",
  "event_28d", "time_28d_days", "icu_intime"
)
if (
  nrow(results) != 5L ||
  !all(required_columns %in% names(results)) ||
  any(names(results) %in% prohibited_output_fields) ||
  !all(results$N == 3772L) ||
  !all(results$deaths == sum(first_data$event)) ||
  !all(results$profile_1_n + results$profile_2_n == 3772L) ||
  !all(results$convergence_code == 0L) ||
  !all(is.finite(results$HR)) ||
  !all(results$lower_95_CI < results$HR & results$HR < results$upper_95_CI)
) {
  stop("Final aggregate first-stay output failed structural verification.")
}

temporary_output <- paste0(output_path, ".tmp")
write.csv(results, temporary_output, row.names = FALSE, na = "")
if (!file.rename(temporary_output, output_path)) {
  unlink(temporary_output)
  stop("Could not finalize the aggregate first-stay output.")
}

cat(sprintf(
  paste0(
    "STAGE12L2_QC full_stays=4380 unique_subjects=3772 first_stays=%d ",
    "duplicate_subjects=0 tied_earliest_subjects=%d deaths=%d profile1=%d ",
    "profile2=%d models=%d no_patient_level_export=TRUE\n"
  ),
  nrow(first_data), subjects_with_tied_earliest_time,
  sum(first_data$event), profile_counts[1L], profile_counts[2L], nrow(results)
))
