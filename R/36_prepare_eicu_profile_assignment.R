#!/usr/bin/env Rscript

# Stage 10C-1: eICU-specific MICE and frozen MIMIC profile transfer.
# No mortality, outcome, survival, discharge, or patient-identifier field is read.
# Only a protected local object and aggregate QC tables are finalized.

options(stringsAsFactors = FALSE, warn = 1)

if (!requireNamespace("mice", quietly = TRUE)) {
  stop("Required package 'mice' is not installed; no installation is authorized.")
}

stage_seed <- 20260724L
expected_n <- 1407L
primary_variables <- c(
  "wbc", "platelet_count", "lactate",
  "glucose", "creatinine", "bicarbonate"
)
auxiliary_variables <- c("age", "sex", "cancer_phenotype")
expected_variables <- c(primary_variables, auxiliary_variables)
prohibited_patterns <- c(
  "patientunitstayid", "patienthealthsystemstayid", "uniquepid",
  "subject_id", "stay_id", "hadm_id", "identifier",
  "mortality", "death", "expired", "outcome", "survival",
  "discharge", "follow_up", "followup", "censor", "los"
)

sql_path <- file.path("sql", "36_eicu_profile_assignment_input.sql")
transfer_path <- file.path(
  "data", "processed",
  "inflammatory_metabolic_profile_transfer_logistic_mi10_seed20260723.rds"
)
object_path <- file.path(
  "data", "processed",
  "eicu_inflammatory_metabolic_profile_transfer_mi10_seed20260724.rds"
)
assignment_qc_path <- file.path(
  "outputs", "tables", "36_eicu_profile_assignment_qc.csv"
)
imputation_qc_path <- file.path(
  "outputs", "tables", "37_eicu_imputation_qc.csv"
)

for (path in c(sql_path, transfer_path)) {
  if (!file.exists(path)) stop("Missing required Stage 10C-1 input: ", path)
}
for (path in c(object_path, assignment_qc_path, imputation_qc_path)) {
  if (file.exists(path)) stop("Stage 10C-1 output already exists; refusing to overwrite: ", path)
}
for (path in c(object_path, assignment_qc_path, imputation_qc_path)) {
  if (!dir.exists(dirname(path))) stop("Controlled output directory is missing: ", dirname(path))
}

contains_prohibited_name <- function(x) {
  any(vapply(prohibited_patterns, function(pattern) {
    any(grepl(pattern, x, ignore.case = TRUE))
  }, logical(1L)))
}

# Validate the frozen MIMIC transfer object before any eICU query.
transfer_model <- readRDS(transfer_path)
if (
  !identical(transfer_model$stage, "Stage 10B") ||
  !identical(transfer_model$variables, primary_variables) ||
  !identical(dim(transfer_model$coefficients), c(7L, 10L)) ||
  !identical(dim(transfer_model$centers), c(6L, 10L)) ||
  !identical(dim(transfer_model$scales), c(6L, 10L)) ||
  !identical(rownames(transfer_model$centers), primary_variables) ||
  !identical(rownames(transfer_model$scales), primary_variables) ||
  !identical(as.numeric(transfer_model$classification_threshold), 0.5) ||
  !identical(transfer_model$eicu_data_used_for_training, FALSE) ||
  !identical(transfer_model$mortality_or_survival_used, FALSE) ||
  !identical(transfer_model$lpa_refitted, FALSE)
) {
  stop("The Stage 10B transfer object does not match the frozen specification.")
}
if (
  any(!is.finite(transfer_model$coefficients)) ||
  any(!is.finite(transfer_model$centers)) ||
  any(!is.finite(transfer_model$scales)) ||
  any(transfer_model$scales <= 0)
) {
  stop("The frozen transfer object contains invalid parameters.")
}

# Stream only the nine approved non-outcome columns. The SQL uses the ICU-stay
# key only for deterministic ordering and never returns it. A stdin mode allows
# the host shell to supply the password-file context without exposing or copying
# credentials into the workspace; patient-level CSV remains an in-memory pipe.
if (identical(Sys.getenv("STAGE10C1_INPUT_STDIN", unset = "0"), "1")) {
  csv_lines <- readLines(file("stdin"), warn = FALSE)
} else {
  psql_bin <- Sys.getenv(
    "MIMIC_PSQL_BIN",
    unset = ""
  )
  pgpass_file <- Sys.getenv(
    "MIMIC_PGPASS_FILE",
    unset = ""
  )
  db_host <- Sys.getenv("EICU_DB_HOST", unset = "")
  db_port <- Sys.getenv("EICU_DB_PORT", unset = "")
  db_user <- Sys.getenv("EICU_DB_USER", unset = "")
  db_name <- Sys.getenv("EICU_DB_NAME", unset = "")
  if (!file.exists(psql_bin) || !file.exists(pgpass_file) ||
      any(!nzchar(c(db_host, db_port, db_user, db_name)))) {
    stop("Set all eICU database environment variables to authorized local values.")
  }
  Sys.setenv(PGPASSFILE = pgpass_file)
  psql_args <- c(
    "-h", db_host,
    "-p", db_port,
    "-U", db_user,
    "-d", db_name,
    "-X", "-w", "-v", "ON_ERROR_STOP=1", "-q", "--csv", "-f", sql_path
  )
  csv_lines <- system2(psql_bin, args = psql_args, stdout = TRUE, stderr = "")
  psql_status <- attr(csv_lines, "status")
  if (!is.null(psql_status) && psql_status != 0L) {
    stop("PostgreSQL failed while streaming the Stage 10C-1 input.")
  }
}
if (length(csv_lines) < 2L) stop("The Stage 10C-1 query returned no data rows.")

analysis_data <- utils::read.csv(
  text = paste(csv_lines, collapse = "\n"),
  header = TRUE,
  na.strings = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
rm(csv_lines)

if (!identical(names(analysis_data), expected_variables)) {
  stop("Input columns do not exactly match the frozen nine-variable order.")
}
if (contains_prohibited_name(names(analysis_data))) {
  stop("A prohibited identifier, outcome, discharge, or survival field entered the input.")
}
if (nrow(analysis_data) != expected_n) {
  stop("The eICU target row count does not reconcile to the frozen feasibility cohort.")
}

analysis_data[c(primary_variables, "age")] <- lapply(
  analysis_data[c(primary_variables, "age")],
  as.numeric
)
analysis_data$sex <- factor(analysis_data$sex, levels = c("F", "M"))
analysis_data$cancer_phenotype <- factor(
  analysis_data$cancer_phenotype,
  levels = c("metastatic solid", "nonmetastatic solid")
)
if (anyNA(analysis_data[auxiliary_variables])) {
  stop("An approved eICU auxiliary variable is missing or unexpectedly coded.")
}
if (any(!is.finite(as.matrix(analysis_data[c(primary_variables, "age")][
  !is.na(analysis_data[c(primary_variables, "age")])
])))) {
  stop("A numeric input contains a nonfinite observed value.")
}
for (variable in primary_variables) {
  observed <- analysis_data[[variable]][!is.na(analysis_data[[variable]])]
  if (!length(observed) || any(!is.finite(observed)) || any(observed < 0)) {
    stop("Invalid observed values for transfer variable: ", variable)
  }
}

original_missing_n <- colSums(is.na(analysis_data))
if (!all(original_missing_n[primary_variables] > 0L)) {
  stop("Every frozen eICU imputation target was expected to contain missing values.")
}
if (any(original_missing_n[auxiliary_variables] != 0L)) {
  stop("Predictor-only auxiliaries must be complete under the frozen plan.")
}

method <- setNames(rep("", length(expected_variables)), expected_variables)
method[primary_variables] <- "pmm"
predictor_matrix <- matrix(
  0L,
  nrow = length(expected_variables),
  ncol = length(expected_variables),
  dimnames = list(expected_variables, expected_variables)
)
predictor_matrix[primary_variables, expected_variables] <- 1L
diag(predictor_matrix) <- 0L
if (
  !all(method[primary_variables] == "pmm") ||
  any(nzchar(method[auxiliary_variables])) ||
  any(diag(predictor_matrix) != 0L) ||
  any(predictor_matrix[auxiliary_variables, , drop = FALSE] != 0L)
) {
  stop("Frozen Stage 10C-1 method vector or predictor matrix is invalid.")
}

imputation <- mice::mice(
  data = analysis_data,
  m = 10L,
  method = method,
  predictorMatrix = predictor_matrix,
  maxit = 20L,
  printFlag = FALSE,
  seed = stage_seed,
  donors = 5L,
  remove.constant = FALSE,
  remove.collinear = FALSE
)

gelman_rhat <- function(chain_matrix) {
  chain_matrix <- as.matrix(chain_matrix)
  if (nrow(chain_matrix) < 2L || ncol(chain_matrix) < 2L || anyNA(chain_matrix)) {
    return(NA_real_)
  }
  within_variance <- mean(apply(chain_matrix, 2L, stats::var))
  between_variance <- nrow(chain_matrix) * stats::var(colMeans(chain_matrix))
  if (!is.finite(within_variance) || within_variance < 0) return(NA_real_)
  if (within_variance == 0) return(if (between_variance == 0) 1 else Inf)
  variance_hat <- (
    (nrow(chain_matrix) - 1) / nrow(chain_matrix)
  ) * within_variance + between_variance / nrow(chain_matrix)
  sqrt(variance_hat / within_variance)
}

chain_slice <- function(chain_array, variable) {
  if (
    length(dim(chain_array)) != 3L ||
    is.null(dimnames(chain_array)[[1L]]) ||
    !variable %in% dimnames(chain_array)[[1L]]
  ) {
    stop("Unexpected MICE chain-array structure.")
  }
  result <- chain_array[variable, 11:20, , drop = TRUE]
  if (is.null(dim(result))) stop("Invalid MICE chain slice.")
  result
}

qc_rows <- list()
add_qc <- function(
    qc_domain,
    variable,
    metric,
    observed_value = NA_real_,
    imputed_value = NA_real_,
    diagnostic_value = NA_real_,
    threshold = NA_real_,
    status = "INFO",
    notes = "") {
  qc_rows[[length(qc_rows) + 1L]] <<- data.frame(
    qc_domain = qc_domain,
    variable = variable,
    metric = metric,
    observed_value = observed_value,
    imputed_value = imputed_value,
    diagnostic_value = diagnostic_value,
    threshold = threshold,
    status = status,
    notes = notes,
    stringsAsFactors = FALSE
  )
}

add_qc(
  "input", "all", "cohort_n", observed_value = nrow(analysis_data),
  diagnostic_value = nrow(analysis_data), threshold = expected_n,
  status = if (nrow(analysis_data) == expected_n) "PASS" else "FAIL",
  notes = "Solid tumor + sepsis screening cohort with ICU LOS >=24h."
)
for (variable in primary_variables) {
  add_qc(
    "input_missingness", variable, "missing_n",
    observed_value = original_missing_n[[variable]],
    diagnostic_value = original_missing_n[[variable]],
    notes = "Original first-24-hour median before eICU MICE."
  )
  add_qc(
    "input_missingness", variable, "missing_percent",
    observed_value = 100 * original_missing_n[[variable]] / expected_n,
    diagnostic_value = 100 * original_missing_n[[variable]] / expected_n,
    notes = "Original first-24-hour median before eICU MICE."
  )
}
for (variable in auxiliary_variables) {
  add_qc(
    "auxiliary_completeness", variable, "missing_n",
    observed_value = original_missing_n[[variable]],
    diagnostic_value = original_missing_n[[variable]], threshold = 0,
    status = if (original_missing_n[[variable]] == 0L) "PASS" else "FAIL",
    notes = "Predictor-only auxiliary; not imputed."
  )
}
for (variable in primary_variables) {
  mean_rhat <- gelman_rhat(chain_slice(imputation$chainMean, variable))
  variance_rhat <- gelman_rhat(chain_slice(imputation$chainVar, variable))
  add_qc(
    "convergence", variable, "chain_mean_rhat",
    diagnostic_value = mean_rhat, threshold = 1.10,
    status = if (is.finite(mean_rhat) && mean_rhat <= 1.10) "PASS" else "FAIL",
    notes = "Iterations 11-20 across ten chains."
  )
  add_qc(
    "convergence", variable, "chain_variance_rhat",
    diagnostic_value = variance_rhat, threshold = 1.10,
    status = if (is.finite(variance_rhat) && variance_rhat <= 1.10) "PASS" else "FAIL",
    notes = "Iterations 11-20 across ten chains."
  )
}

logged_events_n <- if (is.null(imputation$loggedEvents)) 0L else nrow(imputation$loggedEvents)
add_qc(
  "engine", "all", "logged_events_n",
  diagnostic_value = logged_events_n, threshold = 0,
  status = if (logged_events_n == 0L) "PASS" else "FAIL",
  notes = "MICE logged events must be absent."
)
add_qc("engine", "all", "m", diagnostic_value = imputation$m, threshold = 10,
       status = if (imputation$m == 10L) "PASS" else "FAIL")
add_qc("engine", "all", "maxit", diagnostic_value = imputation$iteration, threshold = 20,
       status = if (imputation$iteration == 20L) "PASS" else "FAIL")
add_qc("engine", "all", "donors", diagnostic_value = 5, threshold = 5, status = "PASS")
add_qc("engine", "all", "seed", diagnostic_value = stage_seed, threshold = stage_seed, status = "PASS")

completed_datasets <- mice::complete(imputation, action = "all")
if (length(completed_datasets) != 10L) stop("Expected ten completed eICU datasets.")
completed_missing <- vapply(
  completed_datasets,
  function(data) colSums(is.na(data)),
  numeric(length(expected_variables))
)
rownames(completed_missing) <- expected_variables
for (variable in expected_variables) {
  maximum_missing <- max(completed_missing[variable, ])
  total_missing <- sum(completed_missing[variable, ])
  add_qc(
    "missing_after_imputation", variable, "max_missing_n",
    observed_value = original_missing_n[[variable]],
    imputed_value = maximum_missing,
    diagnostic_value = total_missing, threshold = 0,
    status = if (maximum_missing == 0L && total_missing == 0L) "PASS" else "FAIL",
    notes = "Maximum and total residual missingness across ten completed datasets."
  )
}

distribution_probabilities <- c(p10 = 0.10, p25 = 0.25, median = 0.50, p75 = 0.75, p90 = 0.90)
for (variable in primary_variables) {
  observed_values <- analysis_data[[variable]][!is.na(analysis_data[[variable]])]
  imputed_values <- as.numeric(unlist(imputation$imp[[variable]], use.names = FALSE))
  add_qc(
    "distribution", variable, "n",
    observed_value = length(observed_values), imputed_value = length(imputed_values),
    notes = "Observed cells versus pooled imputed cells across ten datasets."
  )
  observed_quantiles <- stats::quantile(
    observed_values, probs = distribution_probabilities, names = FALSE, type = 7
  )
  imputed_quantiles <- stats::quantile(
    imputed_values, probs = distribution_probabilities, names = FALSE, type = 7
  )
  for (index in seq_along(distribution_probabilities)) {
    add_qc(
      "distribution", variable, names(distribution_probabilities)[index],
      observed_value = observed_quantiles[index],
      imputed_value = imputed_quantiles[index],
      notes = "Original measurement scale."
    )
  }
  add_qc(
    "distribution", variable, "IQR",
    observed_value = stats::IQR(observed_values, type = 7),
    imputed_value = stats::IQR(imputed_values, type = 7),
    notes = "Original measurement scale."
  )
  outside_support_n <- sum(
    imputed_values < min(observed_values) |
      imputed_values > max(observed_values)
  )
  add_qc(
    "plausibility", variable, "outside_observed_support_n",
    diagnostic_value = outside_support_n, threshold = 0,
    status = if (outside_support_n == 0L) "PASS" else "FAIL",
    notes = "PMM draws must remain within observed support."
  )
}

imputation_qc <- do.call(rbind, qc_rows)
row.names(imputation_qc) <- NULL
numeric_qc_columns <- c(
  "observed_value", "imputed_value", "diagnostic_value", "threshold"
)
imputation_qc[numeric_qc_columns] <- lapply(
  imputation_qc[numeric_qc_columns],
  function(x) ifelse(is.na(x), NA_real_, signif(x, 8L))
)
imputation_qc_pass <- !any(imputation_qc$status == "FAIL")

imputation_qc_tmp <- tempfile(
  "stage10c1_imputation_qc_", tmpdir = dirname(imputation_qc_path), fileext = ".csv"
)
utils::write.csv(imputation_qc, imputation_qc_tmp, row.names = FALSE, na = "")
if (!imputation_qc_pass) {
  if (!file.rename(imputation_qc_tmp, imputation_qc_path)) unlink(imputation_qc_tmp)
  stop("Stage 10C-1 MICE QC failed; no transfer assignments or controlled object were finalized.")
}

score_with_frozen_transfer <- function(data) {
  x <- as.data.frame(data[, primary_variables, drop = FALSE])
  x[] <- lapply(x, as.numeric)
  if (anyNA(x) || any(!is.finite(as.matrix(x)))) {
    stop("Transfer predictors must be complete and finite.")
  }
  if (any(x$lactate < 0) || any(x$creatinine < 0)) {
    stop("Frozen log1p transfer inputs must be nonnegative.")
  }
  x$lactate <- log1p(x$lactate)
  x$creatinine <- log1p(x$creatinine)
  x_matrix <- as.matrix(x)
  component_probability <- matrix(
    NA_real_, nrow = nrow(x_matrix), ncol = ncol(transfer_model$coefficients)
  )
  for (component in seq_len(ncol(transfer_model$coefficients))) {
    z <- sweep(x_matrix, 2L, transfer_model$centers[, component], "-")
    z <- sweep(z, 2L, transfer_model$scales[, component], "/")
    beta <- transfer_model$coefficients[, component]
    linear_predictor <- as.numeric(beta[1L] + z %*% beta[-1L])
    component_probability[, component] <- stats::plogis(linear_predictor)
  }
  probability <- rowMeans(component_probability)
  assignment <- ifelse(
    probability >= transfer_model$classification_threshold, 2L, 1L
  )
  if (
    anyNA(probability) || any(!is.finite(probability)) ||
    any(probability < 0 | probability > 1) ||
    anyNA(assignment) || !all(assignment %in% 1:2)
  ) {
    stop("Frozen transfer scoring produced an invalid probability or assignment.")
  }
  list(probability = probability, assignment = assignment)
}

profile_probability <- matrix(NA_real_, expected_n, 10L)
profile_assignment <- matrix(NA_integer_, expected_n, 10L)
for (imputation_index in seq_len(10L)) {
  scored <- score_with_frozen_transfer(completed_datasets[[imputation_index]])
  profile_probability[, imputation_index] <- scored$probability
  profile_assignment[, imputation_index] <- scored$assignment
}
dimnames(profile_probability) <- NULL
dimnames(profile_assignment) <- NULL

complete_case_index <- complete.cases(analysis_data[primary_variables])
complete_case_n <- sum(complete_case_index)
if (complete_case_n != 741L) {
  stop("The six-variable complete-case count does not reconcile to Stage 10C preflight.")
}
complete_case_scored <- score_with_frozen_transfer(
  analysis_data[complete_case_index, , drop = FALSE]
)

assignment_rows <- list()
add_assignment_qc <- function(
    section,
    analysis_set,
    imputation,
    metric,
    profile = "",
    value = NA_real_,
    denominator_n = NA_real_,
    percent = NA_real_,
    minimum = NA_real_,
    maximum = NA_real_,
    threshold = NA_real_,
    status = "INFO",
    notes = "") {
  assignment_rows[[length(assignment_rows) + 1L]] <<- data.frame(
    section = section,
    analysis_set = analysis_set,
    imputation = imputation,
    metric = metric,
    profile = profile,
    value = value,
    denominator_n = denominator_n,
    percent = percent,
    minimum = minimum,
    maximum = maximum,
    threshold = threshold,
    status = status,
    notes = notes,
    stringsAsFactors = FALSE
  )
}

add_assignment_qc(
  "input_safeguard", "all", "all", "target_cohort_n",
  value = expected_n, denominator_n = expected_n, percent = 100,
  threshold = expected_n, status = "PASS",
  notes = "Solid tumor + sepsis screening cohort with ICU LOS >=24h."
)
add_assignment_qc(
  "input_safeguard", "all", "all", "approved_input_variable_n",
  value = length(expected_variables), threshold = 9, status = "PASS",
  notes = "Six transfer variables plus age, sex, and cancer phenotype."
)
add_assignment_qc(
  "input_safeguard", "all", "all", "prohibited_input_field_n",
  value = 0, threshold = 0, status = "PASS",
  notes = "No identifier, mortality, outcome, discharge, or survival field entered R."
)
add_assignment_qc(
  "transfer_rule", "MICE", "all", "component_model_n",
  value = ncol(transfer_model$coefficients), threshold = 10, status = "PASS",
  notes = "All components are from the frozen MIMIC transfer object."
)
add_assignment_qc(
  "transfer_rule", "MICE", "all", "classification_threshold",
  value = transfer_model$classification_threshold,
  threshold = transfer_model$classification_threshold, status = "PASS",
  notes = "Probability >=0.5 assigns Profile 2; no eICU recalibration."
)
add_assignment_qc(
  "transfer_rule", "MICE", "all", "eicu_model_refit_n",
  value = 0, threshold = 0, status = "PASS",
  notes = "The MIMIC classifier was applied without retraining."
)
add_assignment_qc(
  "assignment_completeness", "MICE", "all", "complete_assignment_n",
  value = sum(!is.na(profile_assignment)),
  denominator_n = length(profile_assignment),
  percent = 100 * sum(!is.na(profile_assignment)) / length(profile_assignment),
  threshold = length(profile_assignment),
  status = if (!anyNA(profile_assignment)) "PASS" else "FAIL",
  notes = "All stay-by-imputation assignments must be complete."
)
add_assignment_qc(
  "assignment_completeness", "MICE", "all", "finite_probability_n",
  value = sum(is.finite(profile_probability)),
  denominator_n = length(profile_probability),
  percent = 100 * sum(is.finite(profile_probability)) / length(profile_probability),
  threshold = length(profile_probability),
  status = if (all(is.finite(profile_probability))) "PASS" else "FAIL",
  notes = "All stay-by-imputation transfer probabilities must be finite."
)

profile_counts <- matrix(0L, nrow = 2L, ncol = 10L)
for (imputation_index in seq_len(10L)) {
  profile_counts[, imputation_index] <- tabulate(
    profile_assignment[, imputation_index], nbins = 2L
  )
  for (profile_index in 1:2) {
    add_assignment_qc(
      "profile_distribution", "MICE", sprintf("%02d", imputation_index),
      "profile_n", paste("Profile", profile_index),
      value = profile_counts[profile_index, imputation_index],
      denominator_n = expected_n,
      percent = 100 * profile_counts[profile_index, imputation_index] / expected_n,
      status = "INFO",
      notes = "Imputation-specific application of the frozen ten-component ensemble."
    )
  }
}
profile_percent <- 100 * profile_counts / expected_n
for (profile_index in 1:2) {
  add_assignment_qc(
    "profile_distribution_summary", "MICE", "all", "profile_n_across_imputations",
    paste("Profile", profile_index),
    value = mean(profile_counts[profile_index, ]),
    denominator_n = expected_n,
    minimum = min(profile_counts[profile_index, ]),
    maximum = max(profile_counts[profile_index, ]),
    status = "INFO",
    notes = "Mean and range across ten eICU completed datasets."
  )
  add_assignment_qc(
    "profile_distribution_summary", "MICE", "all", "profile_percent_across_imputations",
    paste("Profile", profile_index),
    value = mean(profile_percent[profile_index, ]),
    percent = mean(profile_percent[profile_index, ]),
    minimum = min(profile_percent[profile_index, ]),
    maximum = max(profile_percent[profile_index, ]),
    status = "INFO",
    notes = "Mean percentage and range across ten eICU completed datasets."
  )
}

unanimous <- apply(profile_assignment, 1L, function(x) length(unique(x)) == 1L)
add_assignment_qc(
  "cross_imputation_stability", "MICE", "all", "unanimous_assignment_n",
  value = sum(unanimous), denominator_n = expected_n,
  percent = 100 * mean(unanimous), status = "INFO",
  notes = "Same transferred profile in all ten eICU completed datasets."
)
add_assignment_qc(
  "cross_imputation_stability", "MICE", "all", "variable_assignment_n",
  value = sum(!unanimous), denominator_n = expected_n,
  percent = 100 * mean(!unanimous), status = "INFO",
  notes = "Transferred profile differs across at least one imputation."
)

votes_profile_1 <- rowSums(profile_assignment == 1L)
votes_profile_2 <- rowSums(profile_assignment == 2L)
consensus_profile <- ifelse(
  votes_profile_1 > votes_profile_2, 1L,
  ifelse(
    votes_profile_2 > votes_profile_1, 2L,
    ifelse(rowMeans(profile_probability) >= 0.5, 2L, 1L)
  )
)
consensus_counts <- tabulate(consensus_profile, nbins = 2L)
for (profile_index in 1:2) {
  add_assignment_qc(
    "profile_distribution", "QC_consensus", "all", "profile_n",
    paste("Profile", profile_index), value = consensus_counts[profile_index],
    denominator_n = expected_n,
    percent = 100 * consensus_counts[profile_index] / expected_n,
    status = "INFO",
    notes = "Modal assignment for QC only; future primary analysis remains imputation-specific."
  )
}

complete_case_counts <- tabulate(complete_case_scored$assignment, nbins = 2L)
for (profile_index in 1:2) {
  add_assignment_qc(
    "profile_distribution", "complete_case", "original", "profile_n",
    paste("Profile", profile_index), value = complete_case_counts[profile_index],
    denominator_n = complete_case_n,
    percent = 100 * complete_case_counts[profile_index] / complete_case_n,
    status = "INFO",
    notes = "Observed values only; same frozen classifier and threshold."
  )
}

assignment_qc <- do.call(rbind, assignment_rows)
row.names(assignment_qc) <- NULL
numeric_assignment_columns <- c(
  "value", "denominator_n", "percent", "minimum", "maximum", "threshold"
)
assignment_qc[numeric_assignment_columns] <- lapply(
  assignment_qc[numeric_assignment_columns],
  function(x) ifelse(is.na(x), NA_real_, signif(x, 8L))
)
assignment_qc_pass <- !any(assignment_qc$status == "FAIL")
if (!assignment_qc_pass) {
  unlink(imputation_qc_tmp)
  stop("Stage 10C-1 transfer-assignment QC failed; no output was finalized.")
}

stage10c1_metadata <- list(
  stage = "Stage 10C-1",
  created_date = as.character(Sys.Date()),
  source_database = "eICU-CRD",
  cohort_definition = "Adult solid tumor + diagnostic sepsis screening proxy + ICU LOS >=24h",
  row_n = expected_n,
  variables = expected_variables,
  primary_variables = primary_variables,
  auxiliary_variables = auxiliary_variables,
  row_order_contract = "Deterministic source query order; patient identifier not selected or stored",
  m = imputation$m,
  maxit = imputation$iteration,
  donors = 5L,
  seed = stage_seed,
  method = method,
  predictor_matrix = predictor_matrix,
  transfer_model_path = transfer_path,
  transfer_model_md5 = unname(tools::md5sum(transfer_path)),
  input_sql_md5 = unname(tools::md5sum(sql_path)),
  transfer_model_retrained = FALSE,
  eicu_recalibration_performed = FALSE,
  patient_identifiers_stored = FALSE,
  mortality_or_outcome_fields_read_or_stored = FALSE,
  survival_analysis_run = FALSE,
  patient_level_csv_exported = FALSE,
  imputation_qc_pass = imputation_qc_pass,
  assignment_qc_pass = assignment_qc_pass,
  r_version = R.version.string,
  mice_version = as.character(utils::packageVersion("mice"))
)
attr(imputation, "stage10c1_metadata") <- stage10c1_metadata
attr(imputation, "stage10c1_imputation_qc") <- imputation_qc

controlled_object <- list(
  stage = "Stage 10C-1",
  imputation = imputation,
  profile_probability = profile_probability,
  profile_assignment = profile_assignment,
  metadata = stage10c1_metadata
)

object_tmp <- tempfile("stage10c1_object_", tmpdir = dirname(object_path), fileext = ".rds")
assignment_qc_tmp <- tempfile(
  "stage10c1_assignment_qc_", tmpdir = dirname(assignment_qc_path), fileext = ".csv"
)
saveRDS(controlled_object, object_tmp, compress = "gzip", version = 3L)
utils::write.csv(assignment_qc, assignment_qc_tmp, row.names = FALSE, na = "")

verification_object <- readRDS(object_tmp)
verification_assignment_qc <- utils::read.csv(
  assignment_qc_tmp, check.names = FALSE, stringsAsFactors = FALSE
)
verification_imputation_qc <- utils::read.csv(
  imputation_qc_tmp, check.names = FALSE, stringsAsFactors = FALSE
)
if (
  !identical(names(verification_object), c(
    "stage", "imputation", "profile_probability", "profile_assignment", "metadata"
  )) ||
  !inherits(verification_object$imputation, "mids") ||
  verification_object$imputation$m != 10L ||
  verification_object$imputation$iteration != 20L ||
  nrow(verification_object$imputation$data) != expected_n ||
  !identical(names(verification_object$imputation$data), expected_variables) ||
  contains_prohibited_name(names(verification_object$imputation$data)) ||
  !identical(dim(verification_object$profile_probability), c(expected_n, 10L)) ||
  !identical(dim(verification_object$profile_assignment), c(expected_n, 10L)) ||
  anyNA(verification_object$profile_assignment) ||
  !all(verification_object$profile_assignment %in% 1:2) ||
  nrow(verification_assignment_qc) != nrow(assignment_qc) ||
  nrow(verification_imputation_qc) != nrow(imputation_qc) ||
  any(grepl("patientunitstayid|subject_id|stay_id", names(verification_assignment_qc), ignore.case = TRUE)) ||
  any(grepl("patientunitstayid|subject_id|stay_id", names(verification_imputation_qc), ignore.case = TRUE)) ||
  !identical(verification_object$metadata$patient_identifiers_stored, FALSE) ||
  !identical(verification_object$metadata$mortality_or_outcome_fields_read_or_stored, FALSE) ||
  !identical(verification_object$metadata$survival_analysis_run, FALSE)
) {
  unlink(c(object_tmp, assignment_qc_tmp, imputation_qc_tmp))
  stop("Stage 10C-1 temporary outputs failed privacy or structural verification.")
}

if (!file.rename(object_tmp, object_path)) {
  unlink(c(object_tmp, assignment_qc_tmp, imputation_qc_tmp))
  stop("Could not finalize the protected Stage 10C-1 object.")
}
if (!file.rename(assignment_qc_tmp, assignment_qc_path)) {
  unlink(c(assignment_qc_tmp, imputation_qc_tmp))
  stop("Protected object finalized but assignment QC could not be finalized.")
}
if (!file.rename(imputation_qc_tmp, imputation_qc_path)) {
  unlink(imputation_qc_tmp)
  stop("Protected object and assignment QC finalized but imputation QC could not be finalized.")
}

cat("STAGE10C1_COMPLETE\n")
cat(sprintf(
  "rows=%d complete_cases=%d m=%d maxit=%d donors=%d seed=%d\n",
  expected_n, complete_case_n, imputation$m, imputation$iteration, 5L, stage_seed
))
cat(sprintf(
  "profile1_mean_percent=%.4f profile2_mean_percent=%.4f unanimous_percent=%.4f\n",
  mean(profile_percent[1L, ]), mean(profile_percent[2L, ]), 100 * mean(unanimous)
))
cat(sprintf(
  "imputation_qc_rows=%d assignment_qc_rows=%d identifiers_stored=FALSE outcomes_used=FALSE survival_model_run=FALSE\n",
  nrow(imputation_qc), nrow(assignment_qc)
))
