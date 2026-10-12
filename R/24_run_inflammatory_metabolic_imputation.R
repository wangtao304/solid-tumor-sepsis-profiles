#!/usr/bin/env Rscript

# Stage 7B: execute the frozen solid-tumor inflammatory-metabolic MICE plan.
# Patient-level values are streamed directly from PostgreSQL into R memory.
# Only a protected local mids object and aggregate QC are finalized.

options(stringsAsFactors = FALSE, warn = 1)

if (!requireNamespace("mice", quietly = TRUE)) {
  stop("Required package 'mice' is not installed; no installation is authorized.")
}

expected_variables <- c(
  "wbc",
  "platelet_count",
  "lactate",
  "glucose",
  "creatinine",
  "bicarbonate",
  "age",
  "sex",
  "cancer_phenotype",
  "charlson_comorbidity_index",
  "sofa",
  "sapsii",
  "vasopressor_0_24",
  "mechanical_ventilation_0_24"
)

primary_variables <- expected_variables[1:6]
auxiliary_variables <- expected_variables[7:14]
continuous_variables <- c(
  primary_variables,
  "age",
  "charlson_comorbidity_index",
  "sofa",
  "sapsii"
)

prohibited_variables <- c(
  "subject_id",
  "stay_id",
  "hadm_id",
  "mortality_28d",
  "icu_mortality",
  "hospital_mortality",
  "los_icu",
  "icu_los",
  "hospital_los",
  "los_hospital",
  "neutrophil",
  "neutrophil_abs",
  "lymphocyte",
  "lymphocyte_abs",
  "nlr",
  "anc",
  "alc",
  "anc_reconstructed",
  "alc_reconstructed"
)

psql_bin <- Sys.getenv(
  "MIMIC_PSQL_BIN",
  unset = ""
)
pgpass_file <- Sys.getenv("MIMIC_PGPASS_FILE", unset = "")
db_host <- Sys.getenv("MIMIC_DB_HOST", unset = "")
db_port <- Sys.getenv("MIMIC_DB_PORT", unset = "")
db_name <- Sys.getenv("MIMIC_DB_NAME", unset = "")
db_user <- Sys.getenv("MIMIC_DB_USER", unset = "")

if (!nzchar(psql_bin) || !file.exists(psql_bin) ||
    !nzchar(pgpass_file) || !file.exists(pgpass_file) ||
    any(!nzchar(c(db_host, db_port, db_name, db_user)))) {
  stop("Set all MIMIC database environment variables to authorized local values.")
}
Sys.setenv(PGPASSFILE = pgpass_file)

psql_args <- c(
  "-h", db_host,
  "-p", db_port,
  "-U", db_user,
  "-d", db_name,
  "-X",
  "-w",
  "-v", "ON_ERROR_STOP=1",
  "-q",
  "--csv",
  "-f", "sql/24_im_new_mice_input.sql"
)

csv_lines <- system2(
  psql_bin,
  args = psql_args,
  stdout = TRUE,
  stderr = ""
)
psql_status <- attr(csv_lines, "status")
if (!is.null(psql_status) && psql_status != 0L) {
  stop("psql failed while streaming the Stage 7B input.")
}
if (length(csv_lines) < 2L) {
  stop("psql returned no Stage 7B data rows.")
}

analysis_data <- utils::read.csv(
  text = paste(csv_lines, collapse = "\n"),
  header = TRUE,
  na.strings = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
rm(csv_lines)

if (!identical(names(analysis_data), expected_variables)) {
  stop("Input columns do not exactly match the frozen 14-variable order.")
}
if (any(prohibited_variables %in% names(analysis_data))) {
  stop("A prohibited variable entered the Stage 7B MICE input.")
}
if (nrow(analysis_data) != 4380L) {
  stop("Expected exactly 4380 solid-tumor rows in the Stage 7B input.")
}

analysis_data[continuous_variables] <- lapply(
  analysis_data[continuous_variables],
  as.numeric
)

analysis_data$sex <- factor(analysis_data$sex, levels = c("F", "M"))
analysis_data$cancer_phenotype <- factor(
  analysis_data$cancer_phenotype,
  levels = c("metastatic solid", "nonmetastatic solid")
)

as_binary_factor <- function(x) {
  normalized <- tolower(as.character(x))
  if (!all(normalized %in% c("t", "f", "true", "false", "1", "0"))) {
    stop("Unexpected coding in a binary auxiliary variable.")
  }
  factor(
    ifelse(normalized %in% c("t", "true", "1"), "Yes", "No"),
    levels = c("No", "Yes")
  )
}

analysis_data$vasopressor_0_24 <- as_binary_factor(
  analysis_data$vasopressor_0_24
)
analysis_data$mechanical_ventilation_0_24 <- as_binary_factor(
  analysis_data$mechanical_ventilation_0_24
)

if (anyNA(analysis_data$sex) || anyNA(analysis_data$cancer_phenotype)) {
  stop("Unexpected categorical coding introduced missing values.")
}
if (anyNA(analysis_data[auxiliary_variables])) {
  stop("A frozen Stage 7B auxiliary variable is unexpectedly missing.")
}

original_missing_n <- colSums(is.na(analysis_data))
if (!all(original_missing_n[primary_variables] > 0L)) {
  stop("Every frozen primary target was expected to have at least one missing value.")
}

method <- setNames(rep("", length(expected_variables)), expected_variables)
method[primary_variables] <- "pmm"

predictor_matrix <- matrix(
  0L,
  nrow = length(expected_variables),
  ncol = length(expected_variables),
  dimnames = list(expected_variables, expected_variables)
)
predictor_matrix[primary_variables, ] <- 1L
diag(predictor_matrix) <- 0L

if (!all(method[primary_variables] == "pmm") ||
    any(nzchar(method[auxiliary_variables])) ||
    any(diag(predictor_matrix) != 0L) ||
    any(predictor_matrix[auxiliary_variables, , drop = FALSE] != 0L)) {
  stop("Frozen Stage 7B method vector or predictor matrix is invalid.")
}

imputation <- mice::mice(
  data = analysis_data,
  m = 10L,
  method = method,
  predictorMatrix = predictor_matrix,
  maxit = 20L,
  printFlag = FALSE,
  seed = 20260721L,
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
  if (!is.finite(within_variance) || within_variance < 0) {
    return(NA_real_)
  }
  if (within_variance == 0) {
    return(if (between_variance == 0) 1 else Inf)
  }
  variance_hat <- (
    (nrow(chain_matrix) - 1) / nrow(chain_matrix)
  ) * within_variance + between_variance / nrow(chain_matrix)
  sqrt(variance_hat / within_variance)
}

chain_slice <- function(chain_array, variable) {
  if (length(dim(chain_array)) != 3L ||
      is.null(dimnames(chain_array)[[1L]]) ||
      !variable %in% dimnames(chain_array)[[1L]]) {
    stop("Unexpected MICE chain-array structure.")
  }
  result <- chain_array[variable, 11:20, , drop = TRUE]
  if (is.null(dim(result))) {
    stop("Chain slice did not produce an iteration-by-chain matrix.")
  }
  result
}

qc_rows <- list()
add_qc_row <- function(
    qc_domain,
    variable,
    metric,
    observed_value = NA_real_,
    imputed_value = NA_real_,
    diagnostic_value = NA_real_,
    threshold = NA_real_,
    pass = NA) {
  qc_rows[[length(qc_rows) + 1L]] <<- data.frame(
    qc_domain = qc_domain,
    variable = variable,
    metric = metric,
    observed_value = observed_value,
    imputed_value = imputed_value,
    diagnostic_value = diagnostic_value,
    threshold = threshold,
    pass = if (is.na(pass)) "" else if (isTRUE(pass)) "TRUE" else "FALSE",
    stringsAsFactors = FALSE
  )
}

for (variable in primary_variables) {
  mean_rhat <- gelman_rhat(chain_slice(imputation$chainMean, variable))
  variance_rhat <- gelman_rhat(chain_slice(imputation$chainVar, variable))
  add_qc_row(
    "convergence", variable, "chain_mean_rhat",
    diagnostic_value = round(mean_rhat, 6), threshold = 1.10,
    pass = is.finite(mean_rhat) && mean_rhat <= 1.10
  )
  add_qc_row(
    "convergence", variable, "chain_variance_rhat",
    diagnostic_value = round(variance_rhat, 6), threshold = 1.10,
    pass = is.finite(variance_rhat) && variance_rhat <= 1.10
  )
}

logged_events_n <- if (is.null(imputation$loggedEvents)) {
  0L
} else {
  nrow(imputation$loggedEvents)
}
add_qc_row(
  "engine", "all", "logged_events_n",
  diagnostic_value = logged_events_n, threshold = 0,
  pass = logged_events_n == 0L
)

completed_missing <- vapply(
  seq_len(imputation$m),
  function(index) colSums(is.na(mice::complete(imputation, action = index))),
  numeric(length(expected_variables))
)
rownames(completed_missing) <- expected_variables

for (variable in expected_variables) {
  maximum_missing <- max(completed_missing[variable, ])
  total_missing <- sum(completed_missing[variable, ])
  add_qc_row(
    "missing_after_imputation", variable, "max_missing_n",
    observed_value = original_missing_n[[variable]],
    imputed_value = maximum_missing,
    diagnostic_value = total_missing,
    threshold = 0,
    pass = maximum_missing == 0 && total_missing == 0
  )
}

distribution_probabilities <- c(
  p10 = 0.10,
  p25 = 0.25,
  median = 0.50,
  p75 = 0.75,
  p90 = 0.90
)

for (variable in primary_variables) {
  observed_values <- analysis_data[[variable]][!is.na(analysis_data[[variable]])]
  imputed_values <- as.numeric(unlist(imputation$imp[[variable]], use.names = FALSE))

  add_qc_row(
    "distribution", variable, "n",
    observed_value = length(observed_values),
    imputed_value = length(imputed_values)
  )

  observed_quantiles <- stats::quantile(
    observed_values,
    probs = distribution_probabilities,
    names = FALSE,
    type = 7
  )
  imputed_quantiles <- stats::quantile(
    imputed_values,
    probs = distribution_probabilities,
    names = FALSE,
    type = 7
  )

  for (index in seq_along(distribution_probabilities)) {
    add_qc_row(
      "distribution",
      variable,
      names(distribution_probabilities)[index],
      observed_value = round(observed_quantiles[index], 6),
      imputed_value = round(imputed_quantiles[index], 6)
    )
  }

  observed_iqr <- stats::IQR(observed_values, type = 7)
  imputed_iqr <- stats::IQR(imputed_values, type = 7)
  add_qc_row(
    "distribution", variable, "IQR",
    observed_value = round(observed_iqr, 6),
    imputed_value = round(imputed_iqr, 6)
  )

  outside_observed_support_n <- sum(
    imputed_values < min(observed_values) |
      imputed_values > max(observed_values)
  )
  add_qc_row(
    "plausibility", variable, "outside_observed_support_n",
    diagnostic_value = outside_observed_support_n,
    threshold = 0,
    pass = outside_observed_support_n == 0L
  )
}

qc_summary <- do.call(rbind, qc_rows)
row.names(qc_summary) <- NULL
qc_pass_values <- qc_summary$pass[nzchar(qc_summary$pass)]
qc_all_pass <- length(qc_pass_values) > 0L && all(qc_pass_values == "TRUE")

stage7b_metadata <- list(
  stage = "Stage 7B",
  source_table = "mimic_cancer_sepsis.immune_metabolic_lca_base",
  cohort_filter = c("metastatic solid", "nonmetastatic solid"),
  row_n = nrow(analysis_data),
  variables = expected_variables,
  primary_variables = primary_variables,
  auxiliary_variables = auxiliary_variables,
  prohibited_variables = prohibited_variables,
  m = imputation$m,
  maxit = imputation$iteration,
  donors = 5L,
  seed = 20260721L,
  method = method,
  predictor_matrix = predictor_matrix,
  r_version = R.version.string,
  mice_version = as.character(utils::packageVersion("mice")),
  qc_all_pass = qc_all_pass,
  patient_level_csv_exported = FALSE,
  lpa_run = FALSE,
  survival_analysis_run = FALSE
)

attr(imputation, "stage7b_metadata") <- stage7b_metadata
attr(imputation, "stage7b_qc_summary") <- qc_summary

object_path <- paste0(
  "data/processed/",
  "inflammatory_metabolic_mice_solid_m10_maxit20_seed20260721.rds"
)
qc_path <- "outputs/tables/24_im_new_imputation_qc.csv"

if (file.exists(object_path) || file.exists(qc_path)) {
  stop("A Stage 7B output already exists; refusing to overwrite it.")
}
if (!dir.exists(dirname(object_path)) || !dir.exists(dirname(qc_path))) {
  stop("A required controlled output directory does not exist.")
}

object_tmp <- tempfile("stage7b_mice_", tmpdir = dirname(object_path), fileext = ".rds")
qc_tmp <- tempfile("stage7b_qc_", tmpdir = dirname(qc_path), fileext = ".csv")
saveRDS(imputation, object_tmp, compress = "gzip")
utils::write.csv(qc_summary, qc_tmp, row.names = FALSE, na = "")

verification_object <- readRDS(object_tmp)
if (!inherits(verification_object, "mids") ||
    verification_object$m != 10L ||
    verification_object$iteration != 20L ||
    nrow(verification_object$data) != 4380L ||
    !identical(names(verification_object$data), expected_variables) ||
    any(prohibited_variables %in% names(verification_object$data))) {
  unlink(object_tmp)
  unlink(qc_tmp)
  stop("Temporary Stage 7B MICE object failed structural verification.")
}

if (!qc_all_pass) {
  unlink(object_tmp)
  if (!file.rename(qc_tmp, qc_path)) {
    unlink(qc_tmp)
  }
  stop("Stage 7B aggregate QC failed; the MICE object was not finalized.")
}

if (!file.rename(object_tmp, object_path)) {
  unlink(object_tmp)
  unlink(qc_tmp)
  stop("Could not atomically finalize the local Stage 7B MICE object.")
}
if (!file.rename(qc_tmp, qc_path)) {
  stop("MICE object finalized but aggregate QC CSV could not be finalized.")
}

final_object <- readRDS(object_path)
if (!inherits(final_object, "mids") ||
    !isTRUE(attr(final_object, "stage7b_metadata")$qc_all_pass) ||
    any(prohibited_variables %in% names(final_object$data))) {
  stop("Final Stage 7B MICE object failed structural verification.")
}

cat("STAGE7B_MICE_COMPLETED\n")
cat("rows=4380 m=10 maxit=20 seed=20260721 donors=5\n")
cat("qc_rows=", nrow(qc_summary), " qc_all_pass=", qc_all_pass, "\n", sep = "")
cat("logged_events_n=", logged_events_n, "\n", sep = "")
