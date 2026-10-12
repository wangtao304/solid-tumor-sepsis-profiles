# Stage 9A: aggregate, outcome-blind clinical description of the selected profiles.

options(stringsAsFactors = FALSE)

mice_path <- file.path(
  "data",
  "processed",
  "inflammatory_metabolic_mice_solid_m10_maxit20_seed20260721.rds"
)
lpa_path <- file.path(
  "data",
  "processed",
  "inflammatory_metabolic_lpa_solid_2to6_mi10_seed20260722.rds"
)
profile_characteristics_path <- file.path(
  "outputs",
  "tables",
  "28_lpa_profile_characteristics.csv"
)
output_path <- file.path(
  "outputs",
  "tables",
  "29_profile_baseline_characteristics.csv"
)

required_mice_variables <- c(
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

required_paths <- c(mice_path, lpa_path, profile_characteristics_path)
if (!all(file.exists(required_paths))) {
  stop(
    "Missing required Stage 9A local input(s): ",
    paste(required_paths[!file.exists(required_paths)], collapse = ", ")
  )
}
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

mice_object <- readRDS(mice_path)
lpa_object <- readRDS(lpa_path)
selected_characteristics <- utils::read.csv(
  profile_characteristics_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (
  !inherits(mice_object, "mids") ||
  mice_object$m != 10L ||
  nrow(mice_object$data) != 4380L ||
  !setequal(names(mice_object$data), required_mice_variables)
) {
  stop("Approved MICE object does not match the frozen Stage 9A input schema.")
}
if (
  !identical(lpa_object$stage, "Stage 8A-1") ||
  !identical(lpa_object$configuration$candidate_profiles, 2:6) ||
  !identical(lpa_object$outcome_variables_used, FALSE) ||
  !identical(lpa_object$survival_information_used, FALSE) ||
  !identical(lpa_object$patient_level_profile_labels_exported, FALSE)
) {
  stop("Stage 8A-1 LPA object failed Stage 9A safeguard checks.")
}
if (
  !setequal(unique(selected_characteristics$recommended_profiles), 2L) ||
  !setequal(unique(selected_characteristics$profile), 1:2)
) {
  stop("Stage 8A-2 selected profile metadata is not the expected 2-profile solution.")
}

profile_labels <- vapply(1:2, function(profile) {
  labels <- unique(selected_characteristics$profile_label[
    selected_characteristics$profile == profile
  ])
  if (length(labels) != 1L) stop("Profile label metadata is inconsistent.")
  labels
}, character(1L))

completed_datasets <- mice::complete(mice_object, action = "all")
if (
  length(completed_datasets) != 10L ||
  any(vapply(completed_datasets, function(data) anyNA(data), logical(1L)))
) {
  stop("Approved MICE completed datasets failed residual-missingness checks.")
}

profile_assignments <- matrix(
  NA_integer_,
  nrow = nrow(mice_object$data),
  ncol = length(completed_datasets)
)
mean_aligned_posterior <- matrix(
  0,
  nrow = nrow(mice_object$data),
  ncol = 2L
)
for (imputation in seq_along(completed_datasets)) {
  key <- sprintf("imp%02d_k2", imputation)
  fit_record <- lpa_object$fit_records[[key]]
  alignment <- lpa_object$alignments[[key]]
  if (
    is.null(fit_record) ||
    is.null(alignment) ||
    length(alignment) != 2L ||
    length(fit_record$model$classification) != nrow(mice_object$data) ||
    !identical(dim(fit_record$model$posterior), c(nrow(mice_object$data), 2L))
  ) {
    stop("Missing or malformed aligned 2-profile fit for ", key, ".")
  }
  profile_assignments[, imputation] <- match(
    fit_record$model$classification,
    alignment
  )
  mean_aligned_posterior <- mean_aligned_posterior +
    fit_record$model$posterior[, alignment, drop = FALSE] /
    length(completed_datasets)
}
if (anyNA(profile_assignments)) {
  stop("Aligned profile assignment reconstruction produced missing values.")
}

profile_1_votes <- rowSums(profile_assignments == 1L)
profile_2_votes <- rowSums(profile_assignments == 2L)
consensus_profile <- ifelse(
  profile_1_votes > profile_2_votes,
  1L,
  ifelse(
    profile_2_votes > profile_1_votes,
    2L,
    max.col(mean_aligned_posterior, ties.method = "first")
  )
)
consensus_sizes <- tabulate(consensus_profile, nbins = 2L)
if (
  sum(consensus_sizes) != nrow(mice_object$data) ||
  any(consensus_sizes / sum(consensus_sizes) < 0.05)
) {
  stop("Consensus profile assignment failed aggregate size checks.")
}

format_quantiles <- function(values, digits) {
  sprintf(
    paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"),
    values[2L],
    values[1L],
    values[3L]
  )
}

pooled_quantiles <- function(variable, profile = NULL) {
  estimates <- vapply(completed_datasets, function(data) {
    values <- data[[variable]]
    if (!is.null(profile)) values <- values[consensus_profile == profile]
    stats::quantile(
      values,
      probs = c(0.25, 0.50, 0.75),
      names = FALSE,
      na.rm = FALSE,
      type = 7
    )
  }, numeric(3L))
  rowMeans(estimates)
}

continuous_row <- function(section, label, variable, unit, digits) {
  data.frame(
    section = section,
    variable = label,
    level = "",
    unit = unit,
    statistic = "Median (IQR)",
    overall = format_quantiles(pooled_quantiles(variable), digits),
    profile_1_lower_stress = format_quantiles(
      pooled_quantiles(variable, 1L),
      digits
    ),
    profile_2_higher_stress = format_quantiles(
      pooled_quantiles(variable, 2L),
      digits
    ),
    source_note = paste(
      "Average of the 10 imputation-specific medians and quartiles;",
      "descriptive summary, not Rubin-rule inference."
    )
  )
}

format_count_percent <- function(count, denominator) {
  sprintf("%d (%.1f%%)", count, 100 * count / denominator)
}

categorical_row <- function(section, label, level, indicator, note = NULL) {
  overall_n <- sum(indicator)
  profile_1_n <- sum(indicator & consensus_profile == 1L)
  profile_2_n <- sum(indicator & consensus_profile == 2L)
  if (is.null(note)) {
    note <- "Counts use the outcome-blind consensus profile assignment."
  }
  data.frame(
    section = section,
    variable = label,
    level = level,
    unit = "",
    statistic = "n (%)",
    overall = format_count_percent(overall_n, length(indicator)),
    profile_1_lower_stress = format_count_percent(
      profile_1_n,
      consensus_sizes[1L]
    ),
    profile_2_higher_stress = format_count_percent(
      profile_2_n,
      consensus_sizes[2L]
    ),
    source_note = note
  )
}

first_completed <- completed_datasets[[1L]]
invariant_variables <- c(
  "age",
  "sex",
  "cancer_phenotype",
  "charlson_comorbidity_index",
  "sofa",
  "sapsii",
  "vasopressor_0_24",
  "mechanical_ventilation_0_24"
)
for (variable in invariant_variables) {
  invariant <- vapply(completed_datasets[-1L], function(data) {
    identical(data[[variable]], first_completed[[variable]])
  }, logical(1L))
  if (!all(invariant)) {
    stop("Expected invariant descriptor differs across imputations: ", variable)
  }
}

as_yes <- function(values) {
  if (is.logical(values)) return(values)
  as.character(values) %in% c("1", "TRUE", "Yes", "yes")
}

rows <- list()
rows[[length(rows) + 1L]] <- data.frame(
  section = "Cohort",
  variable = "ICU stays",
  level = "",
  unit = "",
  statistic = "n",
  overall = as.character(nrow(mice_object$data)),
  profile_1_lower_stress = as.character(consensus_sizes[1L]),
  profile_2_higher_stress = as.character(consensus_sizes[2L]),
  source_note = sprintf(
    paste0(
      "Modal vote across 10 aligned assignments; %d ties were resolved by ",
      "higher mean posterior probability. Profile labels: P1=%s; P2=%s."
    ),
    sum(profile_1_votes == profile_2_votes),
    profile_labels[1L],
    profile_labels[2L]
  )
)

rows[[length(rows) + 1L]] <- continuous_row(
  "Demographics", "Age", "age", "years", 1L
)
rows[[length(rows) + 1L]] <- categorical_row(
  "Demographics", "Sex", "Female", as.character(first_completed$sex) == "F"
)
rows[[length(rows) + 1L]] <- categorical_row(
  "Demographics", "Sex", "Male", as.character(first_completed$sex) == "M"
)

metastatic <- as.character(first_completed$cancer_phenotype) == "metastatic solid"
rows[[length(rows) + 1L]] <- categorical_row(
  "Cancer", "Metastatic status", "Metastatic", metastatic
)
rows[[length(rows) + 1L]] <- categorical_row(
  "Cancer", "Metastatic status", "Nonmetastatic", !metastatic
)
rows[[length(rows) + 1L]] <- categorical_row(
  "Cancer",
  "Cancer phenotype",
  "Metastatic solid",
  metastatic
)
rows[[length(rows) + 1L]] <- categorical_row(
  "Cancer",
  "Cancer phenotype",
  "Nonmetastatic solid",
  !metastatic
)

rows[[length(rows) + 1L]] <- continuous_row(
  "Severity", "SOFA", "sofa", "points", 1L
)
rows[[length(rows) + 1L]] <- continuous_row(
  "Severity", "SAPS-II", "sapsii", "points", 1L
)
rows[[length(rows) + 1L]] <- continuous_row(
  "Severity",
  "Charlson Comorbidity Index",
  "charlson_comorbidity_index",
  "points",
  1L
)

rows[[length(rows) + 1L]] <- categorical_row(
  "Treatment intensity",
  "Vasopressor use within first 24h",
  "Yes",
  as_yes(first_completed$vasopressor_0_24),
  "Yes level shown; no is the complement."
)
rows[[length(rows) + 1L]] <- categorical_row(
  "Treatment intensity",
  "Invasive mechanical ventilation within first 24h",
  "Yes",
  as_yes(first_completed$mechanical_ventilation_0_24),
  "Yes level shown; no is the complement."
)
rows[[length(rows) + 1L]] <- data.frame(
  section = "Treatment intensity",
  variable = "CRRT",
  level = "",
  unit = "",
  statistic = "Not available",
  overall = "Not available",
  profile_1_lower_stress = "Not available",
  profile_2_higher_stress = "Not available",
  source_note = paste(
    "CRRT is absent from the approved local MICE/LPA objects;",
    "no database reconnection was performed in Stage 9A."
  )
)

rows[[length(rows) + 1L]] <- continuous_row(
  "Laboratory", "WBC", "wbc", "K/uL", 1L
)
rows[[length(rows) + 1L]] <- continuous_row(
  "Laboratory", "Platelet count", "platelet_count", "K/uL", 0L
)
rows[[length(rows) + 1L]] <- continuous_row(
  "Laboratory", "Lactate", "lactate", "mmol/L", 1L
)
rows[[length(rows) + 1L]] <- continuous_row(
  "Laboratory", "Glucose", "glucose", "mg/dL", 0L
)
rows[[length(rows) + 1L]] <- continuous_row(
  "Laboratory", "Creatinine", "creatinine", "mg/dL", 1L
)
rows[[length(rows) + 1L]] <- continuous_row(
  "Laboratory", "Bicarbonate", "bicarbonate", "mEq/L", 1L
)

table_1 <- do.call(rbind, rows)
row.names(table_1) <- NULL

temporary_output <- tempfile(
  "stage9a_table1_",
  tmpdir = dirname(output_path),
  fileext = ".csv"
)
utils::write.csv(table_1, temporary_output, row.names = FALSE, na = "")
verification <- utils::read.csv(
  temporary_output,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

prohibited_output_fields <- c(
  "subject_id", "stay_id", "profile_assignment", "posterior_probability"
)
if (
  nrow(verification) != 20L ||
  ncol(verification) != 9L ||
  any(names(verification) %in% prohibited_output_fields) ||
  sum(verification$variable == "CRRT") != 1L ||
  verification$statistic[verification$variable == "CRRT"] != "Not available" ||
  verification$overall[verification$variable == "ICU stays"] != "4380" ||
  anyNA(verification)
) {
  stop("Temporary Stage 9A aggregate Table 1 failed structural verification.")
}
if (!file.copy(temporary_output, output_path, overwrite = TRUE)) {
  stop("Could not finalize Stage 9A aggregate Table 1.")
}
unlink(temporary_output)

cat("STAGE9A_PROFILE_DESCRIPTION_COMPLETED\n")
cat(sprintf(
  paste0(
    "profile_1_n=%d profile_2_n=%d unanimous_n=%d tie_n=%d rows=%d ",
    "crrt_available=FALSE no_patient_level_export=TRUE\n"
  ),
  consensus_sizes[1L],
  consensus_sizes[2L],
  sum(pmax(profile_1_votes, profile_2_votes) == 10L),
  sum(profile_1_votes == profile_2_votes),
  nrow(table_1)
))
