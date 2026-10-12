options(stringsAsFactors = FALSE)

input_path <- "data/processed/immune_metabolic_mice_m10_maxit20_seed20260719.rds"
phenotype_path <- "outputs/tables/03_cancer_phenotype_distribution.csv"
eicu_path <- "outputs/tables/04_eicu_validation_feasibility.csv"
output_path <- "outputs/tables/19_solid_tumor_feasibility.csv"

if (file.exists(output_path)) {
  stop("Stage 6C output already exists; refusing to overwrite it.")
}
for (path in c(input_path, phenotype_path, eicu_path)) {
  if (!file.exists(path)) stop(sprintf("Required input is missing: %s", path))
}

imputation_object <- readRDS(input_path)
if (!inherits(imputation_object, "mids")) stop("Stage 5 object is not a mids object.")
if (imputation_object$m != 10L) stop("Expected ten imputations in the Stage 5 object.")

original_data <- imputation_object$data
if (nrow(original_data) != 5904L) stop("Expected 5,904 rows in the Stage 5 object.")
if (any(grepl("subject_id|stay_id|mortality|outcome", names(original_data), ignore.case = TRUE))) {
  stop("Stage 5 object unexpectedly contains an identifier or outcome field.")
}

variable_map <- c(
  neutrophil = "neutrophil_abs",
  lymphocyte = "lymphocyte_abs",
  platelet = "platelet_count",
  lactate = "lactate",
  glucose = "glucose",
  creatinine = "creatinine",
  bicarbonate = "bicarbonate"
)
required_variables <- c("cancer_phenotype", unname(variable_map))
if (!all(required_variables %in% names(original_data))) {
  stop("A required phenotype or primary variable is absent from the Stage 5 object.")
}

phenotype_aggregate <- utils::read.csv(phenotype_path, check.names = FALSE)
expected_phenotypes <- c(
  "hematologic",
  "metastatic solid",
  "nonmetastatic solid",
  "mixed"
)
if (!setequal(phenotype_aggregate$phenotype, expected_phenotypes)) {
  stop("Unexpected phenotype categories in the Stage 2 aggregate.")
}

phenotype_values <- as.character(original_data$cancer_phenotype)
if (anyNA(phenotype_values) || !all(phenotype_values %in% expected_phenotypes)) {
  stop("Unexpected or missing cancer phenotype in the Stage 5 object.")
}
phenotype_counts <- table(factor(phenotype_values, levels = expected_phenotypes))
aggregate_counts <- setNames(phenotype_aggregate$n, phenotype_aggregate$phenotype)
if (!identical(as.integer(phenotype_counts), as.integer(aggregate_counts[expected_phenotypes]))) {
  stop("Stage 5 phenotype counts do not reconcile to the Stage 2 aggregate.")
}

solid_levels <- c("metastatic solid", "nonmetastatic solid")
solid_index <- phenotype_values %in% solid_levels
solid_n <- sum(solid_index)
if (solid_n != sum(aggregate_counts[solid_levels])) {
  stop("Solid-tumor count failed aggregate reconciliation.")
}

blank_row <- function(
  data_source,
  section,
  metric,
  value_n = NA_integer_,
  denominator_n = NA_integer_,
  percent = NA_real_,
  available_n = NA_integer_,
  missing_n = NA_integer_,
  missing_percent = NA_real_,
  assessment,
  notes
) {
  data.frame(
    data_source = data_source,
    section = section,
    metric = metric,
    value_n = value_n,
    denominator_n = denominator_n,
    percent = percent,
    available_n = available_n,
    missing_n = missing_n,
    missing_percent = missing_percent,
    assessment = assessment,
    notes = notes,
    stringsAsFactors = FALSE
  )
}

rows <- list()
add_row <- function(row) rows[[length(rows) + 1L]] <<- row

total_n <- nrow(original_data)
add_row(blank_row(
  "MIMIC-IV",
  "cohort_size",
  "All cancer + sepsis with ICU LOS >=24h",
  value_n = total_n,
  denominator_n = total_n,
  percent = 100,
  assessment = "source_cohort",
  notes = "Fixed Stage 4A/Stage 5 cohort; no new database query."
))
add_row(blank_row(
  "MIMIC-IV",
  "cohort_size",
  "Included solid tumor cohort",
  value_n = solid_n,
  denominator_n = total_n,
  percent = 100 * solid_n / total_n,
  assessment = "included",
  notes = "Metastatic solid plus nonmetastatic solid; sepsis and ICU LOS >=24h already fixed."
))
for (phenotype in c("hematologic", "mixed")) {
  phenotype_n <- unname(phenotype_counts[phenotype])
  add_row(blank_row(
    "MIMIC-IV",
    "cohort_size",
    paste("Excluded", phenotype),
    value_n = phenotype_n,
    denominator_n = total_n,
    percent = 100 * phenotype_n / total_n,
    assessment = "excluded",
    notes = "Excluded to reduce hematologic marrow-related confounding of immune-cell measurements."
  ))
}

for (phenotype in solid_levels) {
  phenotype_n <- unname(phenotype_counts[phenotype])
  add_row(blank_row(
    "MIMIC-IV",
    "cancer_phenotype_distribution",
    phenotype,
    value_n = phenotype_n,
    denominator_n = solid_n,
    percent = 100 * phenotype_n / solid_n,
    assessment = "included_subtype",
    notes = "Distribution within the solid-tumor-restricted cohort."
  ))
}

solid_data <- original_data[solid_index, unname(variable_map), drop = FALSE]
for (display_name in names(variable_map)) {
  source_name <- variable_map[[display_name]]
  missing_n <- sum(is.na(solid_data[[source_name]]))
  available_n <- solid_n - missing_n
  add_row(blank_row(
    "MIMIC-IV",
    "primary_variable_missingness",
    display_name,
    denominator_n = solid_n,
    available_n = available_n,
    missing_n = missing_n,
    missing_percent = 100 * missing_n / solid_n,
    assessment = "original_pre_imputation_coverage",
    notes = sprintf("Original Stage 5 mids$data field: %s; imputed values not counted as available.", source_name)
  ))
}

solid_missing_percent <- vapply(unname(variable_map), function(variable) {
  100 * sum(is.na(solid_data[[variable]])) / solid_n
}, numeric(1))
highest_missing_variable <- names(which.max(setNames(solid_missing_percent, names(variable_map))))
highest_missing_percent <- max(solid_missing_percent)
add_row(blank_row(
  "MIMIC-IV",
  "feasibility_conclusion",
  "Solid-tumor-restricted phenotype discovery",
  value_n = solid_n,
  denominator_n = total_n,
  percent = 100 * solid_n / total_n,
  assessment = "potentially_feasible_requires_prespecified_model_phase",
  notes = sprintf(
    "Substantial cohort size, but original missingness remains material (highest: %s %.2f%%). This audit does not authorize a model.",
    highest_missing_variable,
    highest_missing_percent
  )
))

eicu_aggregate <- utils::read.csv(eicu_path, check.names = FALSE)
eicu_target <- eicu_aggregate[
  eicu_aggregate$row_type == "count" &
    eicu_aggregate$metric == "Cancer + sepsis with ICU LOS >=24h",
]
if (nrow(eicu_target) != 1L) stop("Expected one eICU all-cancer target count.")
add_row(blank_row(
  "eICU-CRD",
  "all_cancer_context",
  "Cancer + sepsis with ICU LOS >=24h",
  value_n = as.integer(eicu_target$value),
  assessment = "context_only_not_solid_specific",
  notes = "Previously exported Stage 2 all-cancer aggregate; cannot be treated as a solid-tumor count."
))

eicu_coverage <- eicu_aggregate[eicu_aggregate$row_type == "variable_coverage", ]
for (display_name in names(variable_map)) {
  source_row <- eicu_coverage[tolower(eicu_coverage$metric) == display_name, ]
  if (nrow(source_row) != 1L) {
    stop(sprintf("Missing eICU context coverage row for %s.", display_name))
  }
  add_row(blank_row(
    "eICU-CRD",
    "all_cancer_variable_context",
    display_name,
    denominator_n = as.integer(source_row$denominator_n),
    available_n = as.integer(source_row$value),
    missing_n = as.integer(source_row$missing_n),
    missing_percent = as.numeric(source_row$missing_percent),
    assessment = "context_only_not_solid_specific",
    notes = "Previously exported Stage 2 all-cancer coverage; no phenotype-stratified denominator is available."
  ))
}

add_row(blank_row(
  "eICU-CRD",
  "feasibility_conclusion",
  "Solid-tumor external-validation feasibility",
  assessment = "not_assessable_from_existing_aggregate_evidence",
  notes = paste(
    "The Stage 2 eICU output combines hematologic and solid malignancies.",
    "A separately reviewed aggregate-only phenotype-stratified query would be required; no database connection was made in Stage 6C."
  )
))

output <- do.call(rbind, rows)
numeric_percent_columns <- c("percent", "missing_percent")
output[numeric_percent_columns] <- lapply(
  output[numeric_percent_columns],
  function(x) ifelse(is.na(x), NA_real_, round(x, 2))
)

if (nrow(output) != 23L) stop("Unexpected Stage 6C output row count.")
if (sum(output$value_n[output$section == "cancer_phenotype_distribution"], na.rm = TRUE) != solid_n) {
  stop("Solid-tumor phenotype rows do not sum to the restricted cohort.")
}
if (any(
  output$available_n[output$section == "primary_variable_missingness"] +
    output$missing_n[output$section == "primary_variable_missingness"] != solid_n
)) {
  stop("A MIMIC variable-coverage row does not reconcile to the solid-tumor denominator.")
}

output_tmp <- tempfile("stage6c_solid_", tmpdir = dirname(output_path), fileext = ".csv")
utils::write.csv(output, output_tmp, row.names = FALSE, na = "")
verification <- utils::read.csv(output_tmp, check.names = FALSE)
if (nrow(verification) != nrow(output) || ncol(verification) != ncol(output)) {
  stop("Temporary Stage 6C CSV failed structural verification.")
}
if (!file.rename(output_tmp, output_path)) {
  stop("Could not finalize the Stage 6C aggregate table.")
}

cat("STAGE6C_SOLID_TUMOR_FEASIBILITY_COMPLETED\n")
cat(sprintf(
  "solid_n=%d of %d; highest_missing=%s %.2f%%; eicu_solid_status=not_assessable\n",
  solid_n,
  total_n,
  highest_missing_variable,
  highest_missing_percent
))
