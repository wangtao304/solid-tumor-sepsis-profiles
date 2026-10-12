#!/usr/bin/env Rscript

# Stage 9B: outcome association analysis for the frozen 2-profile solution.
# Identifiers are used in R memory only for deterministic linkage validation.
# No patient-level profile assignment, outcome record, or prediction is exported.

options(stringsAsFactors = FALSE, warn = 1)

if (!requireNamespace("mice", quietly = TRUE)) {
  stop("Required package 'mice' is not installed; no installation is authorized.")
}
if (!requireNamespace("survival", quietly = TRUE)) {
  stop("Required package 'survival' is not installed; no installation is authorized.")
}
if (!requireNamespace("rstpm2", quietly = TRUE)) {
  stop("Required package 'rstpm2' is not installed; no installation is authorized.")
}
suppressPackageStartupMessages(library(rstpm2))

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
sql_path <- file.path("sql", "30_profile_survival_input.sql")
results_path <- file.path("outputs", "tables", "30_profile_survival_results.csv")
rmst_path <- file.path("outputs", "tables", "31_profile_rmst.csv")
complete_case_path <- file.path(
  "outputs",
  "tables",
  "32_profile_complete_case_sensitivity.csv"
)
km_path <- file.path("outputs", "figures", "figure_profile_km.pdf")

required_paths <- c(mice_path, lpa_path, sql_path)
if (!all(file.exists(required_paths))) {
  stop(
    "Missing required Stage 9B input(s): ",
    paste(required_paths[!file.exists(required_paths)], collapse = ", ")
  )
}
dir.create(dirname(results_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(km_path), recursive = TRUE, showWarnings = FALSE)

expected_mice_variables <- c(
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
primary_profile_variables <- expected_mice_variables[1:6]
forbidden_adjustment_variables <- primary_profile_variables

mice_object <- readRDS(mice_path)
lpa_object <- readRDS(lpa_path)
if (
  !inherits(mice_object, "mids") ||
  mice_object$m != 10L ||
  nrow(mice_object$data) != 4380L ||
  !identical(names(mice_object$data), expected_mice_variables)
) {
  stop("Approved Stage 7B MICE object does not match the frozen Stage 9B schema.")
}
if (
  !identical(lpa_object$stage, "Stage 8A-1") ||
  !identical(lpa_object$configuration$variables, primary_profile_variables) ||
  !identical(lpa_object$outcome_variables_used, FALSE) ||
  !identical(lpa_object$survival_information_used, FALSE) ||
  !identical(lpa_object$patient_level_profile_labels_exported, FALSE)
) {
  stop("Stage 8A-1 LPA object failed Stage 9B safeguard checks.")
}

# Reconstruct the frozen Stage 9A consensus assignment without refitting LPA.
profile_assignments <- matrix(
  NA_integer_,
  nrow = nrow(mice_object$data),
  ncol = 10L
)
mean_aligned_posterior <- matrix(
  0,
  nrow = nrow(mice_object$data),
  ncol = 2L
)
for (imputation in 1:10) {
  key <- sprintf("imp%02d_k2", imputation)
  fit_record <- lpa_object$fit_records[[key]]
  alignment <- lpa_object$alignments[[key]]
  if (
    is.null(fit_record) ||
    is.null(alignment) ||
    length(alignment) != 2L ||
    length(fit_record$model$classification) != nrow(mice_object$data)
  ) {
    stop("Missing or malformed frozen 2-profile fit for ", key, ".")
  }
  profile_assignments[, imputation] <- match(
    fit_record$model$classification,
    alignment
  )
  mean_aligned_posterior <- mean_aligned_posterior +
    fit_record$model$posterior[, alignment, drop = FALSE] / 10
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
if (
  anyNA(consensus_profile) ||
  !identical(as.integer(tabulate(consensus_profile, nbins = 2L)), c(3197L, 1183L)) ||
  sum(profile_1_votes == profile_2_votes) != 26L
) {
  stop("Frozen Stage 9A consensus assignment could not be reproduced.")
}

# Stream the outcome/linkage input directly from PostgreSQL into R memory.
psql_bin <- Sys.getenv(
  "MIMIC_PSQL_BIN",
  unset = ""
)
db_host <- Sys.getenv("MIMIC_DB_HOST", unset = "")
db_port <- Sys.getenv("MIMIC_DB_PORT", unset = "")
db_name <- Sys.getenv("MIMIC_DB_NAME", unset = "")
db_user <- Sys.getenv("MIMIC_DB_USER", unset = "")
pgpass_file <- Sys.getenv(
  "MIMIC_PGPASS_FILE",
  unset = ""
)
if (!file.exists(psql_bin) || !file.exists(pgpass_file) ||
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
  "-f", sql_path
)
csv_lines <- system2(
  psql_bin,
  args = psql_args,
  stdout = TRUE,
  stderr = ""
)
psql_status <- attr(csv_lines, "status")
if (!is.null(psql_status) && psql_status != 0L) {
  stop("Read-only Stage 9B psql query failed.")
}
if (length(csv_lines) < 2L) stop("Stage 9B query returned no data rows.")

outcome_data <- utils::read.csv(
  text = paste(csv_lines, collapse = "\n"),
  header = TRUE,
  na.strings = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
rm(csv_lines)

expected_outcome_columns <- c(
  "stay_id",
  "age",
  "sex",
  "cancer_phenotype",
  "charlson_comorbidity_index",
  "sofa",
  "sapsii",
  "vasopressor_0_24",
  "mechanical_ventilation_0_24",
  "event_28d",
  "time_28d_days",
  "event_icu",
  "time_icu_days",
  "event_hospital",
  "time_hospital_days",
  "complete_case_primary"
)
if (
  !identical(names(outcome_data), expected_outcome_columns) ||
  nrow(outcome_data) != nrow(mice_object$data) ||
  anyNA(outcome_data) ||
  anyDuplicated(outcome_data$stay_id) ||
  !all(diff(outcome_data$stay_id) > 0)
) {
  stop("Stage 9B database input failed structural or ordering checks.")
}

normalize_boolean <- function(values) {
  normalized <- tolower(as.character(values))
  if (!all(normalized %in% c("t", "f", "true", "false", "1", "0", "yes", "no"))) {
    stop("Unexpected Boolean coding in Stage 9B database input.")
  }
  normalized %in% c("t", "true", "1", "yes")
}
numeric_equal <- function(database_values, local_values, tolerance = 1e-8) {
  database_values <- as.numeric(database_values)
  local_values <- as.numeric(local_values)
  identical(is.na(database_values), is.na(local_values)) &&
    all(abs(database_values[!is.na(database_values)] -
      local_values[!is.na(local_values)]) <= tolerance)
}

local_data <- mice_object$data
link_checks <- c(
  age = numeric_equal(outcome_data$age, local_data$age),
  sex = identical(as.character(outcome_data$sex), as.character(local_data$sex)),
  cancer_phenotype = identical(
    as.character(outcome_data$cancer_phenotype),
    as.character(local_data$cancer_phenotype)
  ),
  charlson = numeric_equal(
    outcome_data$charlson_comorbidity_index,
    local_data$charlson_comorbidity_index
  ),
  sofa = numeric_equal(outcome_data$sofa, local_data$sofa),
  sapsii = numeric_equal(outcome_data$sapsii, local_data$sapsii),
  vasopressor = identical(
    normalize_boolean(outcome_data$vasopressor_0_24),
    as.character(local_data$vasopressor_0_24) == "Yes"
  ),
  mechanical_ventilation = identical(
    normalize_boolean(outcome_data$mechanical_ventilation_0_24),
    as.character(local_data$mechanical_ventilation_0_24) == "Yes"
  ),
  complete_case = identical(
    normalize_boolean(outcome_data$complete_case_primary),
    stats::complete.cases(local_data[primary_profile_variables])
  )
)
if (!all(link_checks)) {
  stop(
    "Stage 9B row-order linkage validation failed for: ",
    paste(names(link_checks)[!link_checks], collapse = ", ")
  )
}

event_columns <- c("event_28d", "event_icu", "event_hospital")
time_columns <- c("time_28d_days", "time_icu_days", "time_hospital_days")
outcome_data[event_columns] <- lapply(
  outcome_data[event_columns],
  as.integer
)
outcome_data[time_columns] <- lapply(
  outcome_data[time_columns],
  as.numeric
)
if (
  !all(unlist(outcome_data[event_columns], use.names = FALSE) %in% 0:1) ||
  any(!is.finite(unlist(outcome_data[time_columns], use.names = FALSE))) ||
  any(unlist(outcome_data[time_columns], use.names = FALSE) <= 0) ||
  any(outcome_data$time_28d_days > 28 + 1e-8) ||
  any(outcome_data$time_28d_days[outcome_data$event_28d == 0L] != 28) ||
  !identical(as.integer(colSums(outcome_data[event_columns])), c(1537L, 711L, 1120L))
) {
  stop("Stage 9B event or follow-up time validation failed.")
}

analysis_data <- data.frame(
  profile = factor(
    consensus_profile,
    levels = 1:2,
    labels = c("Profile 1", "Profile 2")
  ),
  age = as.numeric(local_data$age),
  sex = factor(local_data$sex, levels = c("F", "M")),
  charlson = as.numeric(local_data$charlson_comorbidity_index),
  cancer_phenotype = factor(
    local_data$cancer_phenotype,
    levels = c("metastatic solid", "nonmetastatic solid")
  ),
  sofa = as.numeric(local_data$sofa),
  event_28d = outcome_data$event_28d,
  time_28d_days = outcome_data$time_28d_days,
  event_icu = outcome_data$event_icu,
  time_icu_days = outcome_data$time_icu_days,
  event_hospital = outcome_data$event_hospital,
  time_hospital_days = outcome_data$time_hospital_days,
  complete_case_primary = normalize_boolean(outcome_data$complete_case_primary)
)
rm(outcome_data)
if (
  anyNA(analysis_data) ||
  sum(analysis_data$complete_case_primary) != 2649L
) {
  stop("Stage 9B analysis frame failed completeness checks.")
}
analysis_data$age_centered_per_10_years <-
  (analysis_data$age - mean(analysis_data$age)) / 10
analysis_data$charlson_centered <-
  analysis_data$charlson - mean(analysis_data$charlson)
analysis_data$sofa_centered <-
  analysis_data$sofa - mean(analysis_data$sofa)
analysis_data$sapsii_centered <-
  analysis_data$sapsii - mean(analysis_data$sapsii)

outcome_specifications <- list(
  list(
    outcome = "28-day mortality",
    event = "event_28d",
    time = "time_28d_days",
    time_definition = "ICU admission to death; administratively censored at 28 days"
  ),
  list(
    outcome = "ICU mortality",
    event = "event_icu",
    time = "time_icu_days",
    time_definition = "ICU admission to in-ICU death; censored at ICU discharge"
  ),
  list(
    outcome = "Hospital mortality",
    event = "event_hospital",
    time = "time_hospital_days",
    time_definition = "ICU admission to in-hospital death; censored at hospital discharge"
  )
)
model_specifications <- data.frame(
  model = c("Model 1", "Model 2", "Model 3", "Model 4", "Model 4 sensitivity (SAPS-II)"),
  rhs = c(
    "profile",
    "profile + age_centered_per_10_years + sex",
    paste(
      "profile + age_centered_per_10_years + sex +",
      "charlson_centered + cancer_phenotype"
    ),
    paste(
      "profile + age_centered_per_10_years + sex +",
      "charlson_centered + cancer_phenotype + sofa_centered"
    ),
    paste(
      "profile + age_centered_per_10_years + sex +",
      "charlson_centered + cancer_phenotype + sapsii_centered"
    )
  ),
  stringsAsFactors = FALSE
)
if (any(vapply(forbidden_adjustment_variables, function(variable) {
  any(grepl(paste0("(^|[ +])", variable, "($|[ +])"), model_specifications$rhs))
}, logical(1L)))) {
  stop("A forbidden phenotype indicator entered a Stage 9B adjustment model.")
}

fit_flexible_parametric_models <- function(data, population_label) {
  rows <- list()
  for (outcome_specification in outcome_specifications) {
    for (model_index in seq_len(nrow(model_specifications))) {
      model_name <- model_specifications$model[model_index]
      rhs <- model_specifications$rhs[model_index]
      formula_text <- sprintf(
        "survival::Surv(%s, %s) ~ %s",
        outcome_specification$time,
        outcome_specification$event,
        rhs
      )
      model_formula <- stats::as.formula(formula_text)
      fit <- stpm2(
        model_formula,
        data = data,
        df = 4L,
        scale = "hazard",
        control = list(
          maxit = 10000L,
          use.gr = TRUE,
          reltol.search = 1e-10,
          reltol.final = 1e-12,
          reltol.outer = 1e-12
        )
      )
      fit_details <- fit@details
      coefficient_names <- names(stats::coef(fit))
      profile_term <- grep("^profileProfile 2$", coefficient_names, value = TRUE)
      if (
        length(profile_term) != 1L ||
        !identical(as.integer(fit_details$convergence), 0L) ||
        !is.finite(fit_details$maxgrad) ||
        fit_details$maxgrad > 1e-2
      ) {
        stop(
          sprintf(
            paste0(
              "Flexible parametric model failed convergence/extraction checks: ",
              "%s, %s, %s; profile_terms=%d; convergence=%s; maxgrad=%.8g; ",
              "coefficient_names=%s"
            ),
            population_label,
            outcome_specification$outcome,
            model_name,
            length(profile_term),
            paste(fit_details$convergence, collapse = "|"),
            fit_details$maxgrad,
            paste(coefficient_names, collapse = "|")
          )
        )
      }
      beta <- stats::coef(fit)[profile_term]
      variance <- stats::vcov(fit)[profile_term, profile_term]
      standard_error <- sqrt(variance)
      z_value <- beta / standard_error
      rows[[length(rows) + 1L]] <- data.frame(
        analysis_population = population_label,
        outcome = outcome_specification$outcome,
        model = model_name,
        n = nrow(data),
        events = sum(data[[outcome_specification$event]]),
        reference_profile = "Profile 1: Lower inflammatory-metabolic stress",
        comparison_profile = "Profile 2: Higher inflammatory-metabolic stress",
        effect_measure = "Hazard ratio",
        estimate = exp(beta),
        confidence_interval_lower = exp(beta - stats::qnorm(0.975) * standard_error),
        confidence_interval_upper = exp(beta + stats::qnorm(0.975) * standard_error),
        p_value = 2 * stats::pnorm(-abs(z_value)),
        baseline_spline_df = 4L,
        model_formula = formula_text,
        covariate_parameterization = paste(
          "Age centered and scaled per 10 years;",
          "Charlson and SOFA mean-centered; categorical variables unchanged"
        ),
        time_definition = outcome_specification$time_definition,
        convergence_code = as.integer(fit_details$convergence),
        maximum_absolute_gradient = fit_details$maxgrad,
        aic = stats::AIC(fit),
        stringsAsFactors = FALSE
      )
    }
  }
  output <- do.call(rbind, rows)
  row.names(output) <- NULL
  output
}

main_results <- fit_flexible_parametric_models(
  analysis_data,
  "Full cohort with frozen consensus assignment"
)
complete_case_data <- analysis_data[analysis_data$complete_case_primary, ]
complete_case_results <- fit_flexible_parametric_models(
  complete_case_data,
  "Complete-case phenotype-variable sensitivity cohort"
)

# Unadjusted 28-day Kaplan-Meier curves and RMST at tau=28 days.
km_fit <- survival::survfit(
  survival::Surv(time_28d_days, event_28d) ~ profile,
  data = analysis_data,
  conf.type = "log-log"
)
rmst_summary <- summary(km_fit, rmean = 28)$table
if (
  !all(c("rmean", "se(rmean)") %in% colnames(rmst_summary)) ||
  nrow(rmst_summary) != 2L
) {
  stop("Could not extract 28-day RMST from the Kaplan-Meier fit.")
}
rmst_estimates <- as.numeric(rmst_summary[, "rmean"])
rmst_standard_errors <- as.numeric(rmst_summary[, "se(rmean)"])
rmst_labels <- c(
  "Profile 1: Lower inflammatory-metabolic stress",
  "Profile 2: Higher inflammatory-metabolic stress"
)
rmst_rows <- lapply(1:2, function(profile) {
  data.frame(
    analysis_population = "Full cohort with frozen consensus assignment",
    outcome = "28-day mortality",
    horizon_days = 28,
    estimand = paste0("RMST: ", rmst_labels[profile]),
    estimate_days = rmst_estimates[profile],
    standard_error_days = rmst_standard_errors[profile],
    confidence_interval_lower_days = rmst_estimates[profile] -
      stats::qnorm(0.975) * rmst_standard_errors[profile],
    confidence_interval_upper_days = rmst_estimates[profile] +
      stats::qnorm(0.975) * rmst_standard_errors[profile],
    p_value = NA_real_,
    n = sum(analysis_data$profile == levels(analysis_data$profile)[profile]),
    events = sum(
      analysis_data$event_28d[
        analysis_data$profile == levels(analysis_data$profile)[profile]
      ]
    ),
    method = "Unadjusted Kaplan-Meier area under survival curve",
    stringsAsFactors = FALSE
  )
})
rmst_difference <- rmst_estimates[2L] - rmst_estimates[1L]
rmst_difference_se <- sqrt(sum(rmst_standard_errors^2))
rmst_rows[[3L]] <- data.frame(
  analysis_population = "Full cohort with frozen consensus assignment",
  outcome = "28-day mortality",
  horizon_days = 28,
  estimand = "RMST difference: Profile 2 minus Profile 1",
  estimate_days = rmst_difference,
  standard_error_days = rmst_difference_se,
  confidence_interval_lower_days = rmst_difference -
    stats::qnorm(0.975) * rmst_difference_se,
  confidence_interval_upper_days = rmst_difference +
    stats::qnorm(0.975) * rmst_difference_se,
  p_value = 2 * stats::pnorm(-abs(rmst_difference / rmst_difference_se)),
  n = nrow(analysis_data),
  events = sum(analysis_data$event_28d),
  method = "Unadjusted difference; independent-group Greenwood standard errors",
  stringsAsFactors = FALSE
)
rmst_results <- do.call(rbind, rmst_rows)
row.names(rmst_results) <- NULL

# Publication-ready Kaplan-Meier figure with aggregate risk table.
km_tmp <- tempfile(
  "stage9b_km_",
  tmpdir = dirname(km_path),
  fileext = ".pdf"
)
profile_colors <- c("#0072B2", "#D55E00")
risk_times <- c(0, 7, 14, 21, 28)
risk_summary <- summary(km_fit, times = risk_times, extend = TRUE)
strata_names <- names(km_fit$strata)
risk_matrix <- matrix(
  NA_integer_,
  nrow = 2L,
  ncol = length(risk_times),
  dimnames = list(c("Profile 1", "Profile 2"), as.character(risk_times))
)
for (profile_index in 1:2) {
  rows <- as.character(risk_summary$strata) == strata_names[profile_index]
  risk_matrix[profile_index, ] <- risk_summary$n.risk[rows]
}
if (anyNA(risk_matrix)) stop("Kaplan-Meier risk table extraction failed.")

grDevices::pdf(
  km_tmp,
  width = 8.6,
  height = 6.8,
  family = "sans",
  title = "Kaplan-Meier survival by inflammatory-metabolic profile",
  useDingbats = FALSE
)
graphics::layout(matrix(c(1, 2), nrow = 2L), heights = c(4.5, 1.25))
graphics::par(
  mar = c(1.6, 4.6, 4.5, 1.0),
  mgp = c(2.7, 0.8, 0),
  las = 1,
  xaxs = "i",
  yaxs = "i"
)
graphics::plot(
  km_fit,
  col = profile_colors,
  lwd = 2.4,
  lty = 1,
  conf.int = FALSE,
  mark.time = TRUE,
  xlim = c(0, 28),
  ylim = c(0, 1),
  xaxt = "n",
  xlab = "",
  ylab = "Survival probability",
  main = "28-day survival by inflammatory-metabolic profile",
  cex.main = 1.2,
  font.main = 2,
  col.axis = "#333333",
  col.lab = "#222222"
)
graphics::axis(1, at = risk_times, labels = FALSE)
graphics::abline(h = seq(0, 1, by = 0.2), col = "#E6E6E6", lwd = 0.7)
graphics::mtext(
  "Profile assignment was frozen before outcome access",
  side = 3,
  line = 1.0,
  cex = 0.84,
  col = "#555555"
)
profile_events <- vapply(1:2, function(profile) {
  sum(analysis_data$event_28d[analysis_data$profile == paste("Profile", profile)])
}, integer(1L))
profile_sizes <- as.integer(table(analysis_data$profile))
legend_labels <- sprintf(
  "Profile %d: %s (n=%d; deaths=%d)",
  1:2,
  c("Lower stress", "Higher stress"),
  profile_sizes,
  profile_events
)
graphics::legend(
  "bottomleft",
  legend = legend_labels,
  col = profile_colors,
  lwd = 2.4,
  bty = "n",
  cex = 0.86
)

graphics::par(mar = c(3.5, 4.6, 0.3, 1.0), xaxs = "i", yaxs = "i")
graphics::plot(
  NA,
  xlim = c(0, 28),
  ylim = c(0.3, 2.8),
  axes = FALSE,
  xlab = "Days since ICU admission",
  ylab = ""
)
graphics::text(-1.0, 2.45, "Number at risk", adj = 1, font = 2, xpd = NA, cex = 0.83)
for (profile_index in 1:2) {
  y <- 2.15 - (profile_index - 1L) * 0.65
  graphics::text(
    -1.0,
    y,
    paste0("Profile ", profile_index),
    adj = 1,
    xpd = NA,
    col = profile_colors[profile_index],
    font = 2,
    cex = 0.82
  )
  graphics::text(
    risk_times,
    y,
    format(risk_matrix[profile_index, ], big.mark = ","),
    cex = 0.82
  )
}
graphics::axis(1, at = risk_times, labels = risk_times, col.axis = "#333333")
graphics::mtext(
  "Unadjusted Kaplan-Meier estimates; administrative censoring at 28 days.",
  side = 1,
  line = 2.55,
  adj = 0,
  cex = 0.70,
  col = "#555555"
)
grDevices::dev.off()

numeric_result_columns <- c(
  "estimate",
  "confidence_interval_lower",
  "confidence_interval_upper",
  "p_value",
  "maximum_absolute_gradient",
  "aic"
)
main_results[numeric_result_columns] <- lapply(
  main_results[numeric_result_columns],
  function(values) signif(values, 8L)
)
complete_case_results[numeric_result_columns] <- lapply(
  complete_case_results[numeric_result_columns],
  function(values) signif(values, 8L)
)
numeric_rmst_columns <- c(
  "estimate_days",
  "standard_error_days",
  "confidence_interval_lower_days",
  "confidence_interval_upper_days",
  "p_value"
)
rmst_results[numeric_rmst_columns] <- lapply(
  rmst_results[numeric_rmst_columns],
  function(values) signif(values, 8L)
)

results_tmp <- tempfile(
  "stage9b_results_",
  tmpdir = dirname(results_path),
  fileext = ".csv"
)
rmst_tmp <- tempfile(
  "stage9b_rmst_",
  tmpdir = dirname(rmst_path),
  fileext = ".csv"
)
complete_case_tmp <- tempfile(
  "stage9b_cc_",
  tmpdir = dirname(complete_case_path),
  fileext = ".csv"
)
utils::write.csv(main_results, results_tmp, row.names = FALSE, na = "")
utils::write.csv(rmst_results, rmst_tmp, row.names = FALSE, na = "")
utils::write.csv(
  complete_case_results,
  complete_case_tmp,
  row.names = FALSE,
  na = ""
)

verification_results <- utils::read.csv(
  results_tmp,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
verification_rmst <- utils::read.csv(
  rmst_tmp,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
verification_complete_case <- utils::read.csv(
  complete_case_tmp,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
prohibited_output_fields <- c(
  "subject_id",
  "stay_id",
  "profile_assignment",
  "posterior_probability"
)
if (
  nrow(verification_results) != 12L ||
  nrow(verification_rmst) != 3L ||
  nrow(verification_complete_case) != 12L ||
  any(unlist(lapply(
    list(
      verification_results,
      verification_rmst,
      verification_complete_case
    ),
    function(data) any(names(data) %in% prohibited_output_fields)
  ))) ||
  any(verification_results$convergence_code != 0L) ||
  any(verification_complete_case$convergence_code != 0L) ||
  any(verification_results$n != 4380L) ||
  any(verification_complete_case$n != 2649L) ||
  file.info(km_tmp)$size < 1000
) {
  stop("Temporary Stage 9B aggregate outputs failed structural verification.")
}

if (!file.copy(results_tmp, results_path, overwrite = TRUE)) {
  stop("Could not finalize Stage 9B survival-result table.")
}
if (!file.copy(rmst_tmp, rmst_path, overwrite = TRUE)) {
  stop("Could not finalize Stage 9B RMST table.")
}
if (!file.copy(complete_case_tmp, complete_case_path, overwrite = TRUE)) {
  stop("Could not finalize Stage 9B complete-case sensitivity table.")
}
if (!file.copy(km_tmp, km_path, overwrite = TRUE)) {
  stop("Could not finalize Stage 9B Kaplan-Meier figure.")
}
unlink(c(results_tmp, rmst_tmp, complete_case_tmp, km_tmp))

cat("STAGE9B_SURVIVAL_ANALYSIS_COMPLETED\n")
cat(sprintf(
  paste0(
    "full_n=%d complete_case_n=%d events_28d=%d events_icu=%d ",
    "events_hospital=%d models=%d no_patient_level_export=TRUE\n"
  ),
  nrow(analysis_data),
  nrow(complete_case_data),
  sum(analysis_data$event_28d),
  sum(analysis_data$event_icu),
  sum(analysis_data$event_hospital),
  nrow(main_results) + nrow(complete_case_results)
))
