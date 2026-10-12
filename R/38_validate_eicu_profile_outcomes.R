#!/usr/bin/env Rscript

# Stage 10C-2: eICU hospital-mortality external validation.
# Uses frozen Stage 10C-1 profile assignments without classifier retraining,
# threshold changes, or LPA refitting. patientunitstayid is transient linkage
# state only and is never written to an output.

options(stringsAsFactors = FALSE, warn = 1)

expected_n <- 1407L
expected_complete_case_n <- 741L
primary_variables <- c(
  "wbc", "platelet_count", "lactate",
  "glucose", "creatinine", "bicarbonate"
)
auxiliary_variables <- c("age", "sex", "cancer_phenotype")
approved_variables <- c(primary_variables, auxiliary_variables)

stage10c1_path <- file.path(
  "data", "processed",
  "eicu_inflammatory_metabolic_profile_transfer_mi10_seed20260724.rds"
)
transfer_path <- file.path(
  "data", "processed",
  "inflammatory_metabolic_profile_transfer_logistic_mi10_seed20260723.rds"
)
sql_path <- file.path("sql", "38_eicu_profile_outcome_input.sql")
result_path <- file.path(
  "outputs", "tables", "38_eicu_profile_outcome_results.csv"
)
complete_case_path <- file.path(
  "outputs", "tables", "39_eicu_complete_case_sensitivity.csv"
)
figure_path <- file.path(
  "outputs", "figures", "figure_eicu_external_validation.pdf"
)

for (path in c(stage10c1_path, transfer_path, sql_path)) {
  if (!file.exists(path)) stop("Missing required Stage 10C-2 input: ", path)
}
for (path in c(result_path, complete_case_path, figure_path)) {
  if (file.exists(path)) stop("Stage 10C-2 output already exists; refusing to overwrite: ", path)
  if (!dir.exists(dirname(path))) stop("Controlled output directory is missing: ", dirname(path))
}

# Validate the frozen assignment and transfer contracts before outcome access.
stage10c1 <- readRDS(stage10c1_path)
transfer_model <- readRDS(transfer_path)
if (
  !identical(stage10c1$stage, "Stage 10C-1") ||
  !inherits(stage10c1$imputation, "mids") ||
  stage10c1$imputation$m != 10L ||
  stage10c1$imputation$iteration != 20L ||
  nrow(stage10c1$imputation$data) != expected_n ||
  !identical(names(stage10c1$imputation$data), approved_variables) ||
  !identical(dim(stage10c1$profile_assignment), c(expected_n, 10L)) ||
  anyNA(stage10c1$profile_assignment) ||
  !all(stage10c1$profile_assignment %in% 1:2) ||
  !identical(stage10c1$metadata$transfer_model_retrained, FALSE) ||
  !identical(stage10c1$metadata$eicu_recalibration_performed, FALSE) ||
  !identical(stage10c1$metadata$patient_identifiers_stored, FALSE) ||
  !identical(stage10c1$metadata$mortality_or_outcome_fields_read_or_stored, FALSE) ||
  !identical(transfer_model$variables, primary_variables) ||
  !identical(as.numeric(transfer_model$classification_threshold), 0.5) ||
  !identical(transfer_model$eicu_data_used_for_training, FALSE) ||
  !identical(transfer_model$lpa_refitted, FALSE) ||
  !identical(
    unname(tools::md5sum(transfer_path)),
    unname(stage10c1$metadata$transfer_model_md5)
  )
) {
  stop("Frozen Stage 10C-1 assignment or Stage 10B transfer contract failed validation.")
}

# Stdin mode keeps all patient-level linkage data in an in-memory pipe and lets
# the host shell supply its existing password-file context without copying it.
if (identical(Sys.getenv("STAGE10C2_INPUT_STDIN", unset = "0"), "1")) {
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
    stop("PostgreSQL failed while streaming the Stage 10C-2 linkage input.")
  }
}
if (length(csv_lines) < 2L) stop("The Stage 10C-2 linkage query returned no rows.")

expected_columns <- c(
  "patientunitstayid", approved_variables,
  "hospitaldischargestatus", "hospitaldischargeoffset"
)
database_input <- utils::read.csv(
  text = paste(csv_lines, collapse = "\n"),
  header = TRUE,
  na.strings = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
rm(csv_lines)
if (
  !identical(names(database_input), expected_columns) ||
  nrow(database_input) != expected_n ||
  anyNA(database_input$patientunitstayid) ||
  length(unique(database_input$patientunitstayid)) != expected_n ||
  anyDuplicated(database_input$patientunitstayid)
) {
  stop("Stage 10C-2 cohort uniqueness or input-column validation failed.")
}

database_input[c(primary_variables, "age", "hospitaldischargeoffset")] <- lapply(
  database_input[c(primary_variables, "age", "hospitaldischargeoffset")],
  as.numeric
)
database_input$sex <- factor(database_input$sex, levels = c("F", "M"))
database_input$cancer_phenotype <- factor(
  database_input$cancer_phenotype,
  levels = c("metastatic solid", "nonmetastatic solid")
)
if (anyNA(database_input[c("age", "sex", "cancer_phenotype")])) {
  stop("An approved adjustment variable is missing or invalid.")
}

numeric_equal <- function(database_values, local_values, tolerance = 1e-8) {
  database_values <- as.numeric(database_values)
  local_values <- as.numeric(local_values)
  identical(is.na(database_values), is.na(local_values)) &&
    all(abs(
      database_values[!is.na(database_values)] -
        local_values[!is.na(local_values)]
    ) <= tolerance)
}
local_data <- stage10c1$imputation$data
linkage_checks <- c(
  setNames(vapply(c(primary_variables, "age"), function(variable) {
    numeric_equal(database_input[[variable]], local_data[[variable]])
  }, logical(1L)), c(primary_variables, "age")),
  sex = identical(as.character(database_input$sex), as.character(local_data$sex)),
  cancer_phenotype = identical(
    as.character(database_input$cancer_phenotype),
    as.character(local_data$cancer_phenotype)
  ),
  complete_case = identical(
    complete.cases(database_input[primary_variables]),
    complete.cases(local_data[primary_variables])
  )
)
if (!all(linkage_checks)) {
  stop("Profile-outcome row-order linkage failed: ", paste(names(linkage_checks)[!linkage_checks], collapse = ", "))
}

status_value <- trimws(as.character(database_input$hospitaldischargestatus))
status_value[is.na(status_value) | !nzchar(status_value)] <- NA_character_
unexpected_status <- !is.na(status_value) & !status_value %in% c("Alive", "Expired")
if (any(unexpected_status)) {
  stop("Hospital discharge status contains an unexpected nonmissing category.")
}
outcome_available <- !is.na(status_value)
outcome_available_n <- sum(outcome_available)
outcome_missing_n <- sum(!outcome_available)
hospital_mortality <- ifelse(
  outcome_available,
  as.integer(status_value == "Expired"),
  NA_integer_
)
hospital_death_n <- sum(hospital_mortality, na.rm = TRUE)
if (hospital_death_n <= 0L || hospital_death_n >= outcome_available_n) {
  stop("Hospital mortality outcome has an invalid event count.")
}

# eICU has no standalone death timestamp. hospitaldischargeoffset is a discharge
# disposition time, so it is not treated as a reliable death time for survival
# modeling. The aggregate below documents the observable approximation only.
observed_expired_by_28d_n <- sum(
  hospital_mortality == 1L &
    !is.na(database_input$hospitaldischargeoffset) &
    database_input$hospitaldischargeoffset > 0 &
    database_input$hospitaldischargeoffset <= 28 * 24 * 60,
  na.rm = TRUE
)
reliable_death_time_available <- FALSE

# Delete identifiers and disposition strings immediately after linkage/outcome
# construction. They are never incorporated into a model frame or output.
database_input$patientunitstayid <- NULL
database_input$hospitaldischargestatus <- NULL

model_specification <- data.frame(
  model = c("Model 1", "Model 2", "Model 3"),
  rhs = c(
    "profile",
    "profile + age10 + sex",
    "profile + age10 + sex + cancer_phenotype"
  ),
  adjustment = c(
    "Profile only",
    "+ age + sex",
    "+ cancer phenotype"
  ),
  stringsAsFactors = FALSE
)
if (any(vapply(primary_variables, function(variable) {
  any(grepl(variable, model_specification$rhs, fixed = TRUE))
}, logical(1L)))) {
  stop("A prohibited profile-defining laboratory variable entered an outcome model.")
}

fit_logistic <- function(profile, row_index = which(outcome_available)) {
  if (length(profile) == expected_n) {
    profile_value <- profile[row_index]
  } else if (length(profile) == length(row_index)) {
    profile_value <- profile
  } else {
    stop("Profile vector length does not match the logistic model index.")
  }
  data <- data.frame(
    event = hospital_mortality[row_index],
    profile = factor(
      profile_value,
      levels = 1:2,
      labels = c("Profile 1", "Profile 2")
    ),
    age10 = (as.numeric(local_data$age[row_index]) - mean(as.numeric(local_data$age[row_index]))) / 10,
    sex = factor(local_data$sex[row_index], levels = c("F", "M")),
    cancer_phenotype = factor(
      local_data$cancer_phenotype[row_index],
      levels = c("metastatic solid", "nonmetastatic solid")
    )
  )
  if (anyNA(data)) stop("A logistic model frame contains missing values.")
  fits <- vector("list", nrow(model_specification))
  for (model_index in seq_len(nrow(model_specification))) {
    formula <- stats::as.formula(paste("event ~", model_specification$rhs[model_index]))
    fit <- stats::glm(formula, data = data, family = stats::binomial())
    term <- grep("^profileProfile 2$", names(stats::coef(fit)), value = TRUE)
    if (
      !isTRUE(fit$converged) || length(term) != 1L ||
      any(!is.finite(stats::coef(fit))) ||
      !is.finite(stats::vcov(fit)[term, term]) ||
      stats::vcov(fit)[term, term] <= 0
    ) {
      stop("A prespecified logistic model failed convergence or coefficient validation.")
    }
    fits[[model_index]] <- list(
      estimate = unname(stats::coef(fit)[term]),
      variance = unname(stats::vcov(fit)[term, term]),
      aic = stats::AIC(fit)
    )
  }
  fits
}

pool_scalar <- function(estimates, variances) {
  m <- length(estimates)
  qbar <- mean(estimates)
  ubar <- mean(variances)
  between <- stats::var(estimates)
  total <- ubar + (1 + 1 / m) * between
  if (!is.finite(total) || total <= 0) stop("Invalid Rubin pooled variance.")
  if (!is.finite(between) || between < 1e-15) {
    degrees_freedom <- Inf
    critical <- stats::qnorm(0.975)
    p_value <- 2 * stats::pnorm(-abs(qbar / sqrt(total)))
  } else {
    degrees_freedom <- (m - 1) * (1 + ubar / ((1 + 1 / m) * between))^2
    critical <- stats::qt(0.975, df = degrees_freedom)
    p_value <- 2 * stats::pt(-abs(qbar / sqrt(total)), df = degrees_freedom)
  }
  list(
    estimate = qbar,
    total_variance = total,
    between_variance = between,
    degrees_freedom = degrees_freedom,
    fraction_missing_information = ((1 + 1 / m) * between) / total,
    confidence_lower = qbar - critical * sqrt(total),
    confidence_upper = qbar + critical * sqrt(total),
    p_value = p_value
  )
}

mi_fits <- lapply(seq_len(10L), function(imputation_index) {
  fit_logistic(stage10c1$profile_assignment[, imputation_index])
})

result_rows <- vector("list", nrow(model_specification))
for (model_index in seq_len(nrow(model_specification))) {
  estimates <- vapply(mi_fits, function(x) x[[model_index]]$estimate, numeric(1L))
  variances <- vapply(mi_fits, function(x) x[[model_index]]$variance, numeric(1L))
  aic_values <- vapply(mi_fits, function(x) x[[model_index]]$aic, numeric(1L))
  pooled <- pool_scalar(estimates, variances)
  result_rows[[model_index]] <- data.frame(
    section = "association",
    outcome = "Hospital mortality",
    analysis_set = "MICE primary",
    model = model_specification$model[model_index],
    adjustment = model_specification$adjustment[model_index],
    model_family = "Logistic regression; Rubin pooled across 10 imputations",
    n = outcome_available_n,
    events = hospital_death_n,
    reference = "Profile 1",
    exposure = "Profile 2",
    effect_measure = "Odds ratio",
    estimate = exp(pooled$estimate),
    confidence_interval_lower = exp(pooled$confidence_lower),
    confidence_interval_upper = exp(pooled$confidence_upper),
    p_value = pooled$p_value,
    between_imputation_variance = pooled$between_variance,
    fraction_missing_information = pooled$fraction_missing_information,
    degrees_freedom = pooled$degrees_freedom,
    mean_aic = mean(aic_values),
    status = "ESTIMATED",
    notes = sprintf(
      "Frozen imputation-specific profile assignment; %d record(s) with unknown hospital status excluded; no classifier refit, recalibration, threshold change, or laboratory adjustment.",
      outcome_missing_n
    ),
    stringsAsFactors = FALSE
  )
}
primary_results <- do.call(rbind, result_rows)
not_estimable_row <- data.frame(
  section = "outcome_time_feasibility",
  outcome = "28-day in-hospital mortality with live-discharge censoring",
  analysis_set = "MICE sensitivity",
  model = "Not fitted",
  adjustment = "Not applicable",
  model_family = "Not estimable with a reliable death time",
  n = outcome_available_n,
  events = observed_expired_by_28d_n,
  reference = "Profile 1",
  exposure = "Profile 2",
  effect_measure = "Not estimated",
  estimate = NA_real_,
  confidence_interval_lower = NA_real_,
  confidence_interval_upper = NA_real_,
  p_value = NA_real_,
  between_imputation_variance = NA_real_,
  fraction_missing_information = NA_real_,
  degrees_freedom = NA_real_,
  mean_aic = NA_real_,
  status = "NOT_ESTIMABLE",
  notes = paste(
    "eICU provides hospital discharge status and offset but no standalone death timestamp.",
    "The event count is an aggregate disposition-offset approximation only; no model was fitted,",
    "and this endpoint is not 28-day all-cause mortality."
  ),
  stringsAsFactors = FALSE
)
primary_output <- rbind(primary_results, not_estimable_row)

six_variable_complete_index <- complete.cases(local_data[primary_variables])
if (sum(six_variable_complete_index) != expected_complete_case_n) {
  stop("Complete-case denominator does not reconcile to Stage 10C-1.")
}
complete_case_assignments <- stage10c1$profile_assignment[
  six_variable_complete_index, , drop = FALSE
]
if (!all(apply(complete_case_assignments, 1L, function(x) length(unique(x)) == 1L))) {
  stop("A complete-case profile assignment varies across imputations.")
}
complete_case_profile <- complete_case_assignments[, 1L]
if (!identical(tabulate(complete_case_profile, nbins = 2L), c(500L, 241L))) {
  stop("Complete-case profile distribution does not reconcile to Stage 10C-1 QC.")
}
complete_case_model_index <- six_variable_complete_index & outcome_available
complete_case_model_n <- sum(complete_case_model_index)
complete_case_fits <- fit_logistic(
  stage10c1$profile_assignment[complete_case_model_index, 1L],
  row_index = which(complete_case_model_index)
)
complete_case_results <- do.call(rbind, lapply(
  seq_len(nrow(model_specification)),
  function(model_index) {
    fit <- complete_case_fits[[model_index]]
    standard_error <- sqrt(fit$variance)
    data.frame(
      section = "association",
      outcome = "Hospital mortality",
      analysis_set = "Six-variable complete case",
      model = model_specification$model[model_index],
      adjustment = model_specification$adjustment[model_index],
      model_family = "Logistic regression",
      n = complete_case_model_n,
      events = sum(hospital_mortality[complete_case_model_index]),
      reference = "Profile 1",
      exposure = "Profile 2",
      effect_measure = "Odds ratio",
      estimate = exp(fit$estimate),
      confidence_interval_lower = exp(fit$estimate - stats::qnorm(0.975) * standard_error),
      confidence_interval_upper = exp(fit$estimate + stats::qnorm(0.975) * standard_error),
      p_value = 2 * stats::pnorm(-abs(fit$estimate / standard_error)),
      between_imputation_variance = 0,
      fraction_missing_information = 0,
      degrees_freedom = Inf,
      mean_aic = fit$aic,
      status = "ESTIMATED",
      notes = sprintf(
        "Observed six-variable complete cases with known hospital status; %d record(s) with unknown status excluded from the 741 laboratory-complete stays; frozen Stage 10C-1 assignment and 0.5 threshold.",
        sum(six_variable_complete_index & !outcome_available)
      ),
      stringsAsFactors = FALSE
    )
  }
))

numeric_result_columns <- c(
  "estimate", "confidence_interval_lower", "confidence_interval_upper",
  "p_value", "between_imputation_variance",
  "fraction_missing_information", "degrees_freedom", "mean_aic"
)
primary_output[numeric_result_columns] <- lapply(
  primary_output[numeric_result_columns],
  function(x) ifelse(is.na(x), NA_real_, signif(x, 8L))
)
complete_case_results[numeric_result_columns] <- lapply(
  complete_case_results[numeric_result_columns],
  function(x) ifelse(is.na(x), NA_real_, signif(x, 8L))
)

if (
  nrow(primary_output) != 4L ||
  nrow(complete_case_results) != 3L ||
  sum(primary_output$status == "ESTIMATED") != 3L ||
  sum(primary_output$status == "NOT_ESTIMABLE") != 1L ||
  any(primary_output$estimate[primary_output$status == "ESTIMATED"] <= 0) ||
  any(complete_case_results$estimate <= 0) ||
  any(primary_output$confidence_interval_lower[primary_output$status == "ESTIMATED"] >
      primary_output$estimate[primary_output$status == "ESTIMATED"]) ||
  any(primary_output$confidence_interval_upper[primary_output$status == "ESTIMATED"] <
      primary_output$estimate[primary_output$status == "ESTIMATED"]) ||
  reliable_death_time_available
) {
  stop("Stage 10C-2 aggregate result validation failed.")
}

# Publication-ready forest plot. It displays only aggregate odds ratios.
plot_mi <- primary_output[primary_output$status == "ESTIMATED", ]
plot_cc <- complete_case_results
plot_order <- 3:1
plot_mi <- plot_mi[plot_order, ]
plot_cc <- plot_cc[plot_order, ]
base_y <- 1:3
y_mi <- base_y + 0.12
y_cc <- base_y - 0.12
all_lower <- c(plot_mi$confidence_interval_lower, plot_cc$confidence_interval_lower)
all_upper <- c(plot_mi$confidence_interval_upper, plot_cc$confidence_interval_upper)
x_min <- min(0.5, all_lower) * 0.88
x_max <- max(2.0, all_upper) * 1.12
tick_candidates <- c(0.5, 0.75, 1, 1.5, 2, 3, 4, 6)
ticks <- tick_candidates[tick_candidates >= x_min & tick_candidates <= x_max]
if (!1 %in% ticks) ticks <- sort(unique(c(ticks, 1)))

figure_tmp <- tempfile(
  "stage10c2_figure_", tmpdir = dirname(figure_path), fileext = ".pdf"
)
grDevices::pdf(figure_tmp, width = 9.2, height = 6.4, useDingbats = FALSE)
graphics::par(mar = c(7.0, 12.5, 5.6, 2.0), las = 1, family = "sans")
graphics::plot(
  NA_real_, NA_real_,
  xlim = log(c(x_min, x_max)), ylim = c(0.45, 3.65),
  xaxt = "n", yaxt = "n", xlab = "Odds ratio (log scale)", ylab = "",
  main = "eICU external validation of inflammatory-metabolic profiles"
)
graphics::mtext(
  "Hospital mortality: Profile 2 versus Profile 1",
  side = 3, line = 0.8, cex = 0.95, col = "#374151"
)
graphics::axis(1, at = log(ticks), labels = format(ticks, trim = TRUE))
graphics::axis(
  2, at = base_y,
  labels = rev(c(
    "Model 1: Profile only",
    "Model 2: + age + sex",
    "Model 3: + cancer phenotype"
  )),
  tick = FALSE, line = 0.2, cex.axis = 0.92
)
graphics::abline(v = log(1), col = "#6B7280", lty = 2, lwd = 1.2)
graphics::abline(h = c(1.5, 2.5), col = "#E5E7EB", lwd = 0.8)
draw_estimates <- function(data, y, color, pch) {
  graphics::segments(
    log(data$confidence_interval_lower), y,
    log(data$confidence_interval_upper), y,
    col = color, lwd = 2.2
  )
  graphics::points(log(data$estimate), y, col = color, bg = color, pch = pch, cex = 1.15)
}
draw_estimates(plot_mi, y_mi, "#0072B2", 16)
draw_estimates(plot_cc, y_cc, "#D55E00", 17)
graphics::legend(
  "topright",
  legend = c(
    sprintf("MICE primary (n=%s)", format(outcome_available_n, big.mark = ",")),
    sprintf("Complete case (n=%s)", format(complete_case_model_n, big.mark = ","))
  ),
  col = c("#0072B2", "#D55E00"), pch = c(16, 17),
  pt.cex = 1.1, bty = "n", cex = 0.9
)
graphics::mtext(
  "Frozen MIMIC transfer rule; no retraining, threshold change, LPA refit, or laboratory adjustment.",
  side = 1, line = 4.5, cex = 0.74, col = "#374151"
)
graphics::mtext(
  "28-day in-hospital mortality was not modeled because eICU lacks a standalone reliable death timestamp.",
  side = 1, line = 5.5, cex = 0.70, col = "#6B7280"
)
grDevices::dev.off()

result_tmp <- tempfile("stage10c2_results_", tmpdir = dirname(result_path), fileext = ".csv")
complete_case_tmp <- tempfile(
  "stage10c2_complete_case_", tmpdir = dirname(complete_case_path), fileext = ".csv"
)
utils::write.csv(primary_output, result_tmp, row.names = FALSE, na = "")
utils::write.csv(complete_case_results, complete_case_tmp, row.names = FALSE, na = "")

verified_primary <- utils::read.csv(result_tmp, check.names = FALSE, stringsAsFactors = FALSE)
verified_complete_case <- utils::read.csv(
  complete_case_tmp, check.names = FALSE, stringsAsFactors = FALSE
)
if (
  nrow(verified_primary) != 4L ||
  nrow(verified_complete_case) != 3L ||
  any(grepl(
    "patientunitstayid|subject_id|stay_id|wbc|platelet|lactate|glucose|creatinine|bicarbonate",
    names(verified_primary), ignore.case = TRUE
  )) ||
  any(grepl(
    "patientunitstayid|subject_id|stay_id|wbc|platelet|lactate|glucose|creatinine|bicarbonate",
    names(verified_complete_case), ignore.case = TRUE
  )) ||
  !file.exists(figure_tmp) || file.info(figure_tmp)$size <= 1000
) {
  unlink(c(result_tmp, complete_case_tmp, figure_tmp))
  stop("Stage 10C-2 temporary aggregate outputs failed privacy or structural verification.")
}

if (!file.rename(result_tmp, result_path)) {
  unlink(c(result_tmp, complete_case_tmp, figure_tmp))
  stop("Could not finalize Stage 10C-2 primary results.")
}
if (!file.rename(complete_case_tmp, complete_case_path)) {
  unlink(c(complete_case_tmp, figure_tmp))
  stop("Primary results finalized but complete-case output failed finalization.")
}
if (!file.rename(figure_tmp, figure_path)) {
  unlink(figure_tmp)
  stop("Stage 10C-2 tables finalized but the validation figure failed finalization.")
}

rm(database_input, hospital_mortality, local_data, stage10c1, transfer_model)

cat("STAGE10C2_COMPLETE\n")
cat(sprintf(
  "cohort_n=%d outcome_available_n=%d outcome_missing_n=%d hospital_deaths=%d complete_case_n=%d complete_case_deaths=%d\n",
  expected_n, outcome_available_n, outcome_missing_n, hospital_death_n, complete_case_model_n,
  complete_case_results$events[1L]
))
for (model_index in seq_len(3L)) {
  cat(sprintf(
    "%s MICE_OR=%.4f CI=%.4f-%.4f CC_OR=%.4f CI=%.4f-%.4f\n",
    model_specification$model[model_index],
    primary_results$estimate[model_index],
    primary_results$confidence_interval_lower[model_index],
    primary_results$confidence_interval_upper[model_index],
    complete_case_results$estimate[model_index],
    complete_case_results$confidence_interval_lower[model_index],
    complete_case_results$confidence_interval_upper[model_index]
  ))
}
cat(sprintf(
  "observed_expired_by_28d_n=%d reliable_death_time=FALSE sensitivity_model_run=FALSE\n",
  observed_expired_by_28d_n
))
cat("classifier_retrained=FALSE threshold_changed=FALSE lpa_refitted=FALSE patient_level_export=FALSE\n")
