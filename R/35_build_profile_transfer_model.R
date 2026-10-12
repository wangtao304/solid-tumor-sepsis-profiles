options(stringsAsFactors = FALSE, warn = 1)

if (!requireNamespace("mice", quietly = TRUE)) {
  stop("Required package 'mice' is not installed; no installation is authorized.")
}

stage_seed <- 20260723L
cv_folds <- 10L
bootstrap_replicates <- 1000L
classification_threshold <- 0.5

mice_path <- file.path(
  "data", "processed",
  "inflammatory_metabolic_mice_solid_m10_maxit20_seed20260721.rds"
)
lpa_path <- file.path(
  "data", "processed",
  "inflammatory_metabolic_lpa_solid_2to6_mi10_seed20260722.rds"
)
model_path <- file.path(
  "data", "processed",
  "inflammatory_metabolic_profile_transfer_logistic_mi10_seed20260723.rds"
)
output_path <- file.path(
  "outputs", "tables", "35_profile_transfer_model_performance.csv"
)

for (path in c(mice_path, lpa_path)) {
  if (!file.exists(path)) stop("Missing required Stage 10B input: ", path)
}
dir.create(dirname(model_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

variables <- c(
  "wbc", "platelet_count", "lactate",
  "glucose", "creatinine", "bicarbonate"
)
coefficient_labels <- c(
  "(Intercept)", "WBC", "platelet", "log1p(lactate)",
  "glucose", "log1p(creatinine)", "bicarbonate"
)
prohibited_names <- c(
  "mortality", "death", "survival", "subject_id", "stay_id",
  "patientunitstayid", "outcome"
)

mice_object <- readRDS(mice_path)
lpa_object <- readRDS(lpa_path)
if (
  !inherits(mice_object, "mids") ||
  mice_object$m != 10L ||
  nrow(mice_object$data) != 4380L ||
  !all(variables %in% names(mice_object$data)) ||
  !identical(lpa_object$configuration$variables, variables) ||
  !identical(lpa_object$configuration$lactate_transform, "log1p") ||
  !identical(lpa_object$configuration$creatinine_transform, "log1p")
) {
  stop("Stage 10B input objects do not match the frozen six-variable specification.")
}
if (any(vapply(prohibited_names, function(x) {
  any(grepl(x, names(mice_object$data), ignore.case = TRUE))
}, logical(1L)))) {
  stop("A prohibited identifier or outcome is present in the classifier input.")
}

completed_datasets <- mice::complete(mice_object, action = "all")
if (length(completed_datasets) != 10L) {
  stop("Expected exactly ten completed MIMIC datasets.")
}

# Reconstruct the frozen Stage 9A consensus assignment; no LPA is fitted.
assignment_votes <- matrix(NA_integer_, nrow(mice_object$data), 10L)
mean_aligned_posterior <- matrix(0, nrow(mice_object$data), 2L)
for (imputation in seq_len(10L)) {
  key <- sprintf("imp%02d_k2", imputation)
  fit_record <- lpa_object$fit_records[[key]]
  alignment <- lpa_object$alignments[[key]]
  if (
    is.null(fit_record) || is.null(alignment) ||
    length(fit_record$model$classification) != nrow(mice_object$data) ||
    !identical(dim(fit_record$model$posterior), c(nrow(mice_object$data), 2L)) ||
    length(alignment) != 2L
  ) {
    stop("Incomplete frozen two-profile record: ", key)
  }
  assignment_votes[, imputation] <- match(
    fit_record$model$classification,
    alignment
  )
  mean_aligned_posterior <- mean_aligned_posterior +
    fit_record$model$posterior[, alignment, drop = FALSE] / 10
}
votes_profile_1 <- rowSums(assignment_votes == 1L)
votes_profile_2 <- rowSums(assignment_votes == 2L)
consensus_profile <- ifelse(
  votes_profile_1 > votes_profile_2,
  1L,
  ifelse(
    votes_profile_2 > votes_profile_1,
    2L,
    max.col(mean_aligned_posterior, ties.method = "first")
  )
)
if (
  anyNA(consensus_profile) ||
  !identical(tabulate(consensus_profile, nbins = 2L), c(3197L, 1183L)) ||
  sum(votes_profile_1 == votes_profile_2) != 26L
) {
  stop("Frozen consensus profile reconstruction failed.")
}
y <- as.integer(consensus_profile == 2L)

transform_values <- function(data, centers = NULL, scales = NULL) {
  x <- as.data.frame(data[, variables, drop = FALSE])
  x[] <- lapply(x, as.numeric)
  if (anyNA(x) || any(!is.finite(as.matrix(x)))) {
    stop("Classifier predictors must be complete and finite.")
  }
  if (any(x$lactate < 0) || any(x$creatinine < 0)) {
    stop("log1p inputs must be nonnegative.")
  }
  x$lactate <- log1p(x$lactate)
  x$creatinine <- log1p(x$creatinine)
  if (is.null(centers)) centers <- vapply(x, mean, numeric(1L))
  if (is.null(scales)) scales <- vapply(x, stats::sd, numeric(1L))
  if (
    !identical(names(centers), variables) ||
    !identical(names(scales), variables) ||
    any(!is.finite(centers)) || any(!is.finite(scales)) || any(scales <= 0)
  ) {
    stop("Invalid standardization parameters.")
  }
  z <- sweep(as.matrix(x), 2L, centers, "-")
  z <- sweep(z, 2L, scales, "/")
  colnames(z) <- variables
  list(values = z, centers = centers, scales = scales)
}

fit_component <- function(x, outcome) {
  fit_data <- data.frame(profile_2 = outcome, x, check.names = FALSE)
  fit <- stats::glm(profile_2 ~ ., data = fit_data, family = stats::binomial())
  if (!isTRUE(fit$converged) || any(!is.finite(stats::coef(fit)))) {
    stop("A logistic-regression component failed to converge.")
  }
  fit
}

# Stratified ten-fold cross-validation. The same folds are used across all
# completed datasets so predictions can be averaged into one transport rule.
set.seed(stage_seed)
fold_id <- integer(length(y))
for (class_value in 0:1) {
  index <- which(y == class_value)
  fold_id[index] <- sample(rep(seq_len(cv_folds), length.out = length(index)))
}
if (any(tabulate(fold_id, nbins = cv_folds) == 0L)) {
  stop("Cross-validation fold generation failed.")
}

oof_probability <- rep(NA_real_, length(y))
for (fold in seq_len(cv_folds)) {
  training <- fold_id != fold
  testing <- !training
  component_probability <- matrix(NA_real_, sum(testing), 10L)
  for (imputation in seq_len(10L)) {
    transformed_training <- transform_values(
      completed_datasets[[imputation]][training, , drop = FALSE]
    )
    transformed_testing <- transform_values(
      completed_datasets[[imputation]][testing, , drop = FALSE],
      centers = transformed_training$centers,
      scales = transformed_training$scales
    )
    fit <- fit_component(transformed_training$values, y[training])
    component_probability[, imputation] <- stats::predict(
      fit,
      newdata = as.data.frame(transformed_testing$values),
      type = "response"
    )
  }
  oof_probability[testing] <- rowMeans(component_probability)
}
if (
  anyNA(oof_probability) ||
  any(!is.finite(oof_probability)) ||
  any(oof_probability <= 0 | oof_probability >= 1)
) {
  stop("Out-of-fold probability construction failed.")
}

auc_rank <- function(outcome, probability) {
  n1 <- sum(outcome == 1L)
  n0 <- sum(outcome == 0L)
  ranks <- rank(probability, ties.method = "average")
  (sum(ranks[outcome == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

classification_metrics <- function(outcome, probability, threshold = 0.5) {
  prediction <- as.integer(probability >= threshold)
  tp <- sum(prediction == 1L & outcome == 1L)
  tn <- sum(prediction == 0L & outcome == 0L)
  fp <- sum(prediction == 1L & outcome == 0L)
  fn <- sum(prediction == 0L & outcome == 1L)
  agreement <- (tp + tn) / length(outcome)
  expected_agreement <-
    ((tp + fp) * (tp + fn) + (fn + tn) * (fp + tn)) / length(outcome)^2
  c(
    auroc = auc_rank(outcome, probability),
    brier_score = mean((probability - outcome)^2),
    classification_agreement = agreement,
    cohen_kappa = (agreement - expected_agreement) / (1 - expected_agreement),
    sensitivity_profile_2 = tp / (tp + fn),
    specificity_profile_1 = tn / (tn + fp),
    positive_predictive_value = tp / (tp + fp),
    negative_predictive_value = tn / (tn + fn),
    predicted_profile_2_proportion = mean(prediction == 1L)
  )
}

performance <- classification_metrics(y, oof_probability, classification_threshold)

# Stratified bootstrap intervals for nonparametric performance/agreement metrics.
set.seed(stage_seed + 1L)
positive_index <- which(y == 1L)
negative_index <- which(y == 0L)
bootstrap_metrics <- matrix(
  NA_real_,
  nrow = bootstrap_replicates,
  ncol = length(performance),
  dimnames = list(NULL, names(performance))
)
for (replicate in seq_len(bootstrap_replicates)) {
  index <- c(
    sample(positive_index, length(positive_index), replace = TRUE),
    sample(negative_index, length(negative_index), replace = TRUE)
  )
  bootstrap_metrics[replicate, ] <- classification_metrics(
    y[index],
    oof_probability[index],
    classification_threshold
  )
}
bootstrap_lower <- apply(bootstrap_metrics, 2L, stats::quantile, probs = 0.025)
bootstrap_upper <- apply(bootstrap_metrics, 2L, stats::quantile, probs = 0.975)

clipped_probability <- pmin(pmax(oof_probability, 1e-8), 1 - 1e-8)
linear_predictor <- stats::qlogis(clipped_probability)
calibration_intercept_fit <- stats::glm(
  y ~ 1L + offset(linear_predictor),
  family = stats::binomial()
)
calibration_slope_fit <- stats::glm(y ~ linear_predictor, family = stats::binomial())
calibration_intercept <- stats::coef(calibration_intercept_fit)[1L]
calibration_intercept_se <- sqrt(stats::vcov(calibration_intercept_fit)[1L, 1L])
calibration_slope <- stats::coef(calibration_slope_fit)[2L]
calibration_slope_se <- sqrt(stats::vcov(calibration_slope_fit)[2L, 2L])

# Fit the final ten-component ensemble on all MIMIC rows. Only aggregate model
# parameters and preprocessing constants are saved; no row-level labels or
# probabilities are stored in the model object.
coefficient_matrix <- matrix(
  NA_real_, 7L, 10L,
  dimnames = list(coefficient_labels, sprintf("imputation_%02d", seq_len(10L)))
)
centers_matrix <- matrix(
  NA_real_, 6L, 10L,
  dimnames = list(variables, sprintf("imputation_%02d", seq_len(10L)))
)
scales_matrix <- centers_matrix
for (imputation in seq_len(10L)) {
  transformed <- transform_values(completed_datasets[[imputation]])
  stored_preprocessing <- lpa_object$preprocessing[[imputation]]
  if (
    max(abs(transformed$centers - stored_preprocessing$centers)) > 1e-10 ||
    max(abs(transformed$scales - stored_preprocessing$scales)) > 1e-10
  ) {
    stop("Classifier preprocessing does not reproduce frozen LPA preprocessing.")
  }
  fit <- fit_component(transformed$values, y)
  coefficient_matrix[, imputation] <- stats::coef(fit)
  centers_matrix[, imputation] <- transformed$centers
  scales_matrix[, imputation] <- transformed$scales
}

coefficient_rows <- data.frame(
  section = "coefficient",
  metric = rownames(coefficient_matrix),
  estimate = rowMeans(coefficient_matrix),
  lower_bound = apply(coefficient_matrix, 1L, min),
  upper_bound = apply(coefficient_matrix, 1L, max),
  bound_type = "Range across 10 final component models",
  evaluation = "Full MIMIC cohort; final deployable ensemble",
  threshold = NA_real_,
  n = length(y),
  profile_2_n = sum(y),
  notes = c(
    "Mean intercept across component logistic models",
    rep("Log-odds coefficient per 1-SD increase after frozen preprocessing", 6L)
  ),
  stringsAsFactors = FALSE
)

performance_labels <- c(
  auroc = "AUROC",
  brier_score = "Brier score",
  classification_agreement = "Classification agreement",
  cohen_kappa = "Cohen kappa",
  sensitivity_profile_2 = "Sensitivity for Profile 2",
  specificity_profile_1 = "Specificity for Profile 1",
  positive_predictive_value = "Positive predictive value",
  negative_predictive_value = "Negative predictive value",
  predicted_profile_2_proportion = "Predicted Profile 2 proportion"
)
performance_sections <- c(
  auroc = "discrimination",
  brier_score = "calibration",
  classification_agreement = "classification_agreement",
  cohen_kappa = "classification_agreement",
  sensitivity_profile_2 = "classification_agreement",
  specificity_profile_1 = "classification_agreement",
  positive_predictive_value = "classification_agreement",
  negative_predictive_value = "classification_agreement",
  predicted_profile_2_proportion = "classification_agreement"
)
performance_rows <- data.frame(
  section = unname(performance_sections[names(performance)]),
  metric = unname(performance_labels[names(performance)]),
  estimate = as.numeric(performance),
  lower_bound = as.numeric(bootstrap_lower[names(performance)]),
  upper_bound = as.numeric(bootstrap_upper[names(performance)]),
  bound_type = "95% stratified bootstrap interval",
  evaluation = "MIMIC stratified 10-fold out-of-fold ensemble predictions",
  threshold = c(NA_real_, NA_real_, rep(classification_threshold, 7L)),
  n = length(y),
  profile_2_n = sum(y),
  notes = "Profile 1 reference; Profile 2 exposure; frozen profile assignment",
  stringsAsFactors = FALSE
)
calibration_rows <- data.frame(
  section = "calibration",
  metric = c("Calibration intercept", "Calibration slope"),
  estimate = c(calibration_intercept, calibration_slope),
  lower_bound = c(
    calibration_intercept - stats::qnorm(0.975) * calibration_intercept_se,
    calibration_slope - stats::qnorm(0.975) * calibration_slope_se
  ),
  upper_bound = c(
    calibration_intercept + stats::qnorm(0.975) * calibration_intercept_se,
    calibration_slope + stats::qnorm(0.975) * calibration_slope_se
  ),
  bound_type = "95% Wald confidence interval",
  evaluation = "MIMIC stratified 10-fold out-of-fold ensemble predictions",
  threshold = NA_real_,
  n = length(y),
  profile_2_n = sum(y),
  notes = c("Ideal value 0", "Ideal value 1"),
  stringsAsFactors = FALSE
)

results <- rbind(coefficient_rows, performance_rows, calibration_rows)
numeric_columns <- c("estimate", "lower_bound", "upper_bound", "threshold")
results[numeric_columns] <- lapply(results[numeric_columns], function(x) signif(x, 8L))
if (
  nrow(results) != 18L ||
  anyNA(results[, setdiff(names(results), "threshold")]) ||
  any(!is.finite(unlist(results[c("estimate", "lower_bound", "upper_bound")])) ) ||
  any(results$lower_bound > results$estimate) ||
  any(results$upper_bound < results$estimate)
) {
  stop("Stage 10B aggregate result validation failed.")
}

transfer_model <- list(
  stage = "Stage 10B",
  created_date = as.character(Sys.Date()),
  model_type = "Ten-component logistic-regression probability ensemble",
  target = "Frozen MIMIC consensus Profile 2 versus Profile 1",
  reference_profile = "Profile 1: lower inflammatory-metabolic stress",
  exposure_profile = "Profile 2: higher inflammatory-metabolic stress",
  variables = variables,
  transformations = c(lactate = "log1p", creatinine = "log1p"),
  standardization = "Component-specific MIMIC mean and standard deviation",
  probability_rule = "Average the 10 component probabilities",
  classification_threshold = classification_threshold,
  classification_rule = "Probability >=0.5 assigns Profile 2; otherwise Profile 1",
  coefficients = coefficient_matrix,
  centers = centers_matrix,
  scales = scales_matrix,
  cross_validation = list(
    folds = cv_folds,
    stratified_by_profile = TRUE,
    seed = stage_seed,
    performance = performance,
    calibration_intercept = unname(calibration_intercept),
    calibration_slope = unname(calibration_slope)
  ),
  training_counts = c(profile_1 = sum(y == 0L), profile_2 = sum(y == 1L)),
  training_source = "MIMIC only; approved Stage 7B MICE and Stage 8A LPA objects",
  eicu_data_used_for_training = FALSE,
  mortality_or_survival_used = FALSE,
  lpa_refitted = FALSE,
  patient_level_labels_or_probabilities_stored = FALSE,
  input_md5 = unname(tools::md5sum(c(mice_path, lpa_path))),
  r_version = R.version.string
)
if (any(vapply(prohibited_names, function(x) {
  grepl(x, paste(names(transfer_model), collapse = "|"), ignore.case = TRUE)
}, logical(1L)))) {
  # Explicit negative audit fields contain mortality/survival terms by design;
  # confirm their values are FALSE and save no corresponding data.
  if (
    !identical(transfer_model$mortality_or_survival_used, FALSE) ||
    !identical(transfer_model$patient_level_labels_or_probabilities_stored, FALSE)
  ) stop("Prohibited-data audit flag failed.")
}

model_tmp <- tempfile("stage10b_model_", tmpdir = dirname(model_path), fileext = ".rds")
csv_tmp <- tempfile("stage10b_results_", tmpdir = dirname(output_path), fileext = ".csv")
saveRDS(transfer_model, model_tmp, version = 3L)
utils::write.csv(results, csv_tmp, row.names = FALSE, na = "")

verified_model <- readRDS(model_tmp)
verified_results <- utils::read.csv(csv_tmp, check.names = FALSE)
if (
  !identical(dim(verified_model$coefficients), c(7L, 10L)) ||
  !identical(dim(verified_model$centers), c(6L, 10L)) ||
  !identical(dim(verified_model$scales), c(6L, 10L)) ||
  !identical(verified_model$training_counts, c(profile_1 = 3197L, profile_2 = 1183L)) ||
  nrow(verified_results) != 18L ||
  any(grepl("subject_id|stay_id|patientunitstayid", names(verified_results), ignore.case = TRUE))
) {
  unlink(c(model_tmp, csv_tmp))
  stop("Stage 10B finalized-object verification failed.")
}
if (!file.copy(model_tmp, model_path, overwrite = TRUE)) {
  unlink(c(model_tmp, csv_tmp))
  stop("Could not finalize the controlled transfer-model object.")
}
if (!file.copy(csv_tmp, output_path, overwrite = TRUE)) {
  unlink(c(model_tmp, csv_tmp))
  stop("Could not finalize the aggregate Stage 10B table.")
}
unlink(c(model_tmp, csv_tmp))

cat(sprintf(
  paste0(
    "STAGE10B_COMPLETE n=%d profile2=%d auc=%.6f calibration_intercept=%.6f ",
    "calibration_slope=%.6f agreement=%.6f kappa=%.6f\n"
  ),
  length(y), sum(y), performance[["auroc"]], calibration_intercept,
  calibration_slope, performance[["classification_agreement"]],
  performance[["cohen_kappa"]]
))
