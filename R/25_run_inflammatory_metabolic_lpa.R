#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

if (!requireNamespace("mice", quietly = TRUE)) {
  stop("Required package 'mice' is not installed; no installation is authorized.")
}

input_path <- "data/processed/inflammatory_metabolic_mice_solid_m10_maxit20_seed20260721.rds"
model_path <- "data/processed/inflammatory_metabolic_lpa_solid_2to6_mi10_seed20260722.rds"
fit_path <- "outputs/tables/25_lpa_model_fit.csv"
size_path <- "outputs/tables/26_lpa_profile_size.csv"
checkpoint_path <- "data/interim/stage8a_lpa_fit_checkpoint_mi10_seed20260722.rds"

primary_variables <- c(
  "wbc",
  "platelet_count",
  "lactate",
  "glucose",
  "creatinine",
  "bicarbonate"
)
inflammatory_variables <- c("wbc", "platelet_count")
metabolic_variables <- c("lactate", "glucose", "creatinine", "bicarbonate")
candidate_profiles <- 2:6
imputation_n <- 10L
random_starts <- 30L
maximum_iterations <- 1000L
convergence_tolerance <- 1e-8
variance_floor <- 0.0025
stage8a_seed <- 20260722L

creatinine_skewness_threshold <- 1.0
minimum_profile_proportion <- 0.05
maximum_centroid_deviation <- 1.0
maximum_profile_range_pp <- 10
minimum_centroid_separation <- 0.5
minimum_domain_contrast <- 0.5

adjusted_skewness <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 3L) return(NA_real_)
  standard_deviation <- stats::sd(x)
  if (!is.finite(standard_deviation) || standard_deviation <= 0) return(0)
  standardized <- (x - mean(x)) / standard_deviation
  n * sum(standardized^3) / ((n - 1) * (n - 2))
}

log_sum_exp_rows <- function(log_matrix) {
  row_maximum <- log_matrix[, 1L]
  if (ncol(log_matrix) > 1L) {
    for (column_index in 2:ncol(log_matrix)) {
      row_maximum <- pmax(row_maximum, log_matrix[, column_index])
    }
  }
  row_maximum + log(rowSums(exp(log_matrix - row_maximum)))
}

expectation_step <- function(x, proportions, means, variances) {
  n <- nrow(x)
  profiles <- length(proportions)
  log_probabilities <- matrix(NA_real_, nrow = n, ncol = profiles)

  for (profile_index in seq_len(profiles)) {
    profile_variances <- variances[profile_index, ]
    residuals <- sweep(x, 2L, means[profile_index, ], FUN = "-")
    quadratic <- rowSums(sweep(residuals^2, 2L, profile_variances, FUN = "/"))
    normalizing_constant <- sum(log(2 * pi * profile_variances))
    log_probabilities[, profile_index] <-
      log(proportions[profile_index]) - 0.5 * (normalizing_constant + quadratic)
  }

  log_denominator <- log_sum_exp_rows(log_probabilities)
  posterior <- exp(log_probabilities - log_denominator)
  list(log_likelihood = sum(log_denominator), posterior = posterior)
}

initialize_parameters <- function(x, profiles) {
  warning_messages <- character()
  clustering <- withCallingHandlers(
    stats::kmeans(
      x,
      centers = profiles,
      iter.max = 100L,
      nstart = 1L,
      algorithm = "Lloyd"
    )$cluster,
    warning = function(warning) {
      warning_messages <<- c(warning_messages, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )

  counts <- tabulate(clustering, nbins = profiles)
  if (any(counts < 2L)) stop("Initialization produced a profile with fewer than two rows.")

  means <- rowsum(x, clustering, reorder = TRUE) / counts
  variances <- matrix(
    NA_real_,
    nrow = profiles,
    ncol = ncol(x),
    dimnames = list(NULL, colnames(x))
  )
  for (profile_index in seq_len(profiles)) {
    residuals <- sweep(
      x[clustering == profile_index, , drop = FALSE],
      2L,
      means[profile_index, ],
      FUN = "-"
    )
    variances[profile_index, ] <- colMeans(residuals^2)
  }

  list(
    proportions = counts / nrow(x),
    means = means,
    variances = pmax(variances, variance_floor),
    initialization_warning_n = length(warning_messages),
    initialization_warnings = unique(warning_messages)
  )
}

fit_single_start <- function(x, profiles) {
  parameters <- initialize_parameters(x, profiles)
  previous_log_likelihood <- -Inf
  converged <- FALSE
  iteration <- 0L
  log_likelihood_decrease_n <- 0L

  for (iteration in seq_len(maximum_iterations)) {
    expectation <- expectation_step(
      x,
      parameters$proportions,
      parameters$means,
      parameters$variances
    )
    current_log_likelihood <- expectation$log_likelihood

    if (
      is.finite(previous_log_likelihood) &&
      current_log_likelihood < previous_log_likelihood -
        convergence_tolerance * (1 + abs(previous_log_likelihood))
    ) {
      log_likelihood_decrease_n <- log_likelihood_decrease_n + 1L
    }

    if (
      is.finite(previous_log_likelihood) &&
      abs(current_log_likelihood - previous_log_likelihood) <=
        convergence_tolerance * (1 + abs(previous_log_likelihood))
    ) {
      converged <- TRUE
      break
    }

    previous_log_likelihood <- current_log_likelihood
    posterior <- expectation$posterior
    effective_n <- colSums(posterior)
    if (any(effective_n < 2)) {
      return(list(converged = FALSE, log_likelihood = -Inf))
    }

    parameters$proportions <- effective_n / nrow(x)
    parameters$means <- sweep(
      crossprod(posterior, x),
      1L,
      effective_n,
      FUN = "/"
    )

    updated_variances <- matrix(
      NA_real_,
      nrow = profiles,
      ncol = ncol(x),
      dimnames = list(NULL, colnames(x))
    )
    for (profile_index in seq_len(profiles)) {
      residuals <- sweep(x, 2L, parameters$means[profile_index, ], FUN = "-")
      updated_variances[profile_index, ] <-
        colSums(residuals^2 * posterior[, profile_index]) / effective_n[profile_index]
    }
    parameters$variances <- pmax(updated_variances, variance_floor)
  }

  final_expectation <- expectation_step(
    x,
    parameters$proportions,
    parameters$means,
    parameters$variances
  )
  posterior <- final_expectation$posterior
  classification <- max.col(posterior, ties.method = "first")
  classified_sizes <- tabulate(classification, nbins = profiles)

  list(
    converged = converged && log_likelihood_decrease_n == 0L && all(classified_sizes > 0L),
    iterations = iteration,
    log_likelihood = final_expectation$log_likelihood,
    proportions = parameters$proportions,
    means = parameters$means,
    variances = parameters$variances,
    posterior = posterior,
    classification = classification,
    boundary_variance_n = sum(parameters$variances <= variance_floor * (1 + 1e-8)),
    log_likelihood_decrease_n = log_likelihood_decrease_n,
    initialization_warning_n = parameters$initialization_warning_n,
    initialization_warnings = parameters$initialization_warnings
  )
}

fit_multistart <- function(x, profiles, seed) {
  set.seed(seed)
  start_fits <- vector("list", random_starts)
  start_log_likelihoods <- rep(-Inf, random_starts)
  start_converged <- rep(FALSE, random_starts)
  start_boundary_n <- rep(NA_integer_, random_starts)
  start_warning_n <- rep(0L, random_starts)

  for (start_index in seq_len(random_starts)) {
    candidate <- tryCatch(
      fit_single_start(x, profiles),
      error = function(error) {
        list(
          converged = FALSE,
          log_likelihood = -Inf,
          error = conditionMessage(error),
          initialization_warning_n = 0L
        )
      }
    )
    start_fits[[start_index]] <- candidate
    start_log_likelihoods[start_index] <- candidate$log_likelihood
    start_converged[start_index] <- isTRUE(candidate$converged)
    if (!is.null(candidate$boundary_variance_n)) {
      start_boundary_n[start_index] <- candidate$boundary_variance_n
    }
    start_warning_n[start_index] <- candidate$initialization_warning_n
  }

  converged_indices <- which(start_converged & is.finite(start_log_likelihoods))
  if (length(converged_indices) == 0L) {
    stop("No random start converged for a candidate LPA model.")
  }

  best_index <- converged_indices[which.max(start_log_likelihoods[converged_indices])]
  best_fit <- start_fits[[best_index]]
  replication_tolerance <- 1e-6 * (1 + abs(best_fit$log_likelihood))
  best_fit$best_start <- best_index
  best_fit$best_log_likelihood_replications <- sum(
    start_converged &
      abs(start_log_likelihoods - best_fit$log_likelihood) <= replication_tolerance
  )
  best_fit$start_diagnostics <- data.frame(
    start = seq_len(random_starts),
    converged = start_converged,
    log_likelihood = start_log_likelihoods,
    boundary_variance_n = start_boundary_n,
    initialization_warning_n = start_warning_n
  )
  best_fit$initialization_warning_starts <- sum(start_warning_n > 0L)
  best_fit
}

normalized_entropy <- function(posterior) {
  profiles <- ncol(posterior)
  if (profiles <= 1L) return(1)
  clipped <- pmax(posterior, .Machine$double.eps)
  1 + sum(clipped * log(clipped)) / (nrow(posterior) * log(profiles))
}

minimum_centroid_distance <- function(means) {
  distance_matrix <- as.matrix(stats::dist(means))
  min(distance_matrix[upper.tri(distance_matrix)])
}

all_permutations <- function(values) {
  if (length(values) == 1L) return(matrix(values, nrow = 1L))
  do.call(
    rbind,
    lapply(seq_along(values), function(index) {
      cbind(values[index], all_permutations(values[-index]))
    })
  )
}

minimum_alignment <- function(reference_means, target_means, permutations) {
  costs <- apply(permutations, 1L, function(permutation) {
    sum((target_means[permutation, , drop = FALSE] - reference_means)^2)
  })
  permutations[which.min(costs), ]
}

distribution_summary <- function(x, variable, imputation, transformation) {
  quantiles <- stats::quantile(x, probs = c(0.25, 0.5, 0.75), names = FALSE)
  data.frame(
    variable = variable,
    imputation = imputation,
    transformation = transformation,
    n = length(x),
    minimum = min(x),
    q1 = quantiles[1L],
    median = quantiles[2L],
    q3 = quantiles[3L],
    maximum = max(x),
    mean = mean(x),
    standard_deviation = stats::sd(x),
    skewness = adjusted_skewness(x)
  )
}

preprocess_dataset <- function(data, creatinine_log1p) {
  model_data <- data[, primary_variables, drop = FALSE]
  if (anyNA(model_data)) stop("A completed dataset contains missing primary values.")
  if (any(!vapply(model_data, is.numeric, logical(1)))) {
    stop("An LPA variable is not numeric.")
  }
  if (any(!is.finite(as.matrix(model_data)))) {
    stop("An LPA variable contains a non-finite value.")
  }
  if (any(model_data$lactate < 0)) stop("Lactate contains a negative value.")
  if (creatinine_log1p && any(model_data$creatinine < 0)) {
    stop("Creatinine contains a negative value.")
  }

  model_data$lactate <- log1p(model_data$lactate)
  if (creatinine_log1p) model_data$creatinine <- log1p(model_data$creatinine)

  centers <- vapply(model_data, mean, numeric(1))
  scales <- vapply(model_data, stats::sd, numeric(1))
  if (any(!is.finite(scales)) || any(scales <= 0)) {
    stop("An LPA variable has an invalid scale.")
  }

  standardized <- sweep(
    sweep(as.matrix(model_data), 2L, centers, FUN = "-"),
    2L,
    scales,
    FUN = "/"
  )
  list(data = standardized, centers = centers, scales = scales)
}

# Internal synthetic recovery test for the heteroskedastic diagonal Gaussian EM.
set.seed(8021L)
self_test_x <- rbind(
  cbind(stats::rnorm(250, -1.1, 0.4), stats::rnorm(250, -0.7, 0.5)),
  cbind(stats::rnorm(250, 1.1, 0.9), stats::rnorm(250, 0.7, 1.0))
)
colnames(self_test_x) <- c("x1", "x2")
self_test_fit <- fit_multistart(self_test_x, profiles = 2L, seed = 8021L)
if (
  !self_test_fit$converged ||
  !is.finite(self_test_fit$log_likelihood) ||
  length(unique(self_test_fit$classification)) != 2L ||
  any(!is.finite(self_test_fit$variances))
) {
  stop("Stage 8A LPA EM self-test failed.")
}
rm(self_test_x, self_test_fit)

if (identical(Sys.getenv("STAGE8A_SELF_TEST_ONLY"), "1")) {
  cat("STAGE8A_LPA_EM_SELF_TEST_OK\n")
  quit(save = "no", status = 0L)
}

if (!file.exists(input_path)) stop("The approved Stage 7B MICE object is missing.")
imputation_object <- readRDS(input_path)
if (!inherits(imputation_object, "mids") || imputation_object$m != imputation_n) {
  stop("The Stage 7B input is not the expected 10-imputation mids object.")
}
if (!all(primary_variables %in% names(imputation_object$data))) {
  stop("A frozen inflammatory-metabolic variable is missing from the MICE object.")
}
if (any(grepl(
  "subject_id|stay_id|hadm_id|mortality|death|survival|outcome|los",
  names(imputation_object$data),
  ignore.case = TRUE
))) {
  stop("The MICE input unexpectedly contains an identifier, outcome, survival, or LOS field.")
}

completed_datasets <- lapply(
  seq_len(imputation_n),
  function(imputation_index) mice::complete(imputation_object, action = imputation_index)
)
if (any(vapply(completed_datasets, nrow, integer(1)) != nrow(imputation_object$data))) {
  stop("Completed datasets do not preserve the expected row count.")
}
if (any(vapply(completed_datasets, function(data) {
  anyNA(data[, primary_variables, drop = FALSE])
}, logical(1)))) {
  stop("A completed dataset retains missing primary values.")
}

distribution_rows <- list()
for (imputation_index in seq_len(imputation_n)) {
  completed_data <- completed_datasets[[imputation_index]]
  for (variable in primary_variables) {
    distribution_rows[[length(distribution_rows) + 1L]] <- distribution_summary(
      completed_data[[variable]],
      variable,
      imputation_index,
      "raw"
    )
  }
  distribution_rows[[length(distribution_rows) + 1L]] <- distribution_summary(
    log1p(completed_data$lactate),
    "lactate",
    imputation_index,
    "log1p"
  )
  if (any(completed_data$creatinine < 0)) {
    stop("Creatinine contains a negative value and cannot be evaluated with log1p.")
  }
  distribution_rows[[length(distribution_rows) + 1L]] <- distribution_summary(
    log1p(completed_data$creatinine),
    "creatinine",
    imputation_index,
    "log1p"
  )
}
distribution_audit <- do.call(rbind, distribution_rows)
row.names(distribution_audit) <- NULL

creatinine_raw_skewness <- distribution_audit$skewness[
  distribution_audit$variable == "creatinine" &
    distribution_audit$transformation == "raw"
]
creatinine_log1p_skewness <- distribution_audit$skewness[
  distribution_audit$variable == "creatinine" &
    distribution_audit$transformation == "log1p"
]
creatinine_log1p <-
  stats::median(creatinine_raw_skewness) > creatinine_skewness_threshold &&
  stats::median(abs(creatinine_log1p_skewness)) <
    stats::median(abs(creatinine_raw_skewness))
creatinine_transform <- if (creatinine_log1p) "log1p" else "none"

cat(sprintf(
  paste0(
    "CREATININE_TRANSFORM=%s median_raw_skewness=%.6f ",
    "median_log1p_skewness=%.6f threshold=%.1f\n"
  ),
  creatinine_transform,
  stats::median(creatinine_raw_skewness),
  stats::median(creatinine_log1p_skewness),
  creatinine_skewness_threshold
))

if (identical(Sys.getenv("STAGE8A_TRANSFORM_AUDIT_ONLY"), "1")) {
  cat("STAGE8A_TRANSFORM_AUDIT_OK\n")
  quit(save = "no", status = 0L)
}

if (identical(Sys.getenv("STAGE8A_PROFILE_ONE"), "1")) {
  profile_data <- preprocess_dataset(completed_datasets[[1L]], creatinine_log1p)$data
  profile_start <- proc.time()[["elapsed"]]
  profile_fit <- fit_single_start(profile_data, profiles = 2L)
  profile_elapsed <- proc.time()[["elapsed"]] - profile_start
  cat(sprintf(
    "STAGE8A_PROFILE_ONE iterations=%d converged=%s boundary=%d elapsed_seconds=%.3f\n",
    profile_fit$iterations,
    profile_fit$converged,
    profile_fit$boundary_variance_n,
    profile_elapsed
  ))
  quit(save = "no", status = 0L)
}

if (file.exists(model_path) || file.exists(fit_path) || file.exists(size_path)) {
  stop("A Stage 8A output already exists; refusing to overwrite it.")
}

fit_records <- list()
fit_metrics <- list()
preprocessing_records <- vector("list", imputation_n)

checkpoint_signature <- list(
  input_path = input_path,
  input_md5 = unname(tools::md5sum(input_path)),
  variables = primary_variables,
  candidate_profiles = candidate_profiles,
  imputation_n = imputation_n,
  random_starts = random_starts,
  maximum_iterations = maximum_iterations,
  convergence_tolerance = convergence_tolerance,
  variance_floor = variance_floor,
  seed = stage8a_seed,
  lactate_transform = "log1p",
  creatinine_transform = creatinine_transform
)

if (file.exists(checkpoint_path)) {
  checkpoint <- readRDS(checkpoint_path)
  if (!identical(checkpoint$signature, checkpoint_signature)) {
    stop("Existing Stage 8A checkpoint does not match the frozen configuration.")
  }
  fit_records <- checkpoint$fit_records
  fit_metrics <- checkpoint$fit_metrics
  cat(sprintf("resumed checkpoint fits=%d\n", length(fit_records)))
}

for (imputation_index in seq_len(imputation_n)) {
  preprocessing <- preprocess_dataset(
    completed_datasets[[imputation_index]],
    creatinine_log1p
  )
  preprocessing_records[[imputation_index]] <- list(
    centers = preprocessing$centers,
    scales = preprocessing$scales
  )

  for (profiles in candidate_profiles) {
    record_key <- sprintf("imp%02d_k%d", imputation_index, profiles)
    if (record_key %in% names(fit_records)) {
      cat(sprintf("checkpoint imputation=%d profiles=%d\n", imputation_index, profiles))
      next
    }
    fit_seed <- stage8a_seed + imputation_index * 100L + profiles
    model <- fit_multistart(preprocessing$data, profiles, fit_seed)
    parameter_n <- 2L * profiles * length(primary_variables) + profiles - 1L
    aic <- -2 * model$log_likelihood + 2 * parameter_n
    bic <- -2 * model$log_likelihood + parameter_n * log(nrow(preprocessing$data))
    entropy <- normalized_entropy(model$posterior)
    profile_sizes <- tabulate(model$classification, nbins = profiles)
    fit_records[[record_key]] <- list(
      imputation = imputation_index,
      profiles = profiles,
      seed = fit_seed,
      model = model
    )
    fit_metrics[[record_key]] <- data.frame(
      imputation = imputation_index,
      profiles = profiles,
      n = nrow(preprocessing$data),
      parameters = parameter_n,
      converged = model$converged,
      iterations = model$iterations,
      log_likelihood = model$log_likelihood,
      aic = aic,
      bic = bic,
      entropy = entropy,
      smallest_profile_proportion = min(profile_sizes) / sum(profile_sizes),
      min_centroid_distance = minimum_centroid_distance(model$means),
      boundary_variance_n = model$boundary_variance_n,
      best_loglik_replications = model$best_log_likelihood_replications,
      initialization_warning_starts = model$initialization_warning_starts,
      seed = fit_seed,
      lactate_transform = "log1p",
      creatinine_transform = creatinine_transform
    )
    saveRDS(
      list(
        signature = checkpoint_signature,
        fit_records = fit_records,
        fit_metrics = fit_metrics
      ),
      checkpoint_path,
      compress = "gzip"
    )
    cat(sprintf("completed imputation=%d profiles=%d\n", imputation_index, profiles))
  }
}

fit_metrics_data <- do.call(rbind, fit_metrics)
row.names(fit_metrics_data) <- NULL
fit_metrics_data$bic_rank_within_imputation <- ave(
  fit_metrics_data$bic,
  fit_metrics_data$imputation,
  FUN = function(x) rank(x, ties.method = "min")
)

alignment_records <- list()
profile_size_rows <- list()
profile_mean_rows <- list()
profile_variance_rows <- list()
medoid_rows <- list()

for (profiles in candidate_profiles) {
  permutations <- all_permutations(seq_len(profiles))
  means_by_imputation <- lapply(seq_len(imputation_n), function(imputation_index) {
    fit_records[[sprintf("imp%02d_k%d", imputation_index, profiles)]]$model$means
  })

  pair_cost <- matrix(0, nrow = imputation_n, ncol = imputation_n)
  for (first_index in seq_len(imputation_n - 1L)) {
    for (second_index in (first_index + 1L):imputation_n) {
      alignment <- minimum_alignment(
        means_by_imputation[[first_index]],
        means_by_imputation[[second_index]],
        permutations
      )
      cost <- sum(
        (means_by_imputation[[second_index]][alignment, , drop = FALSE] -
          means_by_imputation[[first_index]])^2
      )
      pair_cost[first_index, second_index] <- cost
      pair_cost[second_index, first_index] <- cost
    }
  }
  medoid_imputation <- which.min(rowSums(pair_cost))
  medoid_means_original <- means_by_imputation[[medoid_imputation]]
  deterministic_order <- do.call(
    order,
    as.data.frame(
      medoid_means_original[, c(
        "lactate",
        "creatinine",
        "glucose",
        "wbc",
        "platelet_count",
        "bicarbonate"
      ), drop = FALSE]
    )
  )
  reference_means <- medoid_means_original[deterministic_order, , drop = FALSE]
  medoid_rows[[as.character(profiles)]] <- data.frame(
    profiles = profiles,
    medoid_imputation = medoid_imputation,
    summed_pair_cost = sum(pair_cost[medoid_imputation, ])
  )

  for (imputation_index in seq_len(imputation_n)) {
    record_key <- sprintf("imp%02d_k%d", imputation_index, profiles)
    model <- fit_records[[record_key]]$model
    alignment <- if (imputation_index == medoid_imputation) {
      deterministic_order
    } else {
      minimum_alignment(reference_means, model$means, permutations)
    }
    alignment_records[[record_key]] <- alignment

    for (aligned_profile in seq_len(profiles)) {
      original_profile <- alignment[aligned_profile]
      profile_n <- sum(model$classification == original_profile)
      profile_proportion <- profile_n / length(model$classification)
      profile_size_rows[[length(profile_size_rows) + 1L]] <- data.frame(
        imputation = imputation_index,
        profiles = profiles,
        profile = aligned_profile,
        n = profile_n,
        profile_proportion = profile_proportion,
        profile_percent = 100 * profile_proportion,
        meets_minimum_5_percent = profile_proportion >= minimum_profile_proportion,
        medoid_imputation = medoid_imputation
      )
      for (variable in primary_variables) {
        profile_mean_rows[[length(profile_mean_rows) + 1L]] <- data.frame(
          imputation = imputation_index,
          profiles = profiles,
          profile = aligned_profile,
          variable = variable,
          mean_z = model$means[original_profile, variable]
        )
        profile_variance_rows[[length(profile_variance_rows) + 1L]] <- data.frame(
          imputation = imputation_index,
          profiles = profiles,
          profile = aligned_profile,
          variable = variable,
          variance_z = model$variances[original_profile, variable]
        )
      }
    }
  }
}

profile_size_data <- do.call(rbind, profile_size_rows)
profile_mean_data <- do.call(rbind, profile_mean_rows)
profile_variance_data <- do.call(rbind, profile_variance_rows)
medoid_data <- do.call(rbind, medoid_rows)
row.names(medoid_data) <- NULL

profile_consensus <- aggregate(
  mean_z ~ profiles + profile + variable,
  data = profile_mean_data,
  FUN = mean
)
names(profile_consensus)[names(profile_consensus) == "mean_z"] <- "consensus_mean_z"

centroid_deviation_rows <- list()
for (profiles in candidate_profiles) {
  for (imputation_index in seq_len(imputation_n)) {
    for (profile in seq_len(profiles)) {
      observed <- profile_mean_data[
        profile_mean_data$profiles == profiles &
          profile_mean_data$imputation == imputation_index &
          profile_mean_data$profile == profile,
        c("variable", "mean_z")
      ]
      consensus <- profile_consensus[
        profile_consensus$profiles == profiles &
          profile_consensus$profile == profile,
        c("variable", "consensus_mean_z")
      ]
      matched <- merge(observed, consensus, by = "variable", sort = FALSE)
      centroid_deviation_rows[[length(centroid_deviation_rows) + 1L]] <- data.frame(
        imputation = imputation_index,
        profiles = profiles,
        profile = profile,
        centroid_deviation = sqrt(sum((matched$mean_z - matched$consensus_mean_z)^2))
      )
    }
  }
}
centroid_deviation_data <- do.call(rbind, centroid_deviation_rows)

cross_imputation_summary_rows <- lapply(candidate_profiles, function(profiles) {
  fit_subset <- fit_metrics_data[fit_metrics_data$profiles == profiles, ]
  size_subset <- profile_size_data[profile_size_data$profiles == profiles, ]
  deviation_subset <- centroid_deviation_data[
    centroid_deviation_data$profiles == profiles,
  ]
  consensus_subset <- profile_consensus[profile_consensus$profiles == profiles, ]

  proportion_ranges <- vapply(seq_len(profiles), function(profile) {
    profile_percent <- size_subset$profile_percent[size_subset$profile == profile]
    max(profile_percent) - min(profile_percent)
  }, numeric(1))
  consensus_ranges <- vapply(primary_variables, function(variable) {
    variable_means <- consensus_subset$consensus_mean_z[
      consensus_subset$variable == variable
    ]
    max(variable_means) - min(variable_means)
  }, numeric(1))
  inflammatory_contrast_n <- sum(
    consensus_ranges[inflammatory_variables] >= minimum_domain_contrast
  )
  metabolic_contrast_n <- sum(
    consensus_ranges[metabolic_variables] >= minimum_domain_contrast
  )

  data.frame(
    profiles = profiles,
    imputations_n = nrow(fit_subset),
    converged_n = sum(fit_subset$converged),
    boundary_fit_n = sum(fit_subset$boundary_variance_n > 0L),
    mean_aic = mean(fit_subset$aic),
    sd_aic = stats::sd(fit_subset$aic),
    mean_bic = mean(fit_subset$bic),
    sd_bic = stats::sd(fit_subset$bic),
    mean_entropy = mean(fit_subset$entropy),
    min_entropy = min(fit_subset$entropy),
    min_profile_proportion = min(size_subset$profile_proportion),
    max_centroid_deviation = max(deviation_subset$centroid_deviation),
    max_profile_proportion_range_pp = max(proportion_ranges),
    min_centroid_distance = min(fit_subset$min_centroid_distance),
    inflammatory_contrast_n = inflammatory_contrast_n,
    metabolic_contrast_n = metabolic_contrast_n,
    size_eligible = min(size_subset$profile_proportion) >= minimum_profile_proportion,
    structure_stable =
      max(deviation_subset$centroid_deviation) <= maximum_centroid_deviation &&
      max(proportion_ranges) <= maximum_profile_range_pp,
    cross_domain_interpretable =
      inflammatory_contrast_n >= 1L &&
      metabolic_contrast_n >= 1L &&
      min(fit_subset$min_centroid_distance) >= minimum_centroid_separation
  )
})
cross_imputation_summary <- do.call(rbind, cross_imputation_summary_rows)
cross_imputation_summary$admissible_for_later_selection <- with(
  cross_imputation_summary,
  converged_n == imputation_n &
    boundary_fit_n == 0L &
    size_eligible &
    structure_stable &
    cross_domain_interpretable
)

numeric_fit_columns <- c(
  "log_likelihood",
  "aic",
  "bic",
  "entropy",
  "smallest_profile_proportion",
  "min_centroid_distance"
)
fit_metrics_data[numeric_fit_columns] <- lapply(
  fit_metrics_data[numeric_fit_columns],
  function(x) round(x, 6)
)
profile_size_data[c("profile_proportion", "profile_percent")] <- lapply(
  profile_size_data[c("profile_proportion", "profile_percent")],
  function(x) round(x, 6)
)

stage8a_object <- list(
  stage = "Stage 8A-1",
  created_date = as.character(Sys.Date()),
  input_path = input_path,
  input_md5 = unname(tools::md5sum(input_path)),
  configuration = list(
    variables = primary_variables,
    lactate_transform = "log1p",
    creatinine_transform = creatinine_transform,
    creatinine_transform_rule = paste(
      "log1p only when median raw skewness across imputations exceeds",
      creatinine_skewness_threshold,
      "and log1p reduces median absolute skewness"
    ),
    standardization = "within-imputation z-score",
    distribution = "Gaussian latent profile",
    covariance_model = paste(
      "profile-specific means; profile-specific diagonal variances;",
      "zero within-profile covariances"
    ),
    candidate_profiles = candidate_profiles,
    random_starts = random_starts,
    maximum_iterations = maximum_iterations,
    convergence_tolerance = convergence_tolerance,
    variance_floor = variance_floor,
    seed = stage8a_seed,
    selection_performed = FALSE
  ),
  distribution_audit = distribution_audit,
  fit_records = fit_records,
  fit_metrics_by_imputation = fit_metrics_data,
  preprocessing = preprocessing_records,
  medoid_imputations = medoid_data,
  alignments = alignment_records,
  profile_size_by_imputation = profile_size_data,
  centroid_deviation_by_imputation = centroid_deviation_data,
  profile_variance_by_imputation = profile_variance_data,
  profile_consensus = profile_consensus,
  cross_imputation_summary = cross_imputation_summary,
  outcome_variables_used = FALSE,
  survival_information_used = FALSE,
  patient_level_profile_labels_exported = FALSE
)

model_tmp <- tempfile("stage8a_lpa_", tmpdir = dirname(model_path), fileext = ".rds")
fit_tmp <- tempfile("stage8a_fit_", tmpdir = dirname(fit_path), fileext = ".csv")
size_tmp <- tempfile("stage8a_size_", tmpdir = dirname(size_path), fileext = ".csv")

saveRDS(stage8a_object, model_tmp, compress = "gzip")
utils::write.csv(fit_metrics_data, fit_tmp, row.names = FALSE, na = "")
utils::write.csv(profile_size_data, size_tmp, row.names = FALSE, na = "")

verification <- readRDS(model_tmp)
verification_fit <- utils::read.csv(fit_tmp, stringsAsFactors = FALSE)
verification_size <- utils::read.csv(size_tmp, stringsAsFactors = FALSE)
if (
  length(verification$fit_records) != imputation_n * length(candidate_profiles) ||
  nrow(verification_fit) != imputation_n * length(candidate_profiles) ||
  nrow(verification_size) != imputation_n * sum(candidate_profiles) ||
  isTRUE(verification$outcome_variables_used) ||
  isTRUE(verification$survival_information_used) ||
  isTRUE(verification$patient_level_profile_labels_exported)
) {
  stop("Temporary Stage 8A outputs failed structural verification.")
}
size_check <- aggregate(
  cbind(n, profile_proportion) ~ imputation + profiles,
  data = verification_size,
  FUN = sum
)
if (
  any(size_check$n != nrow(imputation_object$data)) ||
  any(abs(size_check$profile_proportion - 1) > 1e-5)
) {
  stop("Profile size/proportion totals failed verification.")
}
forbidden_output_pattern <- "subject_id|stay_id|hadm_id|mortality|death|survival|outcome|los"
if (
  any(grepl(forbidden_output_pattern, names(verification_fit), ignore.case = TRUE)) ||
  any(grepl(forbidden_output_pattern, names(verification_size), ignore.case = TRUE))
) {
  stop("A forbidden patient-level or outcome field entered an aggregate output.")
}

if (!file.rename(model_tmp, model_path)) {
  stop("Could not finalize the protected Stage 8A model object.")
}
if (!file.rename(fit_tmp, fit_path)) {
  stop("Could not finalize the Stage 8A model-fit table.")
}
if (!file.rename(size_tmp, size_path)) {
  stop("Could not finalize the Stage 8A profile-size table.")
}
if (file.exists(checkpoint_path) && !file.remove(checkpoint_path)) {
  warning("Final outputs were written, but the temporary checkpoint could not be removed.")
}

cat("STAGE8A_LPA_FITTING_COMPLETED\n")
cat(sprintf(
  "fits=%d fit_rows=%d profile_size_rows=%d creatinine_transform=%s\n",
  length(fit_records),
  nrow(fit_metrics_data),
  nrow(profile_size_data),
  creatinine_transform
))
cat(sprintf(
  "all_converged=%s boundary_fit_n=%d no_patient_level_export=TRUE\n",
  all(fit_metrics_data$converged),
  sum(fit_metrics_data$boundary_variance_n > 0L)
))
