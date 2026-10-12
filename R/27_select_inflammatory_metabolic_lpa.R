# Stage 8A-2: outcome-blind LPA model selection and aggregate profile description.

options(stringsAsFactors = FALSE)

input_path <- file.path(
  "data",
  "processed",
  "inflammatory_metabolic_lpa_solid_2to6_mi10_seed20260722.rds"
)
selection_path <- file.path("outputs", "tables", "27_lpa_profile_selection.csv")
characteristics_path <- file.path(
  "outputs",
  "tables",
  "28_lpa_profile_characteristics.csv"
)
figure_path <- file.path(
  "outputs",
  "figures",
  "figure_lpa_profile_patterns.pdf"
)

required_variables <- c(
  "wbc",
  "platelet_count",
  "lactate",
  "glucose",
  "creatinine",
  "bicarbonate"
)
display_variables <- c(
  wbc = "WBC",
  platelet_count = "platelet",
  lactate = "lactate",
  glucose = "glucose",
  creatinine = "creatinine",
  bicarbonate = "bicarbonate"
)
candidate_profiles <- 2:6
minimum_profile_proportion <- 0.05

if (!file.exists(input_path)) {
  stop("Approved Stage 8A-1 LPA object was not found: ", input_path)
}
dir.create(dirname(selection_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(figure_path), recursive = TRUE, showWarnings = FALSE)

lpa_object <- readRDS(input_path)
if (
  !identical(lpa_object$stage, "Stage 8A-1") ||
  !identical(lpa_object$configuration$variables, required_variables) ||
  !identical(lpa_object$configuration$candidate_profiles, candidate_profiles) ||
  !identical(lpa_object$configuration$selection_performed, FALSE)
) {
  stop("Stage 8A-1 object does not match the frozen Stage 8A-2 specification.")
}
if (
  !identical(lpa_object$outcome_variables_used, FALSE) ||
  !identical(lpa_object$survival_information_used, FALSE) ||
  !identical(lpa_object$patient_level_profile_labels_exported, FALSE)
) {
  stop("Stage 8A-1 safeguard flags do not permit Stage 8A-2 selection.")
}

fit_metrics <- lpa_object$fit_metrics_by_imputation
profile_sizes <- lpa_object$profile_size_by_imputation
cross_summary <- lpa_object$cross_imputation_summary

if (
  nrow(fit_metrics) != 50L ||
  nrow(profile_sizes) != 200L ||
  !all(fit_metrics$converged) ||
  any(fit_metrics$boundary_variance_n != 0L) ||
  anyNA(fit_metrics) ||
  anyNA(profile_sizes)
) {
  stop("Stage 8A-1 aggregate fit or profile-size QC failed.")
}

best_aic_counts <- table(factor(
  fit_metrics$profiles[
    ave(fit_metrics$aic, fit_metrics$imputation, FUN = min) == fit_metrics$aic
  ],
  levels = candidate_profiles
))
best_bic_counts <- table(factor(
  fit_metrics$profiles[
    ave(fit_metrics$bic, fit_metrics$imputation, FUN = min) == fit_metrics$bic
  ],
  levels = candidate_profiles
))

profile_range_rows <- lapply(
  split(profile_sizes, interaction(
    profile_sizes$profiles,
    profile_sizes$profile,
    drop = TRUE
  )),
  function(data) {
    data.frame(
      profiles = data$profiles[1L],
      profile = data$profile[1L],
      mean_profile_percent = mean(data$profile_percent),
      min_profile_percent = min(data$profile_percent),
      max_profile_percent = max(data$profile_percent),
      profile_proportion_range_pp =
        max(data$profile_percent) - min(data$profile_percent)
    )
  }
)
profile_ranges <- do.call(rbind, profile_range_rows)
row.names(profile_ranges) <- NULL
profile_ranges <- profile_ranges[order(
  profile_ranges$profiles,
  profile_ranges$profile
), ]

model_summary <- cross_summary
model_summary$aic_best_imputations_n <- as.integer(
  best_aic_counts[as.character(model_summary$profiles)]
)
model_summary$bic_best_imputations_n <- as.integer(
  best_bic_counts[as.character(model_summary$profiles)]
)

eligible_profiles <- model_summary$profiles[
  model_summary$admissible_for_later_selection
]
if (length(eligible_profiles) != 1L) {
  stop(
    "Stage 8A-2 requires exactly one candidate to pass all frozen screens; found ",
    length(eligible_profiles),
    "."
  )
}
recommended_profiles <- eligible_profiles[1L]

selection_reason <- vapply(model_summary$profiles, function(profiles) {
  row <- model_summary[model_summary$profiles == profiles, ]
  if (profiles == recommended_profiles) {
    return(paste(
      "Recommended: only candidate meeting the minimum-size,",
      "cross-imputation stability, and cross-domain interpretability screens."
    ))
  }
  failures <- c()
  if (!row$size_eligible) failures <- c(failures, "profile size below 5%")
  if (!row$structure_stable) failures <- c(failures, "cross-imputation instability")
  if (!row$cross_domain_interpretable) {
    failures <- c(failures, "insufficient cross-domain contrast")
  }
  paste0("Not recommended: ", paste(failures, collapse = "; "), ".")
}, character(1L))
model_summary$selection_reason <- selection_reason
model_summary$recommended <- model_summary$profiles == recommended_profiles

selection_output <- merge(
  profile_ranges,
  model_summary[c(
    "profiles",
    "imputations_n",
    "converged_n",
    "boundary_fit_n",
    "mean_aic",
    "sd_aic",
    "mean_bic",
    "sd_bic",
    "mean_entropy",
    "min_entropy",
    "min_profile_proportion",
    "max_centroid_deviation",
    "max_profile_proportion_range_pp",
    "min_centroid_distance",
    "size_eligible",
    "structure_stable",
    "cross_domain_interpretable",
    "admissible_for_later_selection",
    "aic_best_imputations_n",
    "bic_best_imputations_n",
    "recommended",
    "selection_reason"
  )],
  by = "profiles",
  sort = FALSE
)
selection_output <- selection_output[order(
  selection_output$profiles,
  selection_output$profile
), ]
selection_output$minimum_profile_percent_across_candidate <-
  100 * selection_output$min_profile_proportion
selection_output$meets_minimum_5_percent <-
  selection_output$min_profile_percent >= 100 * minimum_profile_proportion
selection_output <- selection_output[c(
  "profiles",
  "profile",
  "imputations_n",
  "converged_n",
  "boundary_fit_n",
  "mean_aic",
  "sd_aic",
  "mean_bic",
  "sd_bic",
  "mean_entropy",
  "min_entropy",
  "mean_profile_percent",
  "min_profile_percent",
  "max_profile_percent",
  "profile_proportion_range_pp",
  "minimum_profile_percent_across_candidate",
  "max_profile_proportion_range_pp",
  "max_centroid_deviation",
  "min_centroid_distance",
  "meets_minimum_5_percent",
  "size_eligible",
  "structure_stable",
  "cross_domain_interpretable",
  "admissible_for_later_selection",
  "aic_best_imputations_n",
  "bic_best_imputations_n",
  "recommended",
  "selection_reason"
)]

aligned_mean_rows <- list()
for (imputation in 1:10) {
  key <- sprintf("imp%02d_k%d", imputation, recommended_profiles)
  fit_record <- lpa_object$fit_records[[key]]
  alignment <- lpa_object$alignments[[key]]
  if (
    is.null(fit_record) ||
    is.null(alignment) ||
    length(alignment) != recommended_profiles
  ) {
    stop("Missing Stage 8A-1 fit or alignment for ", key, ".")
  }
  aligned_means <- fit_record$model$means[alignment, required_variables, drop = FALSE]
  aligned_mean_rows[[length(aligned_mean_rows) + 1L]] <- data.frame(
    imputation = imputation,
    profile = rep(seq_len(recommended_profiles), each = length(required_variables)),
    variable = rep(required_variables, times = recommended_profiles),
    mean_z = as.vector(t(aligned_means))
  )
}
aligned_means <- do.call(rbind, aligned_mean_rows)

characteristic_rows <- lapply(
  split(aligned_means, interaction(
    aligned_means$profile,
    aligned_means$variable,
    drop = TRUE
  )),
  function(data) {
    data.frame(
      profile = data$profile[1L],
      variable_internal = data$variable[1L],
      standardized_mean = mean(data$mean_z),
      standardized_sd_across_imputations = stats::sd(data$mean_z),
      standardized_min = min(data$mean_z),
      standardized_max = max(data$mean_z)
    )
  }
)
characteristics <- do.call(rbind, characteristic_rows)
row.names(characteristics) <- NULL

consensus_check <- merge(
  characteristics,
  lpa_object$profile_consensus[
    lpa_object$profile_consensus$profiles == recommended_profiles,
    c("profile", "variable", "consensus_mean_z")
  ],
  by.x = c("profile", "variable_internal"),
  by.y = c("profile", "variable"),
  all.x = TRUE,
  sort = FALSE
)
if (
  anyNA(consensus_check$consensus_mean_z) ||
  max(abs(
    consensus_check$standardized_mean - consensus_check$consensus_mean_z
  )) > 1e-10
) {
  stop("Reconstructed aligned profile means do not match the stored consensus.")
}
characteristics <- consensus_check

stress_direction <- c(
  wbc = 1,
  platelet_count = -1,
  lactate = 1,
  glucose = 1,
  creatinine = 1,
  bicarbonate = -1
)
stress_scores <- vapply(seq_len(recommended_profiles), function(profile) {
  data <- characteristics[characteristics$profile == profile, ]
  mean(data$standardized_mean * stress_direction[data$variable_internal])
}, numeric(1L))
lower_profile <- which.min(stress_scores)
higher_profile <- which.max(stress_scores)
profile_labels <- rep("Intermediate inflammatory-metabolic stress", recommended_profiles)
profile_labels[lower_profile] <- "Lower inflammatory-metabolic stress"
profile_labels[higher_profile] <- "Higher inflammatory-metabolic stress"

profile_interpretations <- rep(
  "Feature-based inflammatory-metabolic pattern.",
  recommended_profiles
)
profile_interpretations[lower_profile] <- paste(
  "Lower WBC, lactate, glucose, and creatinine; higher bicarbonate;",
  "platelet near the cohort mean."
)
profile_interpretations[higher_profile] <- paste(
  "Higher WBC, lactate, glucose, and creatinine; lower bicarbonate;",
  "modestly lower platelet."
)

selected_size_ranges <- profile_ranges[
  profile_ranges$profiles == recommended_profiles,
]
characteristics <- merge(
  characteristics,
  selected_size_ranges[c(
    "profile",
    "mean_profile_percent",
    "min_profile_percent",
    "max_profile_percent"
  )],
  by = "profile",
  all.x = TRUE,
  sort = FALSE
)
characteristics$recommended_profiles <- recommended_profiles
characteristics$profile_label <- profile_labels[characteristics$profile]
characteristics$variable <- unname(
  display_variables[characteristics$variable_internal]
)
characteristics$profile_stress_score <- stress_scores[characteristics$profile]
characteristics$feature_based_interpretation <-
  profile_interpretations[characteristics$profile]
characteristics <- characteristics[order(
  characteristics$profile,
  match(characteristics$variable_internal, required_variables)
), ]
characteristics_output <- characteristics[c(
  "recommended_profiles",
  "profile",
  "profile_label",
  "variable",
  "standardized_mean",
  "standardized_sd_across_imputations",
  "standardized_min",
  "standardized_max",
  "mean_profile_percent",
  "min_profile_percent",
  "max_profile_percent",
  "profile_stress_score",
  "feature_based_interpretation"
)]

numeric_selection_columns <- names(selection_output)[vapply(
  selection_output,
  is.numeric,
  logical(1L)
)]
selection_output[numeric_selection_columns] <- lapply(
  selection_output[numeric_selection_columns],
  function(x) round(x, 6)
)
numeric_characteristic_columns <- names(characteristics_output)[vapply(
  characteristics_output,
  is.numeric,
  logical(1L)
)]
characteristics_output[numeric_characteristic_columns] <- lapply(
  characteristics_output[numeric_characteristic_columns],
  function(x) round(x, 6)
)

selection_tmp <- tempfile(
  "stage8a2_selection_",
  tmpdir = dirname(selection_path),
  fileext = ".csv"
)
characteristics_tmp <- tempfile(
  "stage8a2_characteristics_",
  tmpdir = dirname(characteristics_path),
  fileext = ".csv"
)
figure_tmp <- tempfile(
  "stage8a2_patterns_",
  tmpdir = dirname(figure_path),
  fileext = ".pdf"
)

utils::write.csv(selection_output, selection_tmp, row.names = FALSE, na = "")
utils::write.csv(
  characteristics_output,
  characteristics_tmp,
  row.names = FALSE,
  na = ""
)

plot_data <- characteristics[order(
  characteristics$profile,
  match(characteristics$variable_internal, required_variables)
), ]
profile_colors <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
plot_y_min <- min(plot_data$standardized_min) - 0.18
plot_y_max <- max(plot_data$standardized_max) + 0.18

grDevices::pdf(
  figure_tmp,
  width = 9.2,
  height = 5.9,
  family = "sans",
  title = "Inflammatory-metabolic latent profile patterns",
  useDingbats = FALSE
)
graphics::par(
  mar = c(6.0, 4.8, 4.2, 1.2),
  mgp = c(2.7, 0.8, 0),
  las = 1,
  xaxs = "i",
  yaxs = "i"
)
graphics::plot(
  NA,
  xlim = c(0.65, length(required_variables) + 0.35),
  ylim = c(plot_y_min, plot_y_max),
  xaxt = "n",
  xlab = "",
  ylab = "Standardized mean (z score)",
  main = "Inflammatory-metabolic latent profile patterns",
  cex.main = 1.25,
  font.main = 2,
  col.axis = "#333333",
  col.lab = "#222222"
)
graphics::mtext(
  "Recommended 2-profile solution across 10 imputed datasets",
  side = 3,
  line = 1.0,
  cex = 0.9,
  col = "#555555"
)
graphics::abline(h = pretty(c(plot_y_min, plot_y_max)), col = "#E6E6E6", lwd = 0.8)
graphics::abline(h = 0, col = "#666666", lwd = 1.1)
graphics::axis(
  1,
  at = seq_along(required_variables),
  labels = unname(display_variables[required_variables]),
  las = 2,
  tick = FALSE,
  line = -0.5,
  col.axis = "#333333"
)
for (profile in seq_len(recommended_profiles)) {
  data <- plot_data[plot_data$profile == profile, ]
  data <- data[match(required_variables, data$variable_internal), ]
  x <- seq_along(required_variables) + (profile - 1.5) * 0.06
  graphics::segments(
    x,
    data$standardized_min,
    x,
    data$standardized_max,
    col = grDevices::adjustcolor(profile_colors[profile], alpha.f = 0.65),
    lwd = 2
  )
  graphics::lines(
    x,
    data$standardized_mean,
    type = "o",
    pch = if (profile == 1L) 16 else 17,
    lwd = 2.4,
    cex = 1.05,
    col = profile_colors[profile]
  )
}
legend_labels <- vapply(seq_len(recommended_profiles), function(profile) {
  size_row <- selected_size_ranges[selected_size_ranges$profile == profile, ]
  sprintf(
    "Profile %d: %s (mean %.1f%%; range %.1f-%.1f%%)",
    profile,
    profile_labels[profile],
    size_row$mean_profile_percent,
    size_row$min_profile_percent,
    size_row$max_profile_percent
  )
}, character(1L))
graphics::legend(
  "topleft",
  inset = c(0.01, 0.02),
  legend = legend_labels,
  col = profile_colors[seq_len(recommended_profiles)],
  lwd = 2.4,
  pch = c(16, 17)[seq_len(recommended_profiles)],
  bty = "n",
  cex = 0.82
)
graphics::mtext(
  paste(
    "Points: mean aligned z scores; vertical bars: min-max across imputations.",
    "Lactate and creatinine were log1p-transformed before standardization."
  ),
  side = 1,
  line = 4.8,
  cex = 0.72,
  col = "#555555",
  adj = 0
)
grDevices::dev.off()

verification_selection <- utils::read.csv(
  selection_tmp,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
verification_characteristics <- utils::read.csv(
  characteristics_tmp,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (
  nrow(verification_selection) != sum(candidate_profiles) ||
  nrow(verification_characteristics) !=
    recommended_profiles * length(required_variables) ||
  !all(verification_selection$recommended ==
    (verification_selection$profiles == recommended_profiles)) ||
  any(verification_characteristics$recommended_profiles != recommended_profiles) ||
  file.info(figure_tmp)$size < 1000
) {
  stop("Stage 8A-2 temporary output verification failed.")
}

if (!file.copy(selection_tmp, selection_path, overwrite = TRUE)) {
  stop("Could not finalize Stage 8A-2 profile-selection table.")
}
if (!file.copy(characteristics_tmp, characteristics_path, overwrite = TRUE)) {
  stop("Could not finalize Stage 8A-2 profile-characteristics table.")
}
if (!file.copy(figure_tmp, figure_path, overwrite = TRUE)) {
  stop("Could not finalize Stage 8A-2 profile-pattern figure.")
}
unlink(c(selection_tmp, characteristics_tmp, figure_tmp))

cat("STAGE8A2_PROFILE_SELECTION_COMPLETED\n")
cat(sprintf(
  paste0(
    "recommended_profiles=%d selection_rows=%d characteristic_rows=%d ",
    "outcome_blind=TRUE no_patient_level_export=TRUE\n"
  ),
  recommended_profiles,
  nrow(selection_output),
  nrow(characteristics_output)
))
