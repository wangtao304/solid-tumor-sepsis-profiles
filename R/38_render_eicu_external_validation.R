#!/usr/bin/env Rscript

# Re-render the Stage 10C-2 aggregate forest plot without database access or
# model fitting. Inputs contain model-level results only.

options(stringsAsFactors = FALSE, warn = 1)

primary_path <- "outputs/tables/38_eicu_profile_outcome_results.csv"
complete_case_path <- "outputs/tables/39_eicu_complete_case_sensitivity.csv"
figure_path <- "outputs/figures/figure_eicu_external_validation.pdf"

primary <- utils::read.csv(primary_path, check.names = FALSE)
complete_case <- utils::read.csv(complete_case_path, check.names = FALSE)
plot_mi <- primary[primary$status == "ESTIMATED", ]
plot_cc <- complete_case[complete_case$status == "ESTIMATED", ]
if (
  nrow(plot_mi) != 3L || nrow(plot_cc) != 3L ||
  any(!is.finite(plot_mi$estimate)) || any(!is.finite(plot_cc$estimate))
) {
  stop("Aggregate Stage 10C-2 plot inputs are invalid.")
}

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
    sprintf("MICE primary (n=%s)", format(plot_mi$n[1L], big.mark = ",")),
    sprintf("Complete case (n=%s)", format(plot_cc$n[1L], big.mark = ","))
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

if (!file.exists(figure_tmp) || file.info(figure_tmp)$size <= 1000) {
  unlink(figure_tmp)
  stop("Re-rendered Stage 10C-2 figure failed file validation.")
}
if (!file.copy(figure_tmp, figure_path, overwrite = TRUE)) {
  unlink(figure_tmp)
  stop("Could not replace the Stage 10C-2 validation figure.")
}
unlink(figure_tmp)
cat("STAGE10C2_FIGURE_RENDERED_FROM_AGGREGATE_RESULTS\n")
