#!/usr/bin/env Rscript

# Re-render R1C figures from aggregate CSVs only. No database connection,
# patient-level data, model fitting, or statistical recalculation occurs.

options(stringsAsFactors = FALSE)
out_dir <- "outputs/revision"
curve_path <- file.path(out_dir, "R1C_time_varying_hr_curve.csv")
points_path <- file.path(out_dir, "R1C_time_varying_hr_points.csv")
if (!all(file.exists(c(curve_path, points_path)))) stop("Required aggregate CSV is missing.")
all_curve <- utils::read.csv(curve_path, check.names = FALSE)
point_rows <- utils::read.csv(points_path, check.names = FALSE)
if (nrow(all_curve) != 324L || nrow(point_rows) != 15L) stop("Unexpected aggregate input dimensions.")

draw_main <- function() {
  d <- all_curve[all_curve$scenario == "midpoint", ]
  fixed <- point_rows[point_rows$scenario == "midpoint", ]
  yr <- range(c(d$CI_lower, d$CI_upper, 1), finite = TRUE)
  plot(NA, xlim = c(1, 28), ylim = yr, log = "y", xaxt = "n",
       xlab = "ICU day", ylab = "Hazard ratio for Profile 2 vs Profile 1",
       main = "Time-varying profile association after the 24-hour landmark",
       las = 1, bty = "l")
  axis(1, at = c(2, 7, 14, 21, 28), labels = c("2", "7", "14", "21", "28"))
  polygon(c(d$ICU_day, rev(d$ICU_day)), c(d$CI_lower, rev(d$CI_upper)),
          col = grDevices::adjustcolor("#2A6FBB", alpha.f = 0.18), border = NA)
  lines(d$ICU_day, d$HR, col = "#1F5A94", lwd = 2.2)
  points(fixed$ICU_day, fixed$HR, pch = 21, bg = "white", col = "#1F5A94", cex = 1.0)
  abline(h = 1, col = "#555555", lty = 2)
  abline(v = 1, col = "#333333", lty = 3)
  text(1.15, exp(log(yr[1]) + 0.92 * diff(log(yr))), "24-h\nlandmark",
       pos = 4, cex = 0.72, col = "#333333")
  legend("topright", legend = c("Conditional HR(t)", "95% pointwise CI", "No association"),
         col = c("#1F5A94", grDevices::adjustcolor("#2A6FBB", alpha.f = 0.35), "#555555"),
         lwd = c(2.2, 8, 1), lty = c(1, 1, 2), bty = "n", cex = 0.86)
  mtext("Midpoint date convention; SOFA-adjusted; patient-clustered robust covariance", side = 3, line = 0.25, cex = 0.78)
}

draw_date <- function() {
  cols <- c(midpoint = "#1F5A94", earliest = "#D97706", latest = "#5B4B8A")
  ltys <- c(midpoint = 1, earliest = 2, latest = 3)
  yr <- range(c(all_curve$HR, 1), finite = TRUE)
  plot(NA, xlim = c(1, 28), ylim = yr, log = "y", xaxt = "n",
       xlab = "ICU day", ylab = "Hazard ratio for Profile 2 vs Profile 1",
       main = "Date-precision sensitivity of the time-varying association",
       las = 1, bty = "l")
  axis(1, at = c(2, 7, 14, 21, 28), labels = c("2", "7", "14", "21", "28"))
  for (s in c("midpoint", "earliest", "latest")) {
    d <- all_curve[all_curve$scenario == s, ]
    lines(d$ICU_day, d$HR, col = cols[s], lty = ltys[s], lwd = 2)
  }
  abline(h = 1, col = "#555555", lty = 2)
  abline(v = 1, col = "#333333", lty = 3)
  text(1.15, exp(log(yr[1]) + 0.92 * diff(log(yr))), "24-h\nlandmark",
       pos = 4, cex = 0.72, col = "#333333")
  legend("topright", legend = c("Midpoint", "Earliest legitimate time", "Latest legitimate time"),
         col = cols, lty = ltys, lwd = 2, bty = "n", cex = 0.84)
  mtext("Scenario curves are date-precision sensitivity analyses, not confidence bounds", side = 3, line = 0.25, cex = 0.78)
}

pdf_main <- file.path(out_dir, "R1C_time_varying_hr.pdf")
tiff_main <- file.path(out_dir, "R1C_time_varying_hr.tiff")
pdf_date <- file.path(out_dir, "R1C_date_precision_sensitivity.pdf")
tiff_date <- file.path(out_dir, "R1C_date_precision_sensitivity.tiff")
grDevices::pdf(pdf_main, width = 7.2, height = 5.4, useDingbats = FALSE,
               family = "sans", paper = "special")
draw_main(); grDevices::dev.off()
grDevices::tiff(tiff_main, width = 7.2, height = 5.4, units = "in",
                res = 300, compression = "lzw", bg = "white")
draw_main(); grDevices::dev.off()
grDevices::pdf(pdf_date, width = 7.2, height = 5.4, useDingbats = FALSE,
               family = "sans", paper = "special")
draw_date(); grDevices::dev.off()
grDevices::tiff(tiff_date, width = 7.2, height = 5.4, units = "in",
                res = 300, compression = "lzw", bg = "white")
draw_date(); grDevices::dev.off()
if (!all(file.info(c(pdf_main, tiff_main, pdf_date, tiff_date))$size > 0)) stop("Figure render failed.")
cat("R1C_FIGURES_REBUILT_FROM_AGGREGATES\n")
