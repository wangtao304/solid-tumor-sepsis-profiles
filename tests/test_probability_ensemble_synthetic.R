#!/usr/bin/env Rscript

# Synthetic unit test: verifies probability averaging and demonstrates that
# averaging probabilities is not equivalent to averaging coefficients.
set.seed(20260729L)
n <- 20L
p <- 6L
components <- 10L
x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
beta <- matrix(stats::rnorm((p + 1L) * components, sd = 0.2), p + 1L, components)

component_probability <- vapply(seq_len(components), function(j) {
  stats::plogis(beta[1L, j] + x %*% beta[-1L, j])
}, numeric(n))
ensemble_probability <- rowMeans(component_probability)

mean_beta_probability <- stats::plogis(rowMeans(beta)[1L] + x %*% rowMeans(beta)[-1L])
stopifnot(
  length(ensemble_probability) == n,
  all(is.finite(ensemble_probability)),
  all(ensemble_probability > 0 & ensemble_probability < 1),
  !isTRUE(all.equal(ensemble_probability, as.numeric(mean_beta_probability)))
)

assignment <- ifelse(ensemble_probability >= 0.5, 2L, 1L)
stopifnot(all(assignment %in% 1:2))
cat("SYNTHETIC_PROBABILITY_ENSEMBLE_TEST_OK\n")
