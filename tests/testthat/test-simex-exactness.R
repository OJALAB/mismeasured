# Gaussian-exactness anchor for the improved MC-SIMEX.
# At lambda = 0 with B = 1 the correction is a deterministic function of
# the naive fit, so for identity-link gaussian models the corrected
# coefficients must match the true-z oracle fit up to sampling noise.
# The logistic companion documents that the correction is only a
# first-order approximation for non-identity links.

test_that("improved MC-SIMEX is exact for gaussian with binary z", {
  set.seed(101)
  n <- 5000
  z <- rbinom(n, 1, 0.4)
  y <- 1 + 2 * z + rnorm(n)
  z_star <- ifelse(z == 1, 1 - rbinom(n, 1, 0.12), rbinom(n, 1, 0.1))
  Pi <- matrix(c(0.9, 0.1, 0.12, 0.88), 2, 2)
  df <- data.frame(y = y, z = factor(z_star))

  fit <- simex(y ~ mc(z, Pi), data = df, lambda = 0, B = 1,
               jackknife = FALSE)
  oracle <- glm(y ~ z)

  expect_lt(max(abs(sort(unname(coef(fit))) - sort(unname(coef(oracle))))),
            0.1)
  # and the naive fit is visibly attenuated, so the anchor has teeth
  naive <- glm(y ~ z_star)
  expect_gt(abs(coef(naive)[2] - coef(oracle)[2]), 0.3)
})

test_that("improved MC-SIMEX is exact for gaussian with K = 3 z", {
  set.seed(202)
  n <- 5000
  z <- sample(0:2, n, replace = TRUE, prob = c(0.5, 0.3, 0.2))
  y <- 1 + 2 * (z == 1) + 3 * (z == 2) + rnorm(n)
  Pi <- matrix(c(0.85, 0.10, 0.05,
                 0.08, 0.84, 0.08,
                 0.05, 0.10, 0.85), 3, 3)
  # misclassify z via the columns of Pi
  z_star <- vapply(z, function(zi)
    sample(0:2, 1, prob = Pi[, zi + 1]), integer(1))
  df <- data.frame(y = y, z = factor(z_star, levels = 0:2))

  fit <- simex(y ~ mc(z, Pi), data = df, method = "improved",
               lambda = 0, B = 1, jackknife = FALSE)
  oracle <- glm(y ~ factor(z, levels = 0:2))

  co <- unname(coef(fit))
  co_or <- unname(coef(oracle))
  # parameter order: dummies first, intercept last vs oracle intercept first
  expect_lt(max(abs(sort(co) - sort(co_or))), 0.15)
})

test_that("improved MC-SIMEX is close but not exact for logistic", {
  set.seed(303)
  n <- 20000
  z <- rbinom(n, 1, 0.4)
  x <- rnorm(n)
  y <- rbinom(n, 1, plogis(1 + 3 * z + 1.5 * x))
  z_star <- ifelse(z == 1, 1 - rbinom(n, 1, 0.10), rbinom(n, 1, 0.10))
  Pi <- matrix(c(0.9, 0.1, 0.10, 0.90), 2, 2)
  df <- data.frame(y = y, z = factor(z_star), x = x)

  fit <- simex(y ~ mc(z, Pi) + x, data = df, family = binomial(),
               method = "improved", lambda = 0, B = 1, jackknife = FALSE)
  oracle <- glm(y ~ z + x, family = binomial())
  naive <- glm(y ~ z_star + x, family = binomial())

  co <- coef(fit)
  z_coef <- co[[1]]  # dummy coefficients come first in the improved fit

  # improved beats naive on the misclassified coefficient...
  expect_lt(abs(z_coef - coef(oracle)[2]), abs(coef(naive)[2] - coef(oracle)[2]))
  # ...but the first-order approximation leaves visible bias somewhere
  expect_gt(max(abs(sort(unname(co)) - sort(unname(coef(oracle))))), 0.15)
})
