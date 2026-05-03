# =========================================================================
# Regression tests for MC design-matrix construction (issue #6)
# =========================================================================

# --- Helper: generate standard MC test data ---
.make_mc_data <- function(n = 1000, seed = 42) {
  set.seed(seed)
  x <- rnorm(n)
  z <- rbinom(n, 1, 0.4)
  y <- 1 + 2 * z + 0.5 * x + rnorm(n)
  z_star <- z
  z_star[z == 0] <- rbinom(sum(z == 0), 1, 0.1)
  z_star[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)
  Pi <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
  list(
    df = data.frame(y = y, z = factor(z_star), x = x),
    Pi = Pi
  )
}


# =========================================================================
# P1: Formula semantics preserved (I(), interactions, no-intercept)
# =========================================================================

test_that("mc() with I(x^2) preserves squared term", {
  d <- .make_mc_data(n = 500)
  df <- d$df
  Pi <- d$Pi

  # Should have coefficients: z1, (Intercept), x, I(x^2) = 4 coefficients

  fit <- simex(y ~ mc(z, Pi) + x + I(x^2), data = df, B = 20, seed = 42)
  expect_length(coef(fit), 4)
  expect_true("I(x^2)" %in% names(coef(fit)))
})

test_that("mc() with interaction among non-mc variables works", {
  set.seed(42)
  n <- 500
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  z <- rbinom(n, 1, 0.4)
  y <- 1 + 2 * z + x1 + x2 + 0.5 * x1 * x2 + rnorm(n)
  z_star <- z
  z_star[z == 0] <- rbinom(sum(z == 0), 1, 0.1)
  z_star[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)
  Pi <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
  df <- data.frame(y = y, z = factor(z_star), x1 = x1, x2 = x2)

  fit <- simex(y ~ mc(z, Pi) + x1 * x2, data = df, B = 20, seed = 42)
  # z1, (Intercept), x1, x2, x1:x2 = 5 coefficients
  expect_length(coef(fit), 5)
  expect_true("x1:x2" %in% names(coef(fit)))
})


# =========================================================================
# P1: predict(fit, newdata) == fitted(fit) for training data
# =========================================================================

test_that("predict on training data equals fitted values (single mc)", {
  d <- .make_mc_data(n = 500)
  fit <- simex(y ~ mc(z, d$Pi) + x, data = d$df, B = 20, seed = 42)

  pred <- predict(fit, newdata = d$df)
  expect_equal(pred, fitted(fit), tolerance = 1e-10)
})

test_that("predict on training data equals fitted values (multi mc)", {
  set.seed(42)
  n <- 500
  z1 <- rbinom(n, 1, 0.4)
  z2 <- rbinom(n, 1, 0.6)
  y <- 1 + z1 + z2 + rnorm(n)
  Pi1 <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
  Pi2 <- matrix(c(0.92, 0.08, 0.12, 0.88), 2, 2)
  df <- data.frame(y = y, z1 = factor(z1), z2 = factor(z2))

  fit <- simex(y ~ mc(z1, Pi1) + mc(z2, Pi2), data = df,
               method = "standard", B = 20, seed = 42)

  pred <- predict(fit, newdata = df)
  expect_equal(pred, fitted(fit), tolerance = 1e-10)
})


# =========================================================================
# P1: predict with subset factor levels
# =========================================================================

test_that("predict works when newdata has only one factor level", {
  d <- .make_mc_data(n = 500)
  fit <- simex(y ~ mc(z, d$Pi) + x, data = d$df, B = 20, seed = 42)

  # Predict for only level "1" — should not treat it as baseline
  newdf <- data.frame(z = factor(c(1, 1)), x = c(0.5, -0.5))
  pred <- predict(fit, newdata = newdf)
  expect_length(pred, 2)

  # Compare with manual computation: the z=1 dummy should be 1
  cf <- coef(fit)
  expected <- cf["(Intercept)"] + cf["1"] * 1 + cf["x"] * c(0.5, -0.5)
  expect_equal(unname(pred), unname(expected), tolerance = 1e-10)
})


# =========================================================================
# P1: refit() consistency
# =========================================================================

test_that("refit produces consistent fitted values (single mc)", {
  d <- .make_mc_data(n = 500)
  fit <- simex(y ~ mc(z, d$Pi) + x, data = d$df,
               method = "standard", B = 50, seed = 42)
  refit_fit <- refit(fit, extrapolation = "linear")

  # Fitted values should match coefficients %*% xi_hat
  expect_equal(length(fitted(refit_fit)), 500)
  expect_false(all(fitted(refit_fit) == fitted(fit)))  # different extrapolation
})

test_that("refit produces consistent fitted values (multi mc)", {
  set.seed(42)
  n <- 500
  z1 <- rbinom(n, 1, 0.4)
  z2 <- rbinom(n, 1, 0.6)
  y <- 1 + z1 + z2 + rnorm(n)
  Pi1 <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
  Pi2 <- matrix(c(0.92, 0.08, 0.12, 0.88), 2, 2)
  df <- data.frame(y = y, z1 = factor(z1), z2 = factor(z2))

  fit <- simex(y ~ mc(z1, Pi1) + mc(z2, Pi2), data = df,
               method = "standard", B = 20, seed = 42)
  refit_fit <- refit(fit, extrapolation = "linear")

  # Verify fitted values are computed from correct design matrix
  xi_hat <- fit$xi.hat
  eta <- as.numeric(xi_hat %*% coef(refit_fit))
  expect_equal(fitted(refit_fit), fit$family$linkinv(eta), tolerance = 1e-10)
})


# =========================================================================
# P2: mc matrix validation
# =========================================================================

test_that("mc matrix with negative entries is rejected", {
  set.seed(42)
  n <- 100
  df <- data.frame(y = rnorm(n), z = factor(rbinom(n, 1, 0.5)))
  Pi_bad <- matrix(c(1.5, -0.5, 0.2, 0.8), 2, 2)

  expect_error(
    simex(y ~ mc(z, Pi_bad), data = df),
    "\\[0, 1\\]"
  )
})

test_that("mc matrix with entries > 1 is rejected", {
  set.seed(42)
  n <- 100
  df <- data.frame(y = rnorm(n), z = factor(rbinom(n, 1, 0.5)))
  Pi_bad <- matrix(c(1.1, -0.1, 0.2, 0.8), 2, 2)

  expect_error(
    simex(y ~ mc(z, Pi_bad), data = df),
    "\\[0, 1\\]"
  )
})


# =========================================================================
# mc() with I(x^2) and predict consistency
# =========================================================================

test_that("predict with I(x^2) matches fitted values", {
  d <- .make_mc_data(n = 500)
  df <- d$df
  Pi <- d$Pi

  fit <- simex(y ~ mc(z, Pi) + x + I(x^2), data = df, B = 20, seed = 42)
  pred <- predict(fit, newdata = df)
  expect_equal(pred, fitted(fit), tolerance = 1e-10)
})
