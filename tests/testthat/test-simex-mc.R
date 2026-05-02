test_that("simex with mc() improved (default) corrects attenuation (linear)", {
  set.seed(42)
  n <- 5000
  x <- rnorm(n)
  z <- rbinom(n, 1, 0.4)
  y <- 1 + 2 * z + 0.5 * x + rnorm(n)

  z_star <- z
  z_star[z == 0] <- rbinom(sum(z == 0), 1, 0.1)
  z_star[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)

  Pi <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
  df <- data.frame(y = y, z = factor(z_star), x = x)

  fit <- simex(y ~ mc(z, Pi) + x, data = df, B = 100, seed = 42)

  expect_s3_class(fit, "simex")
  expect_equal(fit$error.type, "mc")
  expect_equal(fit$method, "improved")

  # Improved MC-SIMEX should correct toward true z coefficient (2.0)
  naive_bias <- abs(fit$naive.coefficients[1] - 2.0)
  simex_bias <- abs(coef(fit)[1] - 2.0)
  expect_lt(simex_bias, naive_bias)
})

test_that("simex with mc() standard method works", {
  set.seed(42)
  n <- 1000
  z <- rbinom(n, 1, 0.4)
  y <- 1 + 2 * z + rnorm(n)
  z_star <- z
  z_star[z == 0] <- rbinom(sum(z == 0), 1, 0.1)
  z_star[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)

  Pi <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
  df <- data.frame(y = y, z = factor(z_star))

  fit <- simex(y ~ mc(z, Pi), data = df, method = "standard", B = 50, seed = 42)
  expect_equal(fit$method, "standard")
  expect_equal(fit$extrapolation.method, "quadratic")
})

test_that("lambda='optimal' works for improved mc()", {
  set.seed(42)
  n <- 1000
  z <- rbinom(n, 1, 0.4)
  y <- 1 + 2 * z + rnorm(n)
  z_star <- z
  z_star[z == 0] <- rbinom(sum(z == 0), 1, 0.1)
  z_star[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)
  Pi <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
  df <- data.frame(y = y, z = factor(z_star))

  fit <- simex(y ~ mc(z, Pi), data = df, lambda = "optimal", B = 50, seed = 42)
  expect_equal(fit$method, "improved")
  expect_true(!is.null(fit$c.lambda))
})

test_that("S3 methods work for mc()-based simex", {
  set.seed(42)
  n <- 500
  z <- rbinom(n, 1, 0.4)
  y <- 1 + 2 * z + rnorm(n)
  z_star <- z
  z_star[z == 0] <- rbinom(sum(z == 0), 1, 0.1)
  z_star[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)
  Pi <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
  df <- data.frame(y = y, z = factor(z_star))

  fit <- simex(y ~ mc(z, Pi), data = df, B = 50, seed = 42)

  expect_length(coef(fit), 2)
  expect_true(is.matrix(vcov(fit)))
  expect_equal(nobs(fit), n)

  s <- summary(fit)
  expect_s3_class(s, "summary.simex")

  ci <- confint(fit)
  expect_equal(nrow(ci), 2)

  newdf <- data.frame(z = factor(c(0, 1)))
  pred <- predict(fit, newdata = newdf)
  expect_length(pred, 2)
})

test_that("refit errors for improved mc()", {
  set.seed(42)
  n <- 500
  z <- rbinom(n, 1, 0.4)
  y <- 1 + 2 * z + rnorm(n)
  z_star <- z
  z_star[z == 0] <- rbinom(sum(z == 0), 1, 0.1)
  z_star[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)
  Pi <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
  df <- data.frame(y = y, z = factor(z_star))

  fit <- simex(y ~ mc(z, Pi), data = df, B = 50, seed = 42)
  expect_error(refit(fit), "not applicable")
})


# =========================================================================
# Multiple mc() terms
# =========================================================================

test_that("simex with two mc() terms (standard) corrects attenuation", {
  set.seed(42)
  n <- 3000
  x <- rnorm(n)
  z1 <- rbinom(n, 1, 0.4)
  z2 <- rbinom(n, 1, 0.6)
  y <- 1 + 2 * z1 + 1.5 * z2 + 0.5 * x + rnorm(n)

  # Misclassify z1
  z1_star <- z1
  z1_star[z1 == 0] <- rbinom(sum(z1 == 0), 1, 0.1)
  z1_star[z1 == 1] <- 1 - rbinom(sum(z1 == 1), 1, 0.15)
  Pi1 <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)

  # Misclassify z2
  z2_star <- z2
  z2_star[z2 == 0] <- rbinom(sum(z2 == 0), 1, 0.08)
  z2_star[z2 == 1] <- 1 - rbinom(sum(z2 == 1), 1, 0.12)
  Pi2 <- matrix(c(0.92, 0.08, 0.12, 0.88), 2, 2)

  df <- data.frame(y = y, z1 = factor(z1_star), z2 = factor(z2_star), x = x)

  fit <- simex(y ~ mc(z1, Pi1) + mc(z2, Pi2) + x, data = df,
               method = "standard", B = 100, seed = 42)

  expect_s3_class(fit, "simex")
  expect_equal(fit$error.type, "mc")
  expect_equal(fit$method, "standard")
  expect_length(fit$mc.terms, 2)

  # Should have coefficients for both z factors + x + intercept
  expect_length(coef(fit), 4)

  # SIMEX correction should reduce bias for z1 coefficient
  naive_z1_bias <- abs(fit$naive.coefficients["z11"] - 2.0)
  simex_z1_bias <- abs(coef(fit)["z11"] - 2.0)
  expect_lt(simex_z1_bias, naive_z1_bias)
})

test_that("multi-mc defaults to standard method", {
  set.seed(42)
  n <- 500
  z1 <- rbinom(n, 1, 0.4)
  z2 <- rbinom(n, 1, 0.6)
  y <- 1 + z1 + z2 + rnorm(n)

  Pi1 <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
  Pi2 <- matrix(c(0.92, 0.08, 0.12, 0.88), 2, 2)
  df <- data.frame(y = y, z1 = factor(z1), z2 = factor(z2))

  # Default method should be "standard" for multi-mc
  fit <- simex(y ~ mc(z1, Pi1) + mc(z2, Pi2), data = df, B = 20, seed = 42)
  expect_equal(fit$method, "standard")
})

test_that("multi-mc errors with method='improved'", {
  set.seed(42)
  n <- 500
  z1 <- rbinom(n, 1, 0.4)
  z2 <- rbinom(n, 1, 0.6)
  y <- 1 + z1 + z2 + rnorm(n)

  Pi1 <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
  Pi2 <- matrix(c(0.92, 0.08, 0.12, 0.88), 2, 2)
  df <- data.frame(y = y, z1 = factor(z1), z2 = factor(z2))

  expect_error(
    simex(y ~ mc(z1, Pi1) + mc(z2, Pi2), data = df,
          method = "improved", B = 20, seed = 42),
    "single mc"
  )
})

test_that("S3 methods work for multi-mc simex", {
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

  expect_length(coef(fit), 3)  # z11 + z21 + intercept
  expect_true(is.matrix(vcov(fit)))
  expect_equal(nobs(fit), n)
  expect_length(residuals(fit), n)
  expect_length(fitted(fit), n)

  s <- summary(fit)
  expect_s3_class(s, "summary.simex")

  ci <- confint(fit)
  expect_equal(nrow(ci), 3)

  pred <- predict(fit, newdata = data.frame(z1 = factor(c(0, 1)),
                                            z2 = factor(c(0, 1))))
  expect_length(pred, 2)
})

test_that("multi-mc with three-level factor works", {
  set.seed(42)
  n <- 2000
  z1 <- factor(sample(c("A", "B", "C"), n, replace = TRUE, prob = c(0.4, 0.3, 0.3)))
  z2 <- rbinom(n, 1, 0.5)
  y <- 1 + ifelse(z1 == "B", 1, ifelse(z1 == "C", 2, 0)) + z2 + rnorm(n)

  Pi1 <- matrix(c(0.8, 0.1, 0.1, 0.1, 0.8, 0.1, 0.1, 0.1, 0.8), 3, 3)
  Pi2 <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)

  # Apply misclassification
  z1_star <- z1
  for (lev in levels(z1)) {
    idx <- z1 == lev
    col_idx <- which(levels(z1) == lev)
    z1_star[idx] <- sample(levels(z1), sum(idx), replace = TRUE,
                           prob = Pi1[, col_idx])
  }

  df <- data.frame(y = y, z1 = z1_star, z2 = factor(z2))

  fit <- simex(y ~ mc(z1, Pi1) + mc(z2, Pi2), data = df,
               method = "standard", B = 50, seed = 42)

  expect_s3_class(fit, "simex")
  # z1B, z1C (2 dummies) + z21 (1 dummy) + intercept = 4

  expect_length(coef(fit), 4)
})
