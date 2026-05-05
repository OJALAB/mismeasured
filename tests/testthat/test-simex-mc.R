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

test_that(".mat_power_r handles asymmetric Pi with complex eigenvalues", {
  # A column-stochastic Pi that is not symmetric and whose eigenvalues
  # are complex; .mat_power_r used to fail with
  # "unimplemented complex function" because it called sign() on
  # complex eigenvalues. Fixed by using the principal branch d^power
  # in complex arithmetic and stripping numerical imaginary parts.
  # Concrete 4x4 column-stochastic Pi whose eigendecomposition contains
  # the complex pair -0.14 +/- 0.118i. The previous .mat_power_r used
  # sign(d) which is undefined on complex values.
  Pi <- matrix(c(0.125, 0.176, 0.270, 0.429,
                 0.075, 0.332, 0.349, 0.244,
                 0.586, 0.058, 0.192, 0.164,
                 0.294, 0.164, 0.329, 0.213),
               nrow = 4, byrow = FALSE)
  Pi <- sweep(Pi, 2, colSums(Pi), "/")  # ensure column-stochastic
  stopifnot(any(abs(Im(eigen(Pi)$values)) > 1e-6))

  # Integer power: should equal Pi %*% Pi exactly.
  P2 <- mismeasured:::.mat_power_r(Pi, 2)
  expect_equal(P2, Pi %*% Pi, tolerance = 1e-12)

  # Fractional power: should be real, finite, and roughly inverse to
  # Pi^(-power) up to numerical precision.
  P_half <- mismeasured:::.mat_power_r(Pi, 0.5)
  expect_true(all(is.finite(P_half)))
  expect_equal(P_half %*% P_half, Pi, tolerance = 1e-6)
})

test_that("improved MC-SIMEX handles K=4 estimated Pi without crashing", {
  # Regression for the bug observed in simulations/simulation_study.R:
  # estimating Pi from a small validation sample gave asymmetric Pi
  # matrices whose eigendecomposition produced complex eigenvalues, and
  # .mat_power_r crashed.
  set.seed(20260707)
  K <- 4L; n <- 2000
  pi_z <- c(.5, .25, .15, .10)
  Pi   <- 0.8 * diag(K) + (0.2 / 3) * (matrix(1, K, K) - diag(K))
  q    <- rnorm(n)
  z    <- sample(0:3, n, TRUE, pi_z)
  zh   <- vapply(z, function(zi) sample(0:3, 1, prob = Pi[, zi + 1L]), 0L)
  y    <- rpois(n, exp(c(0, 1, -0.9, 0.2)[z + 1L] + 0.8 - 0.7 * q))

  n_fail <- 0L
  for (b in 1:10) {
    vi <- sample.int(n, 200)
    Pi_hat <- matrix(0, K, K)
    for (l in seq_len(K)) {
      sel <- z[vi] == (l - 1L)
      if (any(sel)) {
        Pi_hat[, l] <- tabulate(zh[vi][sel] + 1L, nbins = K) / sum(sel)
      } else {
        Pi_hat[, l] <- 1 / K
      }
    }
    Pi_hat <- pmax(Pi_hat, 1e-3)
    Pi_hat <- sweep(Pi_hat, 2, colSums(Pi_hat), "/")
    df <- data.frame(y = y, z = factor(zh, levels = as.character(0:3)), q1 = q)
    res <- tryCatch(
      simex(y ~ mc(z, Pi_hat) + q1, family = poisson(), data = df,
            method = "improved", jackknife = FALSE),
      error = function(e) { n_fail <<- n_fail + 1L; NULL }
    )
  }
  expect_equal(n_fail, 0L)
})

test_that("simex with K-level mc() improved corrects dummy attenuation", {
  set.seed(42)
  n <- 6000
  levs <- c("A", "B", "C")
  x_true <- sample(0:2, n, replace = TRUE, prob = c(0.45, 0.35, 0.20))
  y <- 1 + 1.5 * (x_true == 1) - 1.0 * (x_true == 2) + rnorm(n, sd = 0.8)

  Pi <- matrix(c(
    0.82, 0.10, 0.08,
    0.10, 0.80, 0.10,
    0.08, 0.10, 0.82
  ), 3, 3)

  z_star <- integer(n)
  for (k in 0:2) {
    idx <- which(x_true == k)
    z_star[idx] <- sample(0:2, length(idx), replace = TRUE,
                          prob = Pi[, k + 1L])
  }
  df <- data.frame(y = y, z = factor(levs[z_star + 1L], levels = levs))

  fit <- simex(y ~ mc(z, Pi), data = df, lambda = 0, jackknife = FALSE)

  expect_s3_class(fit, "simex")
  expect_equal(fit$method, "improved")
  expect_length(fit$pi.vec, 3)
  expect_equal(dim(fit$correction.matrix$lambda_0), c(2L, 2L))

  truth <- c(B = 1.5, C = -1.0)
  naive_bias <- abs(fit$naive.coefficients[names(truth)] - truth)
  simex_bias <- abs(coef(fit)[names(truth)] - truth)
  expect_lt(max(simex_bias), max(naive_bias))
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

test_that("lambda='optimal' errors for K-level improved mc()", {
  set.seed(42)
  n <- 200
  z <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  y <- rnorm(n)
  Pi <- matrix(c(
    0.8, 0.1, 0.1,
    0.1, 0.8, 0.1,
    0.1, 0.1, 0.8
  ), 3, 3)
  df <- data.frame(y = y, z = z)

  expect_error(
    simex(y ~ mc(z, Pi), data = df, lambda = "optimal"),
    "binary"
  )
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


# =========================================================================
# Response misclassification: mc(y, Pi) ~ ...
# =========================================================================

test_that("response mc corrects attenuation (logistic)", {
  set.seed(42)
  n <- 3000
  x <- rnorm(n)
  y_true <- rbinom(n, 1, plogis(-0.5 + 1.2 * x))
  y_obs <- y_true
  y_obs[y_true == 0] <- rbinom(sum(y_true == 0), 1, 0.15)
  y_obs[y_true == 1] <- 1 - rbinom(sum(y_true == 1), 1, 0.1)
  Pi_y <- matrix(c(0.85, 0.15, 0.1, 0.9), 2, 2)
  df <- data.frame(y = factor(y_obs), x = x)

  fit <- simex(mc(y, Pi_y) ~ x, family = binomial(), data = df,
               method = "standard", B = 100, seed = 42)

  expect_s3_class(fit, "simex")
  expect_equal(fit$error.type, "mc")
  expect_equal(fit$method, "standard")
  expect_false(is.null(fit$response.mc))

  naive_bias <- abs(fit$naive.coefficients["x"] - 1.2)
  simex_bias <- abs(coef(fit)["x"] - 1.2)
  expect_lt(simex_bias, naive_bias)
})

test_that("response mc + covariate mc works together", {
  set.seed(42)
  n <- 3000
  x <- rnorm(n)
  z <- rbinom(n, 1, 0.4)
  y_true <- rbinom(n, 1, plogis(-0.5 + 1.2 * z + 0.5 * x))
  y_obs <- y_true
  y_obs[y_true == 0] <- rbinom(sum(y_true == 0), 1, 0.15)
  y_obs[y_true == 1] <- 1 - rbinom(sum(y_true == 1), 1, 0.1)
  Pi_y <- matrix(c(0.85, 0.15, 0.1, 0.9), 2, 2)
  z_star <- z
  z_star[z == 0] <- rbinom(sum(z == 0), 1, 0.1)
  z_star[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)
  Pi_z <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
  df <- data.frame(y = factor(y_obs), z = factor(z_star), x = x)

  fit <- simex(mc(y, Pi_y) ~ mc(z, Pi_z) + x, family = binomial(), data = df,
               method = "standard", B = 100, seed = 42)

  expect_s3_class(fit, "simex")
  expect_length(coef(fit), 3)  # intercept + z1 + x
  expect_length(fit$mc.terms, 1)  # only RHS mc terms
  expect_false(is.null(fit$response.mc))
})

test_that("response mc defaults to standard method", {
  set.seed(42)
  n <- 500
  y <- factor(rbinom(n, 1, 0.5))
  x <- rnorm(n)
  Pi_y <- matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2)
  df <- data.frame(y = y, x = x)

  fit <- simex(mc(y, Pi_y) ~ x, family = binomial(), data = df,
               B = 20, seed = 42)
  expect_equal(fit$method, "standard")
})

test_that("response mc errors with method='improved'", {
  set.seed(42)
  n <- 500
  y <- factor(rbinom(n, 1, 0.5))
  x <- rnorm(n)
  Pi_y <- matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2)
  df <- data.frame(y = y, x = x)

  expect_error(
    simex(mc(y, Pi_y) ~ x, family = binomial(), data = df,
          method = "improved", B = 20, seed = 42),
    "single mc.*covariate"
  )
})

test_that("S3 methods work for response mc", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  y <- factor(rbinom(n, 1, plogis(0.5 * x)))
  Pi_y <- matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2)
  df <- data.frame(y = y, x = x)

  fit <- simex(mc(y, Pi_y) ~ x, family = binomial(), data = df,
               method = "standard", B = 20, seed = 42)

  expect_length(coef(fit), 2)
  expect_true(is.matrix(vcov(fit)))
  expect_equal(nobs(fit), n)

  s <- summary(fit)
  expect_s3_class(s, "summary.simex")

  ci <- confint(fit)
  expect_equal(nrow(ci), 2)

  pred <- predict(fit, newdata = data.frame(x = c(0, 1)), type = "response")
  expect_length(pred, 2)
})
