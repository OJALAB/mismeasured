# =========================================================================
# Tests for frequency weights in simex()
# =========================================================================

test_that("weights=1 gives same result as no weights (ME, Gaussian)", {
  set.seed(42)
  n <- 500
  x_true <- rnorm(n)
  y <- 1 + 2 * x_true + rnorm(n, sd = 0.5)
  df <- data.frame(y = y, x = x_true + rnorm(n, sd = 0.5))

  fit_no_wt <- simex(y ~ me(x, 0.5), data = df, B = 50, seed = 42)
  fit_wt1   <- simex(y ~ me(x, 0.5), data = df, B = 50, seed = 42,
                     weights = rep(1, n))

  expect_equal(coef(fit_wt1), coef(fit_no_wt), tolerance = 1e-10)
})

test_that("weights=1 gives same result as no weights (MC, standard)", {
  set.seed(42)
  n <- 500
  z <- rbinom(n, 1, 0.4)
  y <- 1 + 2 * z + rnorm(n)
  z_star <- z
  z_star[z == 0] <- rbinom(sum(z == 0), 1, 0.1)
  z_star[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)
  Pi <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
  df <- data.frame(y = y, z = factor(z_star))

  fit_no_wt <- simex(y ~ mc(z, Pi), data = df, method = "standard",
                     B = 50, seed = 42)
  fit_wt1   <- simex(y ~ mc(z, Pi), data = df, method = "standard",
                     B = 50, seed = 42, weights = rep(1, n))

  expect_equal(coef(fit_wt1), coef(fit_no_wt), tolerance = 1e-10)
})

test_that("expanded data matches weighted aggregated data (ME, Gaussian)", {
  set.seed(42)
  n <- 300
  x_true <- rnorm(n)
  y <- 1 + 2 * x_true + rnorm(n, sd = 0.5)
  x_obs <- x_true + rnorm(n, sd = 0.5)

  # Create duplicated (expanded) data — each obs repeated 1-3 times
  reps <- sample(1:3, n, replace = TRUE)
  df_expanded <- data.frame(
    y = rep(y, times = reps),
    x = rep(x_obs, times = reps)
  )

  # Weighted version uses original data with frequency weights
  df_weighted <- data.frame(y = y, x = x_obs)

  fit_expanded <- simex(y ~ me(x, 0.5), data = df_expanded,
                        B = 100, seed = 42, jackknife = FALSE)
  fit_weighted <- simex(y ~ me(x, 0.5), data = df_weighted,
                        B = 100, seed = 42, jackknife = FALSE,
                        weights = reps)

  expect_equal(coef(fit_weighted), coef(fit_expanded), tolerance = 0.15)
})

test_that("expanded data matches weighted aggregated data (MC, standard)", {
  set.seed(42)
  n <- 500
  z <- rbinom(n, 1, 0.4)
  y <- 1 + 2 * z + rnorm(n)
  z_star <- z
  z_star[z == 0] <- rbinom(sum(z == 0), 1, 0.1)
  z_star[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)
  Pi <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)

  reps <- sample(1:3, n, replace = TRUE)
  df_expanded <- data.frame(
    y = rep(y, times = reps),
    z = factor(rep(z_star, times = reps))
  )
  df_weighted <- data.frame(y = y, z = factor(z_star))

  fit_expanded <- simex(y ~ mc(z, Pi), data = df_expanded,
                        method = "standard", B = 100, seed = 42,
                        jackknife = FALSE)
  fit_weighted <- simex(y ~ mc(z, Pi), data = df_weighted,
                        method = "standard", B = 100, seed = 42,
                        jackknife = FALSE, weights = reps)

  expect_equal(coef(fit_weighted), coef(fit_expanded), tolerance = 0.15)
})

test_that("Poisson with weights works (ME)", {
  set.seed(42)
  n <- 1000
  x_true <- rnorm(n)
  y <- rpois(n, exp(0.5 + 0.8 * x_true))
  df <- data.frame(y = y, x = x_true + rnorm(n, sd = 0.3))

  fit <- simex(y ~ me(x, 0.3), family = poisson(), data = df,
               B = 50, seed = 42, weights = rep(2, n))

  expect_s3_class(fit, "simex")
  expect_equal(family(fit)$family, "poisson")
  expect_length(coef(fit), 2)
})

test_that("binomial with weights works (MC, standard)", {
  set.seed(42)
  n <- 1000
  z <- rbinom(n, 1, 0.4)
  x <- rnorm(n)
  prob <- plogis(0.5 + z + 0.5 * x)
  y <- rbinom(n, 1, prob)
  z_star <- z
  z_star[z == 0] <- rbinom(sum(z == 0), 1, 0.1)
  z_star[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)
  Pi <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
  df <- data.frame(y = y, z = factor(z_star), x = x)

  fit <- simex(y ~ mc(z, Pi) + x, family = binomial(), data = df,
               method = "standard", B = 50, seed = 42,
               weights = sample(1:3, n, replace = TRUE))

  expect_s3_class(fit, "simex")
  expect_equal(family(fit)$family, "binomial")
  expect_length(coef(fit), 3)
})
