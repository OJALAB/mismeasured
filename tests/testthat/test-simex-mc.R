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
