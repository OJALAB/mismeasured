test_that("simex with me() corrects attenuation (Gaussian)", {
  set.seed(42)
  n <- 2000
  x_true <- rnorm(n)
  y <- 1 + 2 * x_true + rnorm(n, sd = 0.5)
  x_obs <- x_true + rnorm(n, sd = 0.5)
  df <- data.frame(y = y, x = x_obs)

  fit <- simex(y ~ me(x, 0.5), data = df, B = 100, seed = 42)

  expect_s3_class(fit, "simex")
  expect_equal(fit$error.type, "me")

  # SIMEX should correct toward true value (2.0)
  naive_bias <- abs(fit$naive.coefficients["x"] - 2.0)
  simex_bias <- abs(coef(fit)["x"] - 2.0)
  expect_lt(simex_bias, naive_bias)
})

test_that("simex with multiple me() terms works", {
  set.seed(42)
  n <- 2000
  x1_true <- rnorm(n)
  x2_true <- rnorm(n)
  y <- 1 + x1_true + 2 * x2_true + rnorm(n, sd = 0.5)
  df <- data.frame(y = y, x1 = x1_true + rnorm(n, sd = 0.3),
                   x2 = x2_true + rnorm(n, sd = 0.5))

  fit <- simex(y ~ me(x1, 0.3) + me(x2, 0.5), data = df, B = 100, seed = 42)

  expect_length(coef(fit), 3)  # intercept + x1 + x2
  expect_equal(fit$error.type, "me")
  expect_length(fit$me.terms, 2)
})

test_that("S3 methods work for me()-based simex", {
  set.seed(42)
  n <- 500
  x_true <- rnorm(n)
  y <- 1 + 2 * x_true + rnorm(n, sd = 0.5)
  df <- data.frame(y = y, x = x_true + rnorm(n, sd = 0.5))

  fit <- simex(y ~ me(x, 0.5), data = df, B = 50, seed = 42)

  expect_length(coef(fit), 2)
  expect_true(is.matrix(vcov(fit)))
  expect_length(residuals(fit), n)
  expect_length(fitted(fit), n)
  expect_equal(nobs(fit), n)
  expect_equal(family(fit)$family, "gaussian")

  s <- summary(fit)
  expect_s3_class(s, "summary.simex")
  expect_equal(ncol(s$coefficients), 4)  # Estimate, SE, t, p

  ci <- confint(fit)
  expect_equal(nrow(ci), 2)

  pred <- predict(fit, newdata = data.frame(x = c(0, 1)))
  expect_length(pred, 2)
})

test_that("refit changes extrapolation for me() simex", {
  set.seed(42)
  n <- 500
  x_true <- rnorm(n)
  y <- 1 + 2 * x_true + rnorm(n, sd = 0.5)
  df <- data.frame(y = y, x = x_true + rnorm(n, sd = 0.5))

  fit_quad <- simex(y ~ me(x, 0.5), data = df, B = 50, seed = 42)
  fit_lin <- refit(fit_quad, extrapolation = "linear")

  expect_false(identical(coef(fit_quad), coef(fit_lin)))
})
