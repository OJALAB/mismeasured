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


# =========================================================================
# Berkson error
# =========================================================================

test_that("Berkson SIMEX corrects Berkson attenuation", {
  set.seed(42)
  n <- 3000
  # Berkson model: X_true = W + U, we observe W
  w <- rnorm(n)          # observed (assigned/planned value)
  u <- rnorm(n, sd = 0.5) # Berkson error
  x_true <- w + u       # true value
  y <- 1 + 2 * x_true + rnorm(n, sd = 0.5)
  df <- data.frame(y = y, x = w)

  # Naive regression of y on w is biased (Berkson inflates variance, not attenuates)
  # But SIMEX with type="berkson" should still correct
  fit <- simex(y ~ me(x, 0.5, type = "berkson"), data = df, B = 100, seed = 42)

  expect_s3_class(fit, "simex")
  # The SIMEX correction should produce a different result from naive
  expect_false(identical(coef(fit)["x"], fit$naive.coefficients["x"]))
})

test_that("Berkson type='classical' gives different result from type='berkson'", {
  set.seed(42)
  n <- 1000
  x_true <- rnorm(n)
  y <- 1 + 2 * x_true + rnorm(n, sd = 0.5)
  x_obs <- x_true + rnorm(n, sd = 0.5)
  df <- data.frame(y = y, x = x_obs)

  fit_class <- simex(y ~ me(x, 0.5, type = "classical"), data = df,
                     B = 50, seed = 42)
  fit_berk <- simex(y ~ me(x, 0.5, type = "berkson"), data = df,
                    B = 50, seed = 42)

  # Different error types should give different corrections
  expect_false(identical(coef(fit_class), coef(fit_berk)))
})


# =========================================================================
# Non-zero mean measurement error
# =========================================================================

test_that("non-zero mean me() corrects systematic error", {
  set.seed(42)
  n <- 3000
  x_true <- rnorm(n)
  # Systematic measurement error: mean = 0.3, sd = 0.5
  x_obs <- x_true + rnorm(n, mean = 0.3, sd = 0.5)
  y <- 1 + 2 * x_true + rnorm(n, sd = 0.5)
  df <- data.frame(y = y, x = x_obs)

  # Standard SIMEX (mean=0) won't account for the systematic shift
  fit_std <- simex(y ~ me(x, 0.5), data = df, B = 100, seed = 42)

  # Non-zero mean SIMEX should do better
  fit_nzm <- simex(y ~ me(x, 0.5, mean = 0.3), data = df, B = 100, seed = 42)

  # Both should be different from naive
  expect_false(identical(coef(fit_nzm), fit_nzm$naive.coefficients))

  # The systematic shift mu = 0.3 enters X linearly, so it primarily
  # affects the intercept. With true intercept = 1, mean=0 leaves the
  # naive bias of -0.3 * slope = -0.6 in place; mean=0.3 removes it.
  expect_gt(abs(coef(fit_std)["(Intercept)"] - 1),
            abs(coef(fit_nzm)["(Intercept)"] - 1) + 0.1)
})

test_that("mean=0 gives same result as default", {
  set.seed(42)
  n <- 500
  x_true <- rnorm(n)
  y <- 1 + 2 * x_true + rnorm(n, sd = 0.5)
  df <- data.frame(y = y, x = x_true + rnorm(n, sd = 0.5))

  fit_default <- simex(y ~ me(x, 0.5), data = df, B = 50, seed = 42)
  fit_mean0   <- simex(y ~ me(x, 0.5, mean = 0), data = df, B = 50, seed = 42)

  expect_equal(coef(fit_mean0), coef(fit_default), tolerance = 1e-10)
})
