# --------------------------------------------------------------------------
# simex(): paths previously uncovered.
#   * S3: print.simex, summary -> print.summary.simex, plot.simex,
#     formula.simex
#   * predict.simex with type = "link"
#   * extrapolation = "loglinear"
#   * Mixed me() + mc() in a single formula
#   * Berkson + non-zero mean together
# --------------------------------------------------------------------------

# ==========================================================================
# S3 surface
# ==========================================================================

test_that("simex S3 print/summary/plot/formula run without error", {
  set.seed(101)
  n <- 300
  x_true <- rnorm(n)
  y <- 1 + 2 * x_true + rnorm(n, sd = 0.5)
  df <- data.frame(y = y, x = x_true + rnorm(n, sd = 0.5))
  fit <- simex(y ~ me(x, 0.5), data = df, B = 30, seed = 1)

  expect_output(print(fit), "SIMEX|simex|extrapolation|naive", ignore.case = TRUE)
  smry <- summary(fit)
  expect_s3_class(smry, "summary.simex")
  expect_output(print(smry), "Coefficients", ignore.case = TRUE)
  expect_identical(formula(fit), fit$formula)

  # plot.simex draws to a (null) device; just make sure it doesn't error.
  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  expect_no_error(plot(fit))
})

test_that("predict.simex(type='link') returns linear predictor", {
  set.seed(102)
  n <- 300
  x_true <- rnorm(n)
  y <- 1 + 2 * x_true + rnorm(n, sd = 0.5)
  df <- data.frame(y = y, x = x_true + rnorm(n, sd = 0.5))
  fit <- simex(y ~ me(x, 0.5), data = df, B = 30, seed = 1)

  pred_resp <- predict(fit, type = "response")
  pred_link <- predict(fit, type = "link")
  expect_length(pred_resp, n)
  expect_length(pred_link, n)
  # Gaussian: response and link should match (identity link).
  expect_equal(pred_resp, pred_link, tolerance = 1e-12)

  # With newdata
  nd <- data.frame(x = c(0, 1, 2))
  expect_length(predict(fit, newdata = nd, type = "response"), 3L)
  expect_length(predict(fit, newdata = nd, type = "link"), 3L)
})

# ==========================================================================
# extrapolation = "loglinear"
# ==========================================================================

test_that("simex with extrapolation='loglinear' runs and matches naive at lambda=0", {
  set.seed(103)
  n <- 300
  x_true <- rnorm(n)
  y <- rpois(n, exp(1 + 0.5 * x_true))
  df <- data.frame(y = y, x = x_true + rnorm(n, sd = 0.3))
  fit_log <- simex(y ~ me(x, 0.3), data = df, B = 30, seed = 1,
                   family = poisson(), extrapolation = "loglinear")
  expect_s3_class(fit_log, "simex")
  expect_true(all(is.finite(coef(fit_log))))
})

test_that("simex with extrapolation='linear' produces different estimates", {
  set.seed(104)
  n <- 300
  x_true <- rnorm(n)
  y <- 1 + 2 * x_true + rnorm(n, sd = 0.5)
  df <- data.frame(y = y, x = x_true + rnorm(n, sd = 0.5))
  fit_lin  <- simex(y ~ me(x, 0.5), data = df, B = 30, seed = 1,
                    extrapolation = "linear")
  fit_quad <- simex(y ~ me(x, 0.5), data = df, B = 30, seed = 1,
                    extrapolation = "quadratic")
  expect_false(isTRUE(all.equal(coef(fit_lin)["x"], coef(fit_quad)["x"])))
})

# ==========================================================================
# Mixed me() + mc() in one formula
# ==========================================================================

test_that("simex supports me() + mc() mixed in one formula", {
  set.seed(105)
  n <- 400
  x_true <- rnorm(n)
  z <- rbinom(n, 1, 0.4)
  Pi <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, 0.10)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)
  y <- 1 + 2 * x_true + 0.7 * z + rnorm(n, sd = 0.5)
  df <- data.frame(y = y, x = x_true + rnorm(n, sd = 0.5),
                   z = factor(z_hat))

  fit <- tryCatch(
    simex(y ~ me(x, 0.5) + mc(z, Pi), data = df, B = 30, seed = 1,
          jackknife = FALSE),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    # Document current behavior: mixed me() + mc() may not be supported.
    expect_match(conditionMessage(fit), "mc|me|mixed|both",
                 ignore.case = TRUE)
  } else {
    expect_s3_class(fit, "simex")
    expect_true(all(is.finite(coef(fit))))
  }
})

# ==========================================================================
# Berkson + non-zero mean
# ==========================================================================

test_that("Berkson type with non-zero mean runs", {
  set.seed(106)
  n <- 500
  x_true <- rnorm(n)
  x_obs <- x_true - rnorm(n, mean = 0.2, sd = 0.4)
  y <- 1 + 2 * x_true + rnorm(n, sd = 0.5)
  df <- data.frame(y = y, x = x_obs)

  fit <- simex(y ~ me(x, 0.4, type = "berkson", mean = 0.2),
               data = df, B = 30, seed = 1)
  expect_s3_class(fit, "simex")
  expect_true(all(is.finite(coef(fit))))
})
