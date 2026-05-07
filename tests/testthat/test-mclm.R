test_that("mclm() matches mcglm(family='gaussian') (formula interface)", {
  set.seed(33)
  n  <- 1000
  Pi <- matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2)
  z  <- rbinom(n, 1, 0.4)
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, 0.1)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.1)
  x1 <- rnorm(n)
  y  <- -0.5 + 0.8 * z + 0.7 * x1 + rnorm(n, 0, 0.5)
  df <- data.frame(y = y, z = factor(z_hat), x1 = x1)

  fit_lm  <- mclm(y ~ mc(z, Pi) + x1, data = df)
  fit_glm <- mcglm(y ~ mc(z, Pi) + x1, data = df, family = "gaussian")

  expect_s3_class(fit_lm, "mcglm")
  expect_equal(fit_lm$family, fit_glm$family)
  expect_equal(fit_lm$coefficients, fit_glm$coefficients)
  expect_equal(fit_lm$vcov, fit_glm$vcov)
  expect_match(deparse(fit_lm$call)[1], "^mclm\\(")
})

test_that("mclm() matrix interface works", {
  set.seed(7)
  n  <- 800
  Pi <- matrix(c(0.85, 0.15, 0.10, 0.90), 2, 2)
  z  <- rbinom(n, 1, 0.5)
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, 0.15)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.10)
  x1 <- rnorm(n)
  y  <- 1 + 0.5 * z - 0.3 * x1 + rnorm(n, 0, 0.4)

  fit <- mclm(y, z_hat = as.integer(z_hat), x = cbind(1, x1),
              Pi = Pi, method = c("naive", "cs"))
  expect_s3_class(fit, "mcglm")
  expect_named(fit$coefficients, c("naive", "cs"))
  for (nm in names(fit$coefficients))
    expect_true(all(is.finite(fit$coefficients[[nm]])), info = nm)
})
