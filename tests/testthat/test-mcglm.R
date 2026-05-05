# --------------------------------------------------------------------------
# Tests for the mcglm() interface and S3 methods
# --------------------------------------------------------------------------

# --- Input validation ---

test_that("mcglm errors on mismatched dimensions", {
  expect_error(
    mcglm( 1:10, z_hat = 1:5, x = matrix(1, 10, 1), family = "poisson",
          method = "naive")
  )
})

test_that("mcglm errors when correction params missing", {
  expect_error(
    mcglm( rpois(10, 1), z_hat = rbinom(10, 1, 0.5),
          x = cbind(1, rnorm(10)), family = "poisson",
          method = c("naive", "bca")),
    "supply"
  )
})

# --- Binary interface ---

test_that("mcglm binary returns correct structure", {
  set.seed(1)
  n <- 500
  x <- cbind(1, rnorm(n))
  z <- rbinom(n, 1, 0.4)
  eta <- 0.8 * z - 0.5 * x[, 1] + 0.7 * x[, 2]
  y <- rpois(n, exp(eta))
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, 0.10)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)

  fit <- mcglm(y, z_hat = z_hat, x = x, family = "poisson",
               method = c("naive", "bca", "bcm", "cs"),
               p01 = 0.10, p10 = 0.15, pi_z = 0.4)

  expect_s3_class(fit, "mcglm")
  expect_named(fit$coefficients, c("naive", "bca", "bcm", "cs"))
  expect_equal(fit$K, 2L)
  expect_equal(fit$n, n)
  expect_equal(fit$p, 3)

  for (nm in names(fit$coefficients)) {
    expect_named(fit$coefficients[[nm]])
    expect_length(fit$coefficients[[nm]], 3)
  }
})

test_that("mcglm auto-detects binary case", {
  set.seed(2)
  n <- 200
  x <- cbind(1, rnorm(n))
  z_hat <- rbinom(n, 1, 0.5)
  y <- rpois(n, exp(0.5 * z_hat + x %*% c(-0.3, 0.4)))

  fit <- mcglm(y, z_hat = z_hat, x = x, family = "poisson", method = "naive")
  expect_equal(fit$K, 2L)
})

# --- Multicategory interface ---

test_that("mcglm multicategory returns correct structure", {
  set.seed(3)
  n <- 500
  K <- 3
  pi_z <- c(0.5, 0.3, 0.2)
  z <- sample(0:(K - 1), n, replace = TRUE, prob = pi_z)
  Pi <- matrix(c(0.8, 0.1, 0.1,
                 0.1, 0.8, 0.1,
                 0.1, 0.1, 0.8), nrow = K, byrow = TRUE)
  z_hat <- integer(n)
  for (i in seq_len(n)) z_hat[i] <- sample(0:(K - 1), 1, prob = Pi[, z[i] + 1])

  x <- cbind(1, rnorm(n))
  gamma <- c(0, 1.0, -0.9)
  alpha <- c(0.8, -0.7)
  eta <- as.numeric(x %*% alpha) + gamma[z + 1]
  y <- rpois(n, exp(eta))

  fit <- mcglm(y, z_hat = z_hat, x = x, family = "poisson",
               method = c("naive", "bca", "bcm", "cs"),
               Pi = Pi, pi_z = pi_z)

  expect_s3_class(fit, "mcglm")
  expect_equal(fit$K, K)
  expect_equal(fit$p, K - 1 + ncol(x))
  expect_named(fit$coefficients, c("naive", "bca", "bcm", "cs"))
})

# --- Subset methods ---

test_that("mcglm accepts subset of methods", {
  set.seed(4)
  n <- 200
  x <- cbind(1, rnorm(n))
  z_hat <- rbinom(n, 1, 0.5)
  y <- rpois(n, exp(0.5 * z_hat + x %*% c(-0.3, 0.4)))

  fit <- mcglm(y, z_hat = z_hat, x = x, family = "poisson",
               method = c("naive", "cs"),
               p01 = 0.1, p10 = 0.1, pi_z = 0.5)

  expect_named(fit$coefficients, c("naive", "cs"))
})

# --- S3 methods ---

test_that("print.mcglm runs without error", {
  set.seed(5)
  n <- 200
  x <- cbind(1, rnorm(n))
  z_hat <- rbinom(n, 1, 0.5)
  y <- rpois(n, exp(0.3 * z_hat + x %*% c(-0.2, 0.4)))

  fit <- mcglm(y, z_hat = z_hat, x = x, family = "poisson",
               method = c("naive", "bca"),
               p01 = 0.1, p10 = 0.1, pi_z = 0.5)

  expect_output(print(fit), "Coefficients")
})

test_that("summary.mcglm runs without error", {
  set.seed(5)
  n <- 200
  x <- cbind(1, rnorm(n))
  z_hat <- rbinom(n, 1, 0.5)
  y <- rpois(n, exp(0.3 * z_hat + x %*% c(-0.2, 0.4)))

  fit <- mcglm(y, z_hat = z_hat, x = x, family = "poisson",
               method = c("naive", "bca"),
               p01 = 0.1, p10 = 0.1, pi_z = 0.5)

  expect_output(print(summary(fit)), "Bias correction")
})

test_that("coef.mcglm returns all or selected method", {
  set.seed(6)
  n <- 200
  x <- cbind(1, rnorm(n))
  z_hat <- rbinom(n, 1, 0.5)
  y <- rpois(n, exp(0.3 * z_hat + x %*% c(-0.2, 0.4)))

  fit <- mcglm(y, z_hat = z_hat, x = x, family = "poisson",
               method = c("naive", "bca", "bcm"),
               p01 = 0.1, p10 = 0.1, pi_z = 0.5)

  all_coefs <- coef(fit)
  expect_type(all_coefs, "list")
  expect_named(all_coefs, c("naive", "bca", "bcm"))

  bca_coef <- coef(fit, method = "bca")
  expect_true(is.numeric(bca_coef))
  expect_length(bca_coef, 3)
})

# --- Iterate flag ---

test_that("iterate flag changes BCA/BCM results", {
  set.seed(7)
  n <- 1000
  x <- cbind(1, rnorm(n))
  z <- rbinom(n, 1, 0.4)
  eta <- 0.8 * z - 0.5 * x[, 1] + 0.7 * x[, 2]
  y <- rpois(n, exp(eta))
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, 0.10)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)

  fit1 <- mcglm(y, z_hat = z_hat, x = x, family = "poisson",
                method = c("naive", "bca", "bcm"),
                p01 = 0.10, p10 = 0.15, pi_z = 0.4,
                iterate = FALSE)
  fit2 <- mcglm(y, z_hat = z_hat, x = x, family = "poisson",
                method = c("naive", "bca", "bcm"),
                p01 = 0.10, p10 = 0.15, pi_z = 0.4,
                iterate = TRUE)

  expect_equal(fit1$coefficients$naive, fit2$coefficients$naive)
  expect_false(isTRUE(all.equal(fit1$coefficients$bca, fit2$coefficients$bca)))
})

# --- Parameter naming ---

test_that("binary parameter names are gamma, alpha0, alpha1, ...", {
  set.seed(8)
  n <- 200
  x <- cbind(1, rnorm(n), rnorm(n))
  z_hat <- rbinom(n, 1, 0.5)
  y <- rpois(n, exp(0.3 * z_hat + x %*% c(-0.2, 0.4, 0.1)))

  fit <- mcglm(y, z_hat = z_hat, x = x, family = "poisson", method = "naive")
  expect_equal(names(fit$coefficients$naive),
               c("gamma", "alpha0", "alpha1", "alpha2"))
})

test_that("multicategory parameter names are gamma1, gamma2, ..., alpha0, ...", {
  set.seed(9)
  n <- 300
  K <- 4
  z_hat <- sample(0:(K - 1), n, replace = TRUE)
  x <- cbind(1, rnorm(n))
  y <- rpois(n, 2)

  Pi <- matrix(0.05, K, K)
  diag(Pi) <- 0.85
  pi_z <- rep(1 / K, K)

  fit <- mcglm(y, z_hat = z_hat, x = x, family = "poisson",
               method = "naive", Pi = Pi, pi_z = pi_z)
  expect_equal(names(fit$coefficients$naive),
               c("gamma1", "gamma2", "gamma3", "alpha0", "alpha1"))
})

# --- Bias reduction (binary Poisson) ---

test_that("BCA/BCM/CS reduce bias for binary Poisson", {
  set.seed(42)
  n <- 5000
  psi0 <- c(0.8, -0.5, 0.7)
  x <- cbind(1, rnorm(n))
  z <- rbinom(n, 1, 0.4)
  y <- rpois(n, exp(psi0[1] * z + x %*% psi0[-1]))

  p01 <- 0.10; p10 <- 0.15
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)

  fit <- mcglm(y, z_hat = z_hat, x = x, family = "poisson",
               method = c("naive", "bca", "bcm", "cs"),
               p01 = p01, p10 = p10, pi_z = 0.4)

  bias_naive <- abs(fit$coefficients$naive[1] - psi0[1])
  bias_bca   <- abs(fit$coefficients$bca[1] - psi0[1])
  bias_cs    <- abs(fit$coefficients$cs[1] - psi0[1])

  expect_lt(bias_bca, bias_naive)
  expect_lt(bias_cs, bias_naive)
})

# --- Binomial family ---

test_that("mcglm works with binomial family", {
  set.seed(99)
  n <- 1000
  x <- cbind(1, rnorm(n))
  z <- rbinom(n, 1, 0.4)
  y <- rbinom(n, 1, plogis(0.8 * z - 0.5 * x[, 1] + 0.7 * x[, 2]))

  p01 <- 0.08; p10 <- 0.12
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)

  fit <- mcglm(y, z_hat = z_hat, x = x, family = "binomial",
               method = c("naive", "bca", "bcm", "cs"),
               p01 = p01, p10 = p10, pi_z = 0.4)

  expect_s3_class(fit, "mcglm")
  for (nm in names(fit$coefficients)) {
    expect_true(all(is.finite(fit$coefficients[[nm]])),
                info = paste("method:", nm))
  }
})

# --- Gaussian family ---

test_that("mcglm works with gaussian family", {
  set.seed(33)
  n <- 1000
  x <- cbind(1, rnorm(n))
  z <- rbinom(n, 1, 0.4)
  y <- 0.8 * z + x %*% c(-0.5, 0.7) + rnorm(n, 0, 0.5)

  p01 <- 0.10; p10 <- 0.10
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)

  fit <- mcglm(y, z_hat = z_hat, x = x, family = "gaussian",
               method = c("naive", "bca", "bcm", "cs"),
               p01 = p01, p10 = p10, pi_z = 0.4)

  expect_s3_class(fit, "mcglm")
  for (nm in names(fit$coefficients)) {
    expect_true(all(is.finite(fit$coefficients[[nm]])),
                info = paste("method:", nm))
  }
})

# --- Frequency weights ---

test_that("unit freq_weights match unweighted fit", {
  set.seed(10)
  n <- 300
  x <- cbind(1, rnorm(n))
  z_hat <- rbinom(n, 1, 0.5)
  y <- rpois(n, exp(0.5 * z_hat + x %*% c(-0.3, 0.4)))

  fit1 <- mcglm(y, z_hat = z_hat, x = x, family = "poisson",
                method = c("naive", "bca"),
                p01 = 0.1, p10 = 0.1, pi_z = 0.5)
  fit2 <- mcglm(y, z_hat = z_hat, x = x, family = "poisson",
                method = c("naive", "bca"),
                p01 = 0.1, p10 = 0.1, pi_z = 0.5,
                weights = rep(1, n))

  expect_equal(fit1$coefficients$naive, fit2$coefficients$naive, tolerance = 1e-8)
  expect_equal(fit1$coefficients$bca, fit2$coefficients$bca, tolerance = 1e-8)
})

# --- One-step estimator (if RTMB available) ---

# --- Inference & glm-style methods ---

test_that("mcglm stores per-method vcov and SE", {
  set.seed(20)
  n <- 800
  x <- cbind(1, rnorm(n))
  z <- rbinom(n, 1, 0.4)
  y <- rpois(n, exp(0.8 * z - 0.5 * x[, 1] + 0.7 * x[, 2]))
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, 0.10)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)

  fit <- mcglm(y, z_hat = z_hat, x = x, family = "poisson",
               method = c("naive", "bca", "bcm", "cs"),
               p01 = 0.10, p10 = 0.15, pi_z = 0.4)

  for (nm in c("naive", "bca", "bcm", "cs")) {
    V <- vcov(fit, method = nm)
    expect_true(is.matrix(V), info = nm)
    expect_equal(dim(V), c(3L, 3L), info = nm)
    expect_true(all(is.finite(V)), info = nm)
    expect_true(all(diag(V) >= 0), info = nm)
    expect_equal(rownames(V), names(coef(fit, method = nm)), info = nm)
  }
  expect_equal(length(fit$se$naive), 3L)
  expect_true(all(fit$se$cs > 0))
})

test_that("summary.mcglm prints a Wald table per method", {
  set.seed(21)
  n <- 500
  x <- cbind(1, rnorm(n))
  z_hat <- rbinom(n, 1, 0.5)
  y <- rpois(n, exp(0.3 * z_hat + x %*% c(-0.2, 0.4)))

  fit <- mcglm(y, z_hat = z_hat, x = x, family = "poisson",
               method = c("naive", "bca", "cs"),
               p01 = 0.1, p10 = 0.1, pi_z = 0.5)
  s <- summary(fit)

  expect_s3_class(s, "summary.mcglm")
  expect_named(s$coefficients, c("naive", "bca", "cs"))
  for (nm in names(s$coefficients)) {
    tab <- s$coefficients[[nm]]
    expect_equal(colnames(tab),
                 c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
    expect_true(all(is.finite(tab)))
  }
  expect_output(print(s), "Pr\\(>\\|z\\|\\)")
})

test_that("confint.mcglm returns Wald intervals", {
  set.seed(22)
  n <- 500
  x <- cbind(1, rnorm(n))
  z_hat <- rbinom(n, 1, 0.5)
  y <- rpois(n, exp(0.3 * z_hat + x %*% c(-0.2, 0.4)))

  fit <- mcglm(y, z_hat = z_hat, x = x, family = "poisson",
               method = c("naive", "bca"),
               p01 = 0.1, p10 = 0.1, pi_z = 0.5)
  ci <- confint(fit, method = "bca", level = 0.9)
  expect_equal(dim(ci), c(3L, 2L))
  expect_equal(rownames(ci), c("gamma", "alpha0", "alpha1"))
  expect_true(all(ci[, 1] < ci[, 2]))
  est <- coef(fit, method = "bca")
  expect_true(all(ci[, 1] <= est & est <= ci[, 2]))
})

test_that("fitted/predict/residuals work for mcglm", {
  set.seed(23)
  n <- 400
  x <- cbind(1, rnorm(n))
  z_hat <- rbinom(n, 1, 0.5)
  y <- rpois(n, exp(0.3 * z_hat + x %*% c(-0.2, 0.4)))

  fit <- mcglm(y, z_hat = z_hat, x = x, family = "poisson",
               method = c("naive", "bca"),
               p01 = 0.1, p10 = 0.1, pi_z = 0.5)

  fv <- fitted(fit, method = "naive")
  expect_length(fv, n)
  expect_true(all(fv > 0))

  pl <- predict(fit, type = "link", method = "bca")
  pr <- predict(fit, type = "response", method = "bca")
  expect_equal(exp(pl), unname(pr), tolerance = 1e-10)

  rr <- residuals(fit, method = "naive", type = "response")
  expect_equal(unname(rr), y - unname(fv), tolerance = 1e-10)

  rp <- residuals(fit, method = "naive", type = "pearson")
  expect_equal(rp, rr / sqrt(unname(fv)), tolerance = 1e-10)

  expect_silent(residuals(fit, method = "naive", type = "deviance"))
})

test_that("logLik / AIC / BIC / nobs / family / model.matrix work", {
  set.seed(24)
  n <- 300
  x <- cbind(1, rnorm(n))
  z_hat <- rbinom(n, 1, 0.5)
  y <- rpois(n, exp(0.3 * z_hat + x %*% c(-0.2, 0.4)))

  fit <- mcglm(y, z_hat = z_hat, x = x, family = "poisson",
               method = c("naive", "bca"),
               p01 = 0.1, p10 = 0.1, pi_z = 0.5)

  ll <- logLik(fit, method = "naive")
  expect_s3_class(ll, "logLik")
  expect_equal(attr(ll, "nobs"), n)

  expect_equal(AIC(fit), AIC(fit$naive_fit))
  expect_true(is.finite(BIC(fit)))
  expect_equal(nobs(fit), n)
  expect_s3_class(family(fit), "family")
  expect_equal(family(fit)$family, "poisson")

  M <- model.matrix(fit)
  expect_equal(dim(M), c(n, 3L))
})

test_that("Wald CI coverage is reasonable for naive Poisson", {
  skip_on_cran()
  set.seed(101)
  n <- 1500
  psi0 <- c(0.6, -0.3, 0.5)
  hits <- 0L
  reps <- 50
  for (i in seq_len(reps)) {
    x <- cbind(1, rnorm(n))
    z <- rbinom(n, 1, 0.4)
    y <- rpois(n, exp(psi0[1] * z + x %*% psi0[-1]))
    p01 <- 0.05; p10 <- 0.05
    z_hat <- z
    z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01)
    z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)
    fit <- mcglm(y, z_hat = z_hat, x = x, family = "poisson",
                 method = c("naive", "bca"),
                 p01 = p01, p10 = p10, pi_z = 0.4)
    ci <- confint(fit, method = "bca")
    if (psi0[1] >= ci[1, 1] && psi0[1] <= ci[1, 2]) hits <- hits + 1L
  }
  expect_gt(hits / reps, 0.7)
})

test_that("mcglm onestep works when RTMB available", {
  skip_if_not_installed("RTMB")

  set.seed(42)
  n <- 1000
  x <- cbind(1, rnorm(n))
  z <- rbinom(n, 1, 0.4)
  y <- rpois(n, exp(0.8 * z - 0.5 * x[, 1] + 0.7 * x[, 2]))

  p01 <- 0.10; p10 <- 0.15
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)

  fit <- mcglm(y, z_hat = z_hat, x = x, family = "poisson",
               method = c("naive", "onestep"),
               p01 = p01, p10 = p10, pi_z = 0.4,
               fix_omega = TRUE)

  expect_named(fit$coefficients, c("naive", "onestep"))
  expect_true(all(is.finite(fit$coefficients$onestep)))
  expect_false(is.null(fit$vcov_onestep))
  expect_equal(nrow(fit$vcov_onestep), 3)
  # vcov.mcglm dispatch should also surface the onestep matrix
  expect_equal(vcov(fit, method = "onestep"), fit$vcov_onestep)
  expect_true(all(se.mcglm(fit, method = "onestep") > 0))
  ll <- logLik(fit, method = "onestep")
  expect_s3_class(ll, "logLik")
})
