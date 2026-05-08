# Tests for the AKN98 corrected-score method (cs_akn).

# ---- helpers ----

gen_binary <- function(n = 4000, gamma = 0.8, alpha = c(-0.5, 0.7),
                       p01 = 0.1, p10 = 0.1, prev = 0.4, family,
                       seed = 42) {
  set.seed(seed)
  Pi <- matrix(c(1 - p01, p01, p10, 1 - p10), 2, 2)
  z  <- rbinom(n, 1, prev)
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)
  x1  <- rnorm(n)
  eta <- gamma * z + alpha[1] + alpha[2] * x1
  y <- switch(family,
              poisson  = rpois(n, exp(eta)),
              binomial = rbinom(n, 1, plogis(eta)),
              gaussian = eta + rnorm(n, 0, 0.4))
  list(df = data.frame(y = y, z = factor(z_hat), x1 = x1),
       Pi = Pi, truth = c(gamma, alpha))
}

gen_multi <- function(n = 4000, K = 3,
                      gamma = c(1.0, -0.9), alpha = c(0.8, -0.7),
                      seed = 7) {
  set.seed(seed)
  Pi <- matrix(0.1, K, K); diag(Pi) <- 0.8
  pi_z <- c(0.5, 0.3, 0.2)
  Z  <- sample(0:(K - 1), n, replace = TRUE, prob = pi_z)
  Zh <- vapply(Z, function(z) sample(0:(K - 1), 1, prob = Pi[, z + 1]), 0L)
  x1 <- rnorm(n)
  eta <- c(0, gamma)[Z + 1] + alpha[1] + alpha[2] * x1
  y  <- rpois(n, exp(eta))
  list(y = y, z_hat = Zh, x = cbind(1, x1),
       Pi = Pi, pi_z = pi_z, truth = c(gamma, alpha))
}

# ---- basic correctness ----

test_that("cs_akn is close to truth for binary Poisson", {
  d <- gen_binary(n = 4000, family = "poisson", seed = 42)
  fit <- mcglm(y ~ mc(z, d$Pi) + x1, data = d$df, family = "poisson",
               method = c("naive", "cs", "cs_akn"))
  cf <- coef(fit, method = "cs_akn")
  expect_true(all(abs(cf - d$truth) < 0.07))
  # cs and cs_akn solve different forms of the same equation -> close.
  cf_drift <- coef(fit, method = "cs")
  expect_true(max(abs(cf - cf_drift)) < 0.05)
})

test_that("cs_akn is close to truth for binary binomial", {
  d <- gen_binary(n = 6000, family = "binomial", seed = 11,
                  gamma = 0.9, alpha = c(-0.4, 0.6), p01 = 0.15, p10 = 0.10)
  fit <- mcglm(y ~ mc(z, d$Pi) + x1, data = d$df, family = "binomial",
               method = c("naive", "cs", "cs_akn"))
  cf <- coef(fit, method = "cs_akn")
  expect_true(all(abs(cf - d$truth) < 0.10))
})

test_that("cs_akn closed-form Gaussian is consistent with the drift CS", {
  d <- gen_binary(n = 4000, family = "gaussian", seed = 33)
  fit <- mcglm(y ~ mc(z, d$Pi) + x1, data = d$df, family = "gaussian",
               method = c("cs", "cs_akn"))
  cf_akn   <- coef(fit, method = "cs_akn")
  cf_drift <- coef(fit, method = "cs")
  # Drift CS and AKN98 surrogate solve different unbiased equations
  # (different finite-sample iterates) but share a population limit.
  # On n = 4000 they agree to a few parts in 1e-4.
  expect_true(max(abs(cf_akn - cf_drift)) < 1e-3)
  expect_true(all(abs(cf_akn - d$truth) < 0.05))
})

test_that("cs_akn works for K=3 Poisson with matrix interface", {
  d <- gen_multi(n = 4000, K = 3, seed = 7)
  fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "poisson",
               method = c("naive", "cs", "cs_akn"),
               Pi = d$Pi, pi_z = d$pi_z)
  cf <- coef(fit, method = "cs_akn")
  expect_named(cf, c("gamma1", "gamma2", "alpha0", "alpha1"))
  expect_true(all(abs(cf - d$truth) < 0.25))
  # cs_akn should agree closely with the drift cs at this n.
  expect_true(max(abs(cf - coef(fit, method = "cs"))) < 0.05)
})

# ---- variance ----

test_that("cs_akn vcov is positive semidefinite and SE is finite", {
  d <- gen_binary(n = 3000, family = "poisson", seed = 5)
  fit <- mcglm(y ~ mc(z, d$Pi) + x1, data = d$df, family = "poisson",
               method = "cs_akn")
  V <- vcov(fit, method = "cs_akn")
  expect_true(all(is.finite(V)))
  expect_true(all(diag(V) > 0))
  ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  expect_true(min(ev) > -1e-10)
  se <- sqrt(diag(V))
  expect_true(all(is.finite(se)))
  expect_named(se, c("gamma", "alpha0", "alpha1"))
})

test_that("cs_akn matrix-interface accepts only Pi (no pi_z required)", {
  d <- gen_multi(n = 2000, K = 3, seed = 21)
  # cs_akn does not need pi_z; passing only Pi must work. The package
  # always returns the naive fit alongside requested methods.
  fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "poisson",
               method = "cs_akn", Pi = d$Pi)
  expect_s3_class(fit, "mcglm")
  expect_true("cs_akn" %in% names(fit$coefficients))
  expect_true(all(is.finite(coef(fit, method = "cs_akn"))))
})

# ---- input validation ----

test_that("cs_akn errors when neither Pi nor (p01,p10) is supplied (binary)", {
  d <- gen_binary(n = 1000, family = "poisson", seed = 1)
  expect_error(
    mcglm(d$df$y, z_hat = as.integer(as.character(d$df$z)),
          x = cbind(1, d$df$x1), family = "poisson", method = "cs_akn"),
    "supply Pi"
  )
})

test_that("cs_akn errors when Pi is missing (multicategory)", {
  d <- gen_multi(n = 1000, K = 3, seed = 2)
  expect_error(
    mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "poisson",
          method = "cs_akn"),
    "supply Pi"
  )
})

test_that("cs_akn rejects multinomial responses", {
  set.seed(8)
  n <- 500
  Pi <- matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2)
  z <- rbinom(n, 1, 0.4); z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, 0.1)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.1)
  x1 <- rnorm(n)
  y <- sample(0:2, n, replace = TRUE)
  expect_error(
    mcglm(y, z_hat = as.integer(z_hat), x = cbind(1, x1),
          family = "multinomial", method = "cs_akn", Pi = Pi),
    "multinomial"
  )
})

test_that("cs_akn helper reports non-identification clearly", {
  expect_error(
    mismeasured:::.akn_build_Q(diag(2), K = 3),
    "K x K"
  )
  expect_warning(
    expect_error(mismeasured:::.akn_build_Q(matrix(0.5, 2, 2), K = 2),
                 "singular"),
    "near-singular"
  )
})

# ---- weights ----

test_that("cs_akn handles frequency weights", {
  d  <- gen_binary(n = 2000, family = "poisson", seed = 9)
  wt <- runif(nrow(d$df), 0.5, 1.5)
  fit <- mcglm(y ~ mc(z, d$Pi) + x1, data = d$df, family = "poisson",
               method = "cs_akn", weights = wt)
  cf <- coef(fit, method = "cs_akn")
  expect_true(all(is.finite(cf)))
  V <- vcov(fit, method = "cs_akn")
  expect_true(all(is.finite(V)))
})

test_that("cs_akn weighted Gaussian uses the closed-form branch", {
  set.seed(11)
  n <- 200
  Pi <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
  z <- rbinom(n, 1, 0.4)
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, 0.1)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)
  x1 <- rnorm(n)
  y <- 0.6 * z - 0.3 + 0.4 * x1 + rnorm(n, sd = 0.5)
  w <- runif(n, 0.5, 1.5)

  fit <- mcglm(y, z_hat = z_hat, x = cbind(1, x1), family = "gaussian",
               method = "cs_akn", Pi = Pi, weights = w)
  expect_true(all(is.finite(coef(fit, method = "cs_akn"))))
})

# ---- agreement with cs as n grows ----

test_that("cs_akn and cs converge to the same point as n grows (Poisson K=2)", {
  d <- gen_binary(n = 20000, family = "poisson", seed = 4242)
  fit <- mcglm(y ~ mc(z, d$Pi) + x1, data = d$df, family = "poisson",
               method = c("cs", "cs_akn"))
  expect_true(max(abs(coef(fit, method = "cs") -
                      coef(fit, method = "cs_akn"))) < 0.01)
})
