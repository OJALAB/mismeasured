# --------------------------------------------------------------------------
# mcglm() onestep estimator: paths previously uncovered.
#  * Multicategory K=3 onestep (.mcglm_fit_onestep_multi)
#  * fix_omega = FALSE (default; mixture weights estimated jointly)
#  * Gaussian family with homoskedastic = TRUE/FALSE
# --------------------------------------------------------------------------

.simulate_K_os <- function(n, K = 3, family = "poisson", seed = 1) {
  set.seed(seed)
  pi_z <- c(0.5, 0.3, 0.2)[seq_len(K)]; pi_z <- pi_z / sum(pi_z)
  Pi <- matrix(0.05, K, K); diag(Pi) <- 1 - 0.05 * (K - 1)
  z <- sample(0:(K - 1), n, replace = TRUE, prob = pi_z)
  z_hat <- vapply(z, function(zi) sample(0:(K - 1), 1, prob = Pi[, zi + 1]),
                  integer(1))
  x <- cbind(1, rnorm(n))
  gf <- c(0, 0.7, -0.4)[seq_len(K)]
  eta <- drop(x %*% c(-0.3, 0.5)) + gf[z + 1]
  y <- switch(family,
              poisson  = rpois(n, exp(eta)),
              binomial = rbinom(n, 1, plogis(eta)),
              gaussian = eta + rnorm(n, sd = 0.5))
  list(y = y, z_hat = z_hat, x = x, Pi = Pi, pi_z = pi_z)
}

.simulate_bin_os <- function(n, family = "gaussian", seed = 1) {
  set.seed(seed)
  x  <- cbind(1, rnorm(n))
  z  <- rbinom(n, 1, 0.4)
  eta <- drop(x %*% c(-0.5, 0.7) + 0.8 * z)
  y <- switch(family,
              gaussian = eta + rnorm(n, sd = 0.5),
              poisson  = rpois(n, exp(eta)),
              binomial = rbinom(n, 1, plogis(eta)))
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, 0.10)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)
  list(y = y, z_hat = z_hat, x = x, p01 = 0.10, p10 = 0.15, pi_z = 0.4)
}

# ==========================================================================
# Multicategory onestep
# ==========================================================================

test_that("K=3 onestep runs with fix_omega=TRUE (Poisson)", {
  skip_if_not_installed("RTMB")
  d <- .simulate_K_os(500, K = 3, family = "poisson", seed = 81)

  fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "poisson",
               method = c("naive", "onestep"),
               Pi = d$Pi, pi_z = d$pi_z, fix_omega = TRUE)
  expect_named(fit$coefficients, c("naive", "onestep"))
  os <- fit$coefficients$onestep
  expect_length(os, (3 - 1) + ncol(d$x))
  expect_true(all(is.finite(os)))
  V <- vcov(fit, method = "onestep")
  expect_equal(dim(V), rep(length(os), 2))
  expect_true(all(is.finite(V)))
})

test_that("K=3 onestep with fix_omega=FALSE estimates mixture jointly", {
  skip_if_not_installed("RTMB")
  d <- .simulate_K_os(700, K = 3, family = "poisson", seed = 82)

  # Joint estimation of K*K-1 mixture weights is harder; tolerate
  # non-convergence / mild Hessian warnings — what we are checking is
  # that the code path runs end-to-end and returns finite coefficients.
  fit <- suppressWarnings(
    mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "poisson",
          method = "onestep", Pi = d$Pi, pi_z = d$pi_z,
          fix_omega = FALSE)
  )
  os <- fit$coefficients$onestep
  expect_true(all(is.finite(os)))
  # The reported coefficient block has length (K-1) + r; the K*K - 1
  # mixture parameters live inside vcov_onestep, not in coefficients.
  expect_length(os, (3 - 1) + ncol(d$x))
})

# ==========================================================================
# fix_omega = FALSE for binary onestep (default; previously untested)
# ==========================================================================

test_that("binary onestep fix_omega=FALSE runs without supplying p01/p10", {
  skip_if_not_installed("RTMB")
  d <- .simulate_bin_os(400, family = "poisson", seed = 83)

  # No p01/p10/pi_z required when mixture weights are estimated.
  fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "poisson",
               method = "onestep", fix_omega = FALSE)
  expect_true(all(is.finite(fit$coefficients$onestep)))
  expect_length(fit$coefficients$onestep, 1 + ncol(d$x))
})

# ==========================================================================
# Gaussian onestep, homoskedastic on/off
# ==========================================================================

test_that("Gaussian onestep with homoskedastic=TRUE recovers parameters", {
  skip_if_not_installed("RTMB")
  d <- .simulate_bin_os(800, family = "gaussian", seed = 84)

  fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "gaussian",
               method = "onestep",
               p01 = d$p01, p10 = d$p10, pi_z = d$pi_z,
               fix_omega = TRUE, homoskedastic = TRUE)
  os <- fit$coefficients$onestep
  expect_true(all(is.finite(os)))
  # gamma should be in the right ballpark of 0.8
  expect_lt(abs(os[1] - 0.8), 0.3)
})

test_that("Gaussian onestep with homoskedastic=FALSE runs", {
  skip_if_not_installed("RTMB")
  d <- .simulate_bin_os(800, family = "gaussian", seed = 85)

  fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "gaussian",
               method = "onestep",
               p01 = d$p01, p10 = d$p10, pi_z = d$pi_z,
               fix_omega = TRUE, homoskedastic = FALSE)
  expect_true(all(is.finite(fit$coefficients$onestep)))
})

# ==========================================================================
# logLik / AIC / BIC are exposed for onestep
# ==========================================================================

test_that("logLik / AIC / BIC work for onestep fits", {
  skip_if_not_installed("RTMB")
  d <- .simulate_bin_os(500, family = "poisson", seed = 86)

  fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "poisson",
               method = "onestep",
               p01 = d$p01, p10 = d$p10, pi_z = d$pi_z, fix_omega = TRUE)
  ll <- logLik(fit, method = "onestep")
  expect_s3_class(ll, "logLik")
  expect_true(is.finite(AIC(fit, method = "onestep")))
  expect_true(is.finite(BIC(fit, method = "onestep")))
})
