# --------------------------------------------------------------------------
# mcglm() under frequency weights
# Strategy: duplication-equivalence — fitting on n rows with frequency
# weights w_i must equal fitting on the row-replicated dataset (each row
# repeated w_i times). Held to numerical precision.
# --------------------------------------------------------------------------

# Helper: replicate rows according to weights and refit unweighted
.expand_by_weights <- function(z, x, y, w) {
  idx <- rep.int(seq_along(w), w)
  list(z = z[idx], x = x[idx, , drop = FALSE], y = y[idx])
}

.simulate_binary <- function(n, family, p01 = 0.10, p10 = 0.15,
                             alpha = c(-0.5, 0.7), gamma = 0.8, pi_z = 0.4,
                             seed = 1) {
  set.seed(seed)
  x  <- cbind(1, rnorm(n))
  z  <- rbinom(n, 1, pi_z)
  eta <- drop(x %*% alpha + gamma * z)
  y <- switch(family,
              binomial = rbinom(n, 1, plogis(eta)),
              poisson  = rpois(n, exp(eta)),
              gaussian = eta + rnorm(n, sd = 0.5))
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)
  list(y = y, z_hat = z_hat, x = x, p01 = p01, p10 = p10, pi_z = pi_z)
}

.simulate_K <- function(n, K = 3, family = "poisson", seed = 1) {
  set.seed(seed)
  pi_vec <- c(0.5, 0.3, 0.2)[seq_len(K)]
  pi_vec <- pi_vec / sum(pi_vec)
  Pi <- matrix(0.05, K, K); diag(Pi) <- 1 - 0.05 * (K - 1)
  z <- sample(0:(K - 1), n, replace = TRUE, prob = pi_vec)
  z_hat <- vapply(z, function(zi) sample(0:(K - 1), 1, prob = Pi[, zi + 1]),
                  integer(1))
  x <- cbind(1, rnorm(n))
  gamma_full <- c(0, 0.7, -0.4)[seq_len(K)]
  alpha <- c(-0.3, 0.5)
  eta <- drop(x %*% alpha) + gamma_full[z + 1]
  y <- switch(family,
              binomial = rbinom(n, 1, plogis(eta)),
              poisson  = rpois(n, exp(eta)),
              gaussian = eta + rnorm(n, sd = 0.5))
  list(y = y, z_hat = z_hat, x = x, Pi = Pi, pi_z = pi_vec)
}

# ==========================================================================
# Binary case: bcm, cs match expansion under integer weights
# ==========================================================================

test_that("binary bcm/cs respect frequency weights (Poisson)", {
  d <- .simulate_binary(150, "poisson", seed = 11)
  w <- sample(1:3, 150, replace = TRUE)
  ex <- .expand_by_weights(d$z_hat, d$x, d$y, w)

  fit_w <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "poisson",
                 method = c("naive", "bca", "bcm", "cs"),
                 p01 = d$p01, p10 = d$p10, pi_z = d$pi_z,
                 weights = w)
  fit_e <- mcglm(ex$y, z_hat = ex$z, x = ex$x, family = "poisson",
                 method = c("naive", "bca", "bcm", "cs"),
                 p01 = d$p01, p10 = d$p10, pi_z = d$pi_z)

  for (m in c("naive", "bca", "bcm", "cs")) {
    expect_equal(fit_w$coefficients[[m]], fit_e$coefficients[[m]],
                 tolerance = 1e-6, info = m)
  }
})

test_that("binary bcm/cs respect frequency weights (binomial)", {
  d <- .simulate_binary(150, "binomial", seed = 12)
  w <- sample(1:3, 150, replace = TRUE)
  ex <- .expand_by_weights(d$z_hat, d$x, d$y, w)

  fit_w <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "binomial",
                 method = c("bcm", "cs"),
                 p01 = d$p01, p10 = d$p10, pi_z = d$pi_z, weights = w)
  fit_e <- mcglm(ex$y, z_hat = ex$z, x = ex$x, family = "binomial",
                 method = c("bcm", "cs"),
                 p01 = d$p01, p10 = d$p10, pi_z = d$pi_z)

  expect_equal(fit_w$coefficients$bcm, fit_e$coefficients$bcm, tolerance = 1e-6)
  expect_equal(fit_w$coefficients$cs,  fit_e$coefficients$cs,  tolerance = 1e-6)
})

# ==========================================================================
# Multicategory case: bca, bcm, cs (analytical & numerical Jacobian)
# ==========================================================================

test_that("multicategory bca/bcm/cs respect frequency weights (analytical)", {
  d <- .simulate_K(120, K = 3, family = "poisson", seed = 21)
  w <- sample(1:3, 120, replace = TRUE)
  ex <- .expand_by_weights(d$z_hat, d$x, d$y, w)

  fit_w <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "poisson",
                 method = c("naive", "bca", "bcm", "cs"),
                 Pi = d$Pi, pi_z = d$pi_z, weights = w,
                 jacobian = "analytical")
  fit_e <- mcglm(ex$y, z_hat = ex$z, x = ex$x, family = "poisson",
                 method = c("naive", "bca", "bcm", "cs"),
                 Pi = d$Pi, pi_z = d$pi_z,
                 jacobian = "analytical")
  for (m in c("naive", "bca", "bcm", "cs")) {
    expect_equal(fit_w$coefficients[[m]], fit_e$coefficients[[m]],
                 tolerance = 1e-6, info = m)
  }
})

test_that("multicategory bcm/cs respect frequency weights (numerical Jacobian)", {
  d <- .simulate_K(120, K = 3, family = "binomial", seed = 22)
  w <- sample(1:3, 120, replace = TRUE)
  ex <- .expand_by_weights(d$z_hat, d$x, d$y, w)

  fit_w <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "binomial",
                 method = c("bcm", "cs"),
                 Pi = d$Pi, pi_z = d$pi_z, weights = w,
                 jacobian = "numerical")
  fit_e <- mcglm(ex$y, z_hat = ex$z, x = ex$x, family = "binomial",
                 method = c("bcm", "cs"),
                 Pi = d$Pi, pi_z = d$pi_z,
                 jacobian = "numerical")
  expect_equal(fit_w$coefficients$bcm, fit_e$coefficients$bcm, tolerance = 1e-5)
  expect_equal(fit_w$coefficients$cs,  fit_e$coefficients$cs,  tolerance = 1e-5)
})

# ==========================================================================
# Onestep: K=2 and K=3 must respect weights too
# ==========================================================================

test_that("binary onestep respects frequency weights", {
  skip_if_not_installed("RTMB")
  d <- .simulate_binary(200, "poisson", seed = 31)
  w <- sample(1:2, 200, replace = TRUE)
  ex <- .expand_by_weights(d$z_hat, d$x, d$y, w)

  fit_w <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "poisson",
                 method = "onestep",
                 p01 = d$p01, p10 = d$p10, pi_z = d$pi_z,
                 fix_omega = TRUE, weights = w)
  fit_e <- mcglm(ex$y, z_hat = ex$z, x = ex$x, family = "poisson",
                 method = "onestep",
                 p01 = d$p01, p10 = d$p10, pi_z = d$pi_z,
                 fix_omega = TRUE)

  expect_equal(fit_w$coefficients$onestep, fit_e$coefficients$onestep,
               tolerance = 1e-3)
})

test_that("multicategory onestep respects frequency weights", {
  skip_if_not_installed("RTMB")
  d <- .simulate_K(150, K = 3, family = "poisson", seed = 32)
  w <- sample(1:2, 150, replace = TRUE)
  ex <- .expand_by_weights(d$z_hat, d$x, d$y, w)

  fit_w <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "poisson",
                 method = "onestep",
                 Pi = d$Pi, pi_z = d$pi_z, fix_omega = TRUE, weights = w)
  fit_e <- mcglm(ex$y, z_hat = ex$z, x = ex$x, family = "poisson",
                 method = "onestep",
                 Pi = d$Pi, pi_z = d$pi_z, fix_omega = TRUE)

  expect_equal(fit_w$coefficients$onestep, fit_e$coefficients$onestep,
               tolerance = 1e-3)
})
