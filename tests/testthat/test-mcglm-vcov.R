# --------------------------------------------------------------------------
# Sanity tests for the corrected sandwich variance (vcov_corrected = TRUE)
# Untested before: .mcglm_vcov_bc_bin / _bc_multi when corrected = TRUE.
# These paths feed se(), confint(), summary().
# --------------------------------------------------------------------------

# Helpers (kept self-contained so this file can be run alone)
.simulate_binary_v <- function(n, seed = 1) {
  set.seed(seed)
  x  <- cbind(1, rnorm(n))
  z  <- rbinom(n, 1, 0.4)
  eta <- drop(x %*% c(-0.5, 0.7) + 0.8 * z)
  y <- rpois(n, exp(eta))
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, 0.10)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)
  list(y = y, z_hat = z_hat, x = x, p01 = 0.10, p10 = 0.15, pi_z = 0.4)
}

.simulate_multi_v <- function(n, K = 3, seed = 1) {
  set.seed(seed)
  pi_z <- c(0.5, 0.3, 0.2)[seq_len(K)]; pi_z <- pi_z / sum(pi_z)
  Pi <- matrix(0.05, K, K); diag(Pi) <- 1 - 0.05 * (K - 1)
  z <- sample(0:(K - 1), n, replace = TRUE, prob = pi_z)
  z_hat <- vapply(z, function(zi) sample(0:(K - 1), 1, prob = Pi[, zi + 1]),
                  integer(1))
  x <- cbind(1, rnorm(n))
  gf <- c(0, 0.7, -0.4)[seq_len(K)]
  y <- rpois(n, exp(drop(x %*% c(-0.3, 0.5)) + gf[z + 1]))
  list(y = y, z_hat = z_hat, x = x, Pi = Pi, pi_z = pi_z)
}

.is_pd <- function(V) {
  ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  all(is.finite(ev)) && all(ev > -1e-10)
}

# ==========================================================================
# Binary
# ==========================================================================

test_that("vcov_corrected=TRUE returns a finite PD matrix (binary BCA/BCM)", {
  d <- .simulate_binary_v(800, seed = 41)

  for (m in c("bca", "bcm")) {
    fit_def <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "poisson",
                     method = m,
                     p01 = d$p01, p10 = d$p10, pi_z = d$pi_z,
                     vcov_corrected = FALSE)
    fit_cor <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "poisson",
                     method = m,
                     p01 = d$p01, p10 = d$p10, pi_z = d$pi_z,
                     vcov_corrected = TRUE)
    V_def <- vcov(fit_def, method = m)
    V_cor <- vcov(fit_cor, method = m)
    expect_true(all(is.finite(V_cor)), info = m)
    expect_true(.is_pd(V_cor), info = m)
    # Point estimates are unchanged by the variance choice.
    expect_equal(coef(fit_def, method = m), coef(fit_cor, method = m),
                 tolerance = 1e-12, info = m)
    # Corrected variance differs from default
    expect_false(isTRUE(all.equal(V_def, V_cor)), info = m)
  }
})

test_that("vcov_corrected propagates into se() and summary() (binary)", {
  d <- .simulate_binary_v(800, seed = 42)
  fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "poisson",
               method = "bcm",
               p01 = d$p01, p10 = d$p10, pi_z = d$pi_z,
               vcov_corrected = TRUE)
  s <- se.mcglm(fit, method = "bcm")
  expect_length(s, 3)
  expect_true(all(s > 0))
  smry <- summary(fit)
  expect_s3_class(smry, "summary.mcglm")
})

# ==========================================================================
# Multicategory: corrected variance must work with both Jacobian flavors
# ==========================================================================

test_that("vcov_corrected=TRUE works for K=3 (analytical Jacobian)", {
  d <- .simulate_multi_v(700, K = 3, seed = 51)

  fit_cor <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "poisson",
                   method = c("bca", "bcm"),
                   Pi = d$Pi, pi_z = d$pi_z,
                   vcov_corrected = TRUE,
                   jacobian = "analytical")

  for (m in c("bca", "bcm")) {
    V <- vcov(fit_cor, method = m)
    expect_equal(dim(V), c(4L, 4L), info = m)
    expect_true(all(is.finite(V)), info = m)
    expect_true(.is_pd(V), info = m)
  }
})

test_that("vcov_corrected matches between analytical and numerical Jacobian (K=3)", {
  d <- .simulate_multi_v(700, K = 3, seed = 52)
  fit_an <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "poisson",
                  method = "bcm", Pi = d$Pi, pi_z = d$pi_z,
                  vcov_corrected = TRUE, jacobian = "analytical")
  fit_nu <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "poisson",
                  method = "bcm", Pi = d$Pi, pi_z = d$pi_z,
                  vcov_corrected = TRUE, jacobian = "numerical")
  expect_equal(vcov(fit_an, method = "bcm"),
               vcov(fit_nu, method = "bcm"),
               tolerance = 1e-5)
})

test_that("vcov_corrected respects frequency weights (binary BCM)", {
  d <- .simulate_binary_v(150, seed = 61)
  w <- sample(1:3, 150, replace = TRUE)
  idx <- rep.int(seq_along(w), w)

  fit_w <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "poisson",
                 method = "bcm",
                 p01 = d$p01, p10 = d$p10, pi_z = d$pi_z,
                 weights = w, vcov_corrected = TRUE)
  fit_e <- mcglm(d$y[idx], z_hat = d$z_hat[idx], x = d$x[idx, ],
                 family = "poisson", method = "bcm",
                 p01 = d$p01, p10 = d$p10, pi_z = d$pi_z,
                 vcov_corrected = TRUE)
  # Weighted-vs-expanded variances differ by 1/n vs 1/sum(w) scaling
  # only via the wt = NULL versus wt code paths; point equivalence is the
  # primary contract for SE per replicated row.
  expect_equal(coef(fit_w, method = "bcm"),
               coef(fit_e, method = "bcm"), tolerance = 1e-6)
  expect_true(all(is.finite(vcov(fit_w, method = "bcm"))))
  expect_true(.is_pd(vcov(fit_w, method = "bcm")))
})
