# Regression pins for the improved MC-SIMEX variance machinery.

.make_mc_df <- function(n = 500, seed = 1) {
  set.seed(seed)
  z <- rbinom(n, 1, 0.4)
  y <- 1 + 2 * z + rnorm(n)
  z_star <- ifelse(z == 1, 1 - rbinom(n, 1, 0.12), rbinom(n, 1, 0.1))
  data.frame(y = y, z = factor(z_star))
}

test_that("improved MC-SIMEX vcov works with multiple lambda incl. 0", {
  Pi <- matrix(c(0.9, 0.1, 0.12, 0.88), 2, 2)
  df <- .make_mc_df()

  fit <- simex(y ~ mc(z, Pi), data = df, lambda = c(0, 0.5, 1), B = 5,
               seed = 42)

  expect_identical(dim(fit$vcov), c(2L, 2L))
  expect_true(all(is.finite(fit$vcov)))
  expect_equal(fit$vcov, t(fit$vcov))
  expect_true(all(eigen(fit$vcov, symmetric = TRUE,
                        only.values = TRUE)$values >= 0))
})

test_that("improved MC-SIMEX SEs do not shrink as B grows", {
  Pi <- matrix(c(0.9, 0.1, 0.12, 0.88), 2, 2)
  df <- .make_mc_df()

  se_small <- sqrt(diag(simex(y ~ mc(z, Pi), data = df, B = 10,
                              seed = 42)$vcov))
  se_large <- sqrt(diag(simex(y ~ mc(z, Pi), data = df, B = 200,
                              seed = 42)$vcov))

  # The pre-fix estimator returned MC error of the mean, so the ratio
  # would be near sqrt(200 / 10) ~ 4.5; a calibrated SE is stable in B.
  expect_true(all(se_small / se_large > 0.7))
  expect_true(all(se_small / se_large < 1.4))
})

test_that(".resample_z_lambda draws from the columns of Pi^lambda", {
  expect_identical(
    mismeasured:::.resample_z_lambda(c(0L, 1L, 2L), diag(3), 0, 3L),
    c(0L, 1L, 2L)
  )

  Pi <- matrix(c(0.8, 0.15, 0.05,
                 0.10, 0.80, 0.10,
                 0.05, 0.15, 0.80), 3, 3)
  lambda <- 0.8
  Pi_lam <- mismeasured:::.mat_power_r(Pi, lambda)

  set.seed(7)
  n_per <- 20000L
  z_hat <- rep(0:2, each = n_per)
  z_star <- mismeasured:::.resample_z_lambda(z_hat, Pi, lambda, 3L)

  expect_true(all(z_star %in% 0:2))
  for (j in 1:3) {
    freq <- tabulate(z_star[z_hat == j - 1L] + 1L, nbins = 3L) / n_per
    expect_equal(freq, Pi_lam[, j], tolerance = 0.02)
  }
})

test_that(".variance_k_improved clips non-PSD results to PSD", {
  set.seed(11)
  # Between-replicate spread far larger than the sampling variance forces
  # V_sampling - S_mc * (1 - 1/B) negative before the PSD guard.
  theta_list <- list(matrix(rnorm(200, sd = 5), ncol = 2))
  V <- mismeasured:::.variance_k_improved(
    theta_list = theta_list,
    vcov_list = list(diag(2) * 1e-6),
    transform_list = list(diag(2)),
    lambda = 1,
    B = 100,
    p = 2
  )
  ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev >= 0))
  expect_true(all(is.finite(V)))
})
