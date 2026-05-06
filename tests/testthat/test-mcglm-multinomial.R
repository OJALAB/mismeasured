# --------------------------------------------------------------------------
# mcglm() multinomial response family
# Covers: naive (via nnet) and onestep (via RTMB) for binary K=2 and K=3
# proxy. Also asserts that fitted/predict/residuals/logLik raise clear
# errors for the multinomial response (not yet implemented).
# --------------------------------------------------------------------------

.simulate_multinom <- function(n, J = 3, K = 2, seed = 1) {
  set.seed(seed)
  x <- cbind(1, rnorm(n))
  pi_z <- rep(1 / K, K)
  Pi <- matrix(0.05, K, K); diag(Pi) <- 1 - 0.05 * (K - 1)
  z <- sample(0:(K - 1), n, replace = TRUE, prob = pi_z)
  z_hat <- vapply(z, function(zi) sample(0:(K - 1), 1, prob = Pi[, zi + 1]),
                  integer(1))

  # Generate y from a J-class multinomial logit with category-specific
  # gamma_j[z] and intercepts so the test data has all classes present.
  # eta_{i,j} = alpha0_j + alpha1_j * x[i,2] + gamma_j * 1{z=1}
  alpha0 <- c(0, 0.4, -0.3)[seq_len(J)]
  alpha1 <- c(0,  0.5,  0.2)[seq_len(J)]
  gam    <- c(0,  0.6, -0.5)[seq_len(J)]
  z_ind <- as.numeric(z == 1)
  eta <- outer(seq_len(n), seq_len(J),
               Vectorize(function(i, j) alpha0[j] + alpha1[j] * x[i, 2] +
                                          gam[j] * z_ind[i]))
  P <- exp(eta); P <- P / rowSums(P)
  y <- vapply(seq_len(n), function(i) sample(0:(J - 1), 1, prob = P[i, ]),
              integer(1))

  list(y = y, z_hat = z_hat, x = x, J = J, K = K, Pi = Pi, pi_z = pi_z)
}

# ==========================================================================
# Naive multinomial path
# ==========================================================================

test_that("multinomial naive fit returns sensible structure (K=2)", {
  skip_if_not_installed("nnet")
  d <- .simulate_multinom(600, J = 3, K = 2, seed = 71)

  fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "multinomial",
               method = "naive", K = d$K, J = d$J)

  expect_s3_class(fit, "mcglm")
  expect_true(fit$is_multinomial)
  expect_equal(fit$J, d$J)
  expect_equal(fit$K, d$K)
  # (J-1) blocks of (s + r) = 1 + 2 = 3 coefficients
  expect_length(fit$coefficients$naive, (d$J - 1) * (1 + ncol(d$x)))
  expect_true(all(grepl("^y[12]:", names(fit$coefficients$naive))))
})

test_that("multinomial naive supports K=3 (multicategory proxy)", {
  skip_if_not_installed("nnet")
  d <- .simulate_multinom(700, J = 3, K = 3, seed = 72)

  fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "multinomial",
               method = "naive", K = d$K, J = d$J)
  # (J-1) * (s + r) = 2 * (2 + 2) = 8 coefficients
  expect_length(fit$coefficients$naive, (d$J - 1) * ((d$K - 1) + ncol(d$x)))
})

# ==========================================================================
# Onestep multinomial path (RTMB)
# ==========================================================================

test_that("multinomial onestep runs and produces vcov (K=2, fix_omega=TRUE)", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("nnet")
  d <- .simulate_multinom(500, J = 3, K = 2, seed = 73)

  fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "multinomial",
               method = c("naive", "onestep"), K = d$K, J = d$J,
               p01 = d$Pi[2, 1], p10 = d$Pi[1, 2], pi_z = d$pi_z[2],
               fix_omega = TRUE)
  expect_named(fit$coefficients, c("naive", "onestep"))
  expect_true(all(is.finite(fit$coefficients$onestep)))
  V <- vcov(fit, method = "onestep")
  expect_true(is.matrix(V) && all(is.finite(V)))
})

test_that("multinomial onestep with default fix_omega=FALSE estimates weights", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("nnet")
  d <- .simulate_multinom(500, J = 3, K = 2, seed = 74)

  fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "multinomial",
               method = "onestep", K = d$K, J = d$J)
  expect_true(all(is.finite(fit$coefficients$onestep)))
})

# ==========================================================================
# Argument validation and unsupported S3 methods
# ==========================================================================

test_that("multinomial rejects bca/bcm/cs methods", {
  skip_if_not_installed("nnet")
  d <- .simulate_multinom(100, J = 3, K = 2, seed = 75)
  expect_error(
    mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "multinomial",
          method = "bca", K = d$K, J = d$J),
    "only 'naive' and 'onestep'"
  )
})

test_that("multinomial rejects out-of-range y", {
  skip_if_not_installed("nnet")
  d <- .simulate_multinom(100, J = 3, K = 2, seed = 76)
  bad_y <- d$y; bad_y[1] <- 99L
  expect_error(
    mcglm(bad_y, z_hat = d$z_hat, x = d$x, family = "multinomial",
          method = "naive", K = d$K, J = d$J),
    "must be integers"
  )
})

test_that("multinomial fits raise on fitted/predict/residuals", {
  skip_if_not_installed("nnet")
  d <- .simulate_multinom(200, J = 3, K = 2, seed = 77)
  fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "multinomial",
               method = "naive", K = d$K, J = d$J)
  expect_error(fitted(fit),     "not yet implemented")
  expect_error(predict(fit),    "not yet implemented")
  expect_error(residuals(fit),  "not yet implemented")
})
