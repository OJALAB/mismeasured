# Subtraction-correction estimator (method = "sub") with known
# misclassification probabilities.

.sub_data <- function(n, K = 2L, family = "poisson", seed = 101L) {
  set.seed(seed)
  if (K == 2L) {
    Pi <- matrix(c(0.88, 0.12, 0.15, 0.85), 2L)
    prev <- c(0.6, 0.4)
    gamma <- c(0, 0.8)
  } else {
    Pi <- matrix(c(0.85, 0.10, 0.05,
                   0.10, 0.80, 0.10,
                   0.05, 0.10, 0.85), 3L)
    prev <- c(0.5, 0.3, 0.2)
    gamma <- c(0, 0.7, -0.6)
  }
  z <- sample.int(K, n, TRUE, prev) - 1L
  z_hat <- vapply(z, function(k) sample.int(K, 1L, prob = Pi[, k + 1L]) - 1L,
                  0L)
  x1 <- rnorm(n)
  eta <- gamma[z + 1L] - 0.3 + 0.5 * x1
  y <- switch(family,
              poisson  = rpois(n, exp(eta)),
              binomial = rbinom(n, 1L, plogis(eta)),
              gaussian = rnorm(n, eta, 0.7))
  list(y = y, z = z, z_hat = z_hat, x = cbind(1, x1), x1 = x1, Pi = Pi,
       prev = prev, psi = c(gamma[-1L], -0.3, 0.5),
       pi_z = if (K == 2L) prev[2L] else prev)
}

test_that("gaussian sub equals its closed-form instrumental-variable solution", {
  for (K in 2:3) {
    d <- .sub_data(800L, K = K, family = "gaussian")
    fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "gaussian",
                 method = "sub", Pi = d$Pi, pi_z = d$pi_z)
    joint <- sweep(d$Pi, 2, d$prev, "*")
    W <- joint / rowSums(joint)
    xi_hat <- fit$xi_hat
    D <- cbind(W[d$z_hat + 1L, -1L, drop = FALSE], d$x)
    closed <- solve(crossprod(xi_hat, D), crossprod(xi_hat, d$y))
    expect_equal(unname(coef(fit, method = "sub")), as.numeric(closed),
                 tolerance = 1e-8)
    expect_true(fit$convergence$sub$converged)
  }
})

test_that("sub recovers the true parameters where naive is biased", {
  skip_on_cran()
  for (fam in c("poisson", "binomial")) {
    d <- .sub_data(30000L, family = fam, seed = 7L)
    fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = fam,
                 method = c("naive", "sub"), Pi = d$Pi, pi_z = d$pi_z)
    expect_lt(max(abs(coef(fit, method = "sub") - d$psi)), 0.08)
    expect_gt(abs(coef(fit, method = "naive")[1] - d$psi[1]), 0.15)
  }
  d <- .sub_data(30000L, K = 3L, family = "poisson", seed = 8L)
  fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "poisson",
               method = "sub", Pi = d$Pi, pi_z = d$pi_z)
  expect_lt(max(abs(coef(fit, method = "sub") - d$psi)), 0.08)
})

test_that("sub is no less efficient than cs (Rao-Blackwell ordering)", {
  skip_on_cran()
  set.seed(11)
  n <- 40000L
  Pi <- matrix(c(0.80, 0.20, 0.25, 0.75), 2L)
  z <- rbinom(n, 1L, 0.4)
  z_hat <- ifelse(z == 1L, rbinom(n, 1L, 0.75), rbinom(n, 1L, 0.20))
  x1 <- rnorm(n)
  y <- rpois(n, exp(-0.3 + 1.5 * z + 0.5 * x1))
  fit <- mcglm(y, z_hat = z_hat, x = cbind(1, x1), family = "poisson",
               method = c("cs", "sub"), Pi = Pi, pi_z = 0.4)
  ratio <- diag(vcov(fit, method = "sub")) / diag(vcov(fit, method = "cs"))
  # Population ordering is exact; sample sandwiches scatter around it.
  expect_true(all(ratio < 1.01))
  # The noise term b = -delta(x)(z_hat - pbar)v(x) scales with x, so the
  # gain is clearest on the slope.
  expect_lt(ratio[3], 0.95)
})

test_that("sub sandwich matches a numerical-Jacobian reference", {
  d <- .sub_data(600L, K = 3L, family = "binomial", seed = 21L)
  fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "binomial",
               method = "sub", Pi = d$Pi, pi_z = d$pi_z)
  psi <- unname(coef(fit, method = "sub"))
  joint <- sweep(d$Pi, 2, d$prev, "*")
  W <- joint / rowSums(joint)
  rows <- function(b) {
    eta <- drop(d$x %*% b[3:4])
    mu <- plogis(outer(eta, c(0, b[1:2]), "+"))
    fit$xi_hat * (d$y - rowSums(W[d$z_hat + 1L, ] * mu))
  }
  expect_lt(max(abs(colMeans(rows(psi)))), 1e-8)
  J <- vapply(seq_along(psi), function(j) {
    h <- 1e-6; e <- replace(numeric(4), j, h)
    (colMeans(rows(psi + e)) - colMeans(rows(psi - e))) / (2 * h)
  }, numeric(4))
  S <- crossprod(rows(psi)) / 600
  V <- solve(J) %*% S %*% t(solve(J)) / 600
  expect_equal(unname(vcov(fit, method = "sub")), V, tolerance = 1e-6)
})

test_that("integer frequency weights equal replicated rows for sub", {
  d <- .sub_data(300L, family = "poisson", seed = 31L)
  w <- rep(c(1L, 2L, 3L), length.out = 300L)
  idx <- rep(seq_len(300L), w)
  fw <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "poisson",
              method = "sub", Pi = d$Pi, pi_z = d$pi_z, weights = w)
  fr <- mcglm(d$y[idx], z_hat = d$z_hat[idx], x = d$x[idx, ],
              family = "poisson", method = "sub", Pi = d$Pi, pi_z = d$pi_z)
  expect_equal(coef(fw, method = "sub"), coef(fr, method = "sub"),
               tolerance = 1e-8)
  expect_equal(vcov(fw, method = "sub"), vcov(fr, method = "sub"),
               tolerance = 1e-8)
})

test_that("sub accepts p01/p10 and the formula interface", {
  d <- .sub_data(500L, family = "binomial", seed = 41L)
  df <- data.frame(y = d$y, z = factor(d$z_hat), x1 = d$x1)
  f_pi <- mcglm(y ~ mc(z, d$Pi) + x1, data = df, family = "binomial",
                method = "sub", pi_z = d$pi_z)
  f_pp <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "binomial",
                method = "sub", p01 = d$Pi[2, 1], p10 = d$Pi[1, 2],
                pi_z = d$pi_z)
  expect_equal(unname(coef(f_pi, method = "sub")),
               unname(coef(f_pp, method = "sub")), tolerance = 1e-10)
  expect_equal(names(coef(f_pi, method = "sub")),
               c("gamma", "(Intercept)", "x1"))
})

test_that("sub without P(Z | Z_hat) information errors informatively", {
  d <- .sub_data(200L, seed = 51L)
  expect_error(
    mcglm(d$y, z_hat = d$z_hat, x = d$x, method = "sub", c1 = 0.05,
          c2 = -0.01),
    "For sub, supply Pi")
  expect_error(
    mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "multinomial",
          method = "sub"),
    "Unsupported: sub")
})
