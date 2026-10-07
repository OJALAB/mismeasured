# Propagating estimated misclassification probabilities into the variance of
# the plug-in estimators cs, bca, bcm and onestep.

.pp_data <- function(n = 1500L, nv = 300L, seed = 1101L, type = "external") {
  set.seed(seed)
  z <- rbinom(n, 1L, 0.4)
  zh <- ifelse(z == 1L, rbinom(n, 1L, 0.86), rbinom(n, 1L, 0.12))
  x <- cbind(1, rnorm(n))
  y <- rpois(n, exp(-0.3 + 0.8 * z + 0.5 * x[, 2L]))
  v <- if (type == "internal") {
    idx <- sample.int(n, nv)
    list(z = z[idx], index = idx)
  } else {
    zv <- rbinom(nv, 1L, 0.4)
    list(z = zv, z_hat = ifelse(zv == 1L, rbinom(nv, 1L, 0.86),
                                rbinom(nv, 1L, 0.12)))
  }
  list(y = y, z = z, z_hat = zh, x = x, v = v)
}

test_that("cs with a design-weighted audit matches an independent stacked reference", {
  d <- .pp_data(n = 900L, seed = 1110L)
  p_incl <- ifelse(d$z_hat == 1L, 0.35, 0.12)
  idx <- which(runif(900L) < p_incl)
  dw <- 1 / p_incl[idx]
  val <- validation_sample(d$z[idx], index = idx, weights = dw)
  fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, method = "cs", validation = val)
  psi <- unname(coef(fit, method = "cs"))
  est <- fit$mc_estimate
  zv <- d$z[idx]
  zhv <- d$z_hat[idx]
  n <- length(d$y)
  xi_hat <- cbind(d$z_hat, d$x)
  rows <- function(f) {
    b <- f[1:3]
    p01 <- f[4]; p11 <- f[5]; pi1 <- f[6]
    eta0 <- drop(d$x %*% b[2:3])
    delta <- exp(eta0 + b[1]) - exp(eta0)
    c1 <- p01 * (1 - pi1)
    c2 <- c1 - (1 - p11) * pi1
    phi <- xi_hat * (d$y - exp(drop(xi_hat %*% b))) +
      cbind(c1 * delta, c2 * delta * d$x)
    G <- cbind(phi, matrix(0, n, 3))
    G[idx, 4:6] <- dw * cbind((zv == 0) * (zhv - p01), (zv == 1) * (zhv - p11),
                              zv - pi1)
    G
  }
  full <- c(psi, unname(est$eta))
  G <- rows(full)
  expect_lt(max(abs(colSums(G))), 1e-6 * n)
  A <- vapply(1:6, function(k) {
    h <- 1e-6; e <- replace(numeric(6), k, h)
    (colSums(rows(full + e)) - colSums(rows(full - e))) / (2 * h)
  }, numeric(6))
  V <- solve(A) %*% crossprod(G) %*% t(solve(A))
  expect_equal(unname(vcov(fit, method = "cs")), V[1:3, 1:3], tolerance = 1e-6)
})

test_that("bca/bcm covariance stacks their estimating equations with eta", {
  for (type in c("bca", "bcm")) {
    d <- .pp_data(n = 900L, seed = 1120L)
    fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, method = c("naive", type),
                 validation = d$v)
    est <- fit$mc_estimate
    n <- length(d$y)
    xi_hat <- cbind(d$z_hat, d$x)
    zv <- d$v$z; zhv <- d$v$z_hat
    drift <- function(b, f) {          # binary drift -(c1, c2 x) delta
      c1 <- f[1] * (1 - f[3])
      c2 <- c1 - (1 - f[2]) * f[3]
      eta0 <- drop(d$x %*% b[2:3])
      delta <- exp(eta0 + b[1]) - exp(eta0)
      cbind(-c1 * delta, -c2 * delta * d$x)
    }
    rows <- function(th) {
      bn <- th[1:3]; b <- th[4:6]; f <- th[7:9]
      mu_n <- exp(drop(xi_hat %*% bn))
      S <- xi_hat * (d$y - mu_n)
      Id <- xi_hat * (mu_n * drop(xi_hat %*% (b - bn)))
      h <- if (type == "bca") Id + drift(bn, f) else {
        e <- 1e-6
        Md <- (drift(bn + e * (b - bn), f) - drift(bn - e * (b - bn), f)) / (2 * e)
        Id + Md - (S - drift(bn, f))
      }
      Gv <- cbind((zv == 0) * (zhv - f[1]), (zv == 1) * (zhv - f[2]), zv - f[3])
      rbind(cbind(S, h, matrix(0, n, 3)), cbind(matrix(0, length(zv), 6), Gv))
    }
    th <- c(unname(coef(fit, method = "naive")), unname(coef(fit, method = type)),
            unname(est$eta))
    G <- rows(th)
    expect_lt(max(abs(colSums(G))), 1e-6 * n)
    A <- vapply(1:9, function(k) {
      h <- 1e-6; e <- replace(numeric(9), k, h)
      (colSums(rows(th + e)) - colSums(rows(th - e))) / (2 * h)
    }, numeric(9))
    V <- solve(A) %*% crossprod(G) %*% t(solve(A))
    expect_equal(unname(vcov(fit, method = type)), V[4:6, 4:6],
                 tolerance = 1e-5, label = type)
  }
})

test_that("onestep with estimated probabilities matches ec", {
  skip_if_not_installed("RTMB")
  d <- .pp_data(n = 900L, seed = 1130L)
  fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, method = c("ec", "onestep"),
               fix_omega = TRUE, validation = d$v)
  expect_equal(unname(coef(fit, method = "onestep")),
               unname(coef(fit, method = "ec")), tolerance = 1e-4)
  expect_equal(unname(vcov(fit, method = "onestep")),
               unname(vcov(fit, method = "ec")), tolerance = 1e-3)
  expect_error(mcglm(d$y, z_hat = d$z_hat, x = d$x, method = "onestep",
                     validation = d$v), "fix_omega = FALSE")
})

test_that("supplied probabilities must equal the audit's estimates", {
  d <- .pp_data(n = 400L, seed = 1140L)
  est <- estimate_mc(d$v)
  ok <- mcglm(d$y, z_hat = d$z_hat, x = d$x, method = "bca",
              validation = d$v, Pi = unname(est$Pi), pi_z = est$pi[[2]])
  expect_equal(coef(ok, method = "bca"),
               coef(mcglm(d$y, z_hat = d$z_hat, x = d$x, method = "bca",
                          validation = d$v), method = "bca"))
  expect_error(mcglm(d$y, z_hat = d$z_hat, x = d$x, method = "bca",
                     validation = d$v, Pi = unname(est$Pi), pi_z = 0.3),
               "empirical estimates")
  expect_error(mcglm(d$y, z_hat = d$z_hat, x = d$x, method = "naive",
                     validation = d$v), "requires a corrected method")
})

test_that("propagated cs/bca/bcm standard errors are calibrated", {
  skip_on_cran()
  reps <- 150L
  for (type in c("external", "internal")) {
    out <- matrix(NA_real_, reps, 6L)
    for (r in seq_len(reps)) {
      d <- .pp_data(type = type, seed = 12000L + r)
      fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x,
                   method = c("bca", "bcm", "cs"), validation = d$v,
                   mc_control = control_mc(on_ill = "none"))
      out[r, ] <- c(sapply(c("bca", "bcm", "cs"), function(m)
        c(coef(fit, method = m)[1], fit$se[[m]][1])))
    }
    for (k in c(1L, 3L, 5L)) {
      ratio <- mean(out[, k + 1L]) / sd(out[, k])
      expect_true(ratio > 0.8 && ratio < 1.25, label = paste(type, k, ratio))
    }
    expect_lt(abs(mean(out[, 5]) - 0.8), 0.03, label = type)
  }
})
