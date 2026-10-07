# Expectation correction ("ec") and induced likelihood ("il"), with known
# probabilities and with internal / external validation samples. References
# are written independently with dpois/dbinom/dnorm and numerical
# derivatives.

.ei_data <- function(n = 900L, K = 2L, family = "poisson", seed = 401L,
                     nv = 250L, type = "external") {
  set.seed(seed)
  if (K == 2L) {
    Pi <- matrix(c(0.86, 0.14, 0.12, 0.88), 2L)
    prev <- c(0.6, 0.4)
    gamma <- c(0, 0.8)
  } else {
    Pi <- matrix(c(0.85, 0.10, 0.05, 0.10, 0.80, 0.10, 0.05, 0.10, 0.85), 3L)
    prev <- c(0.5, 0.3, 0.2)
    gamma <- c(0, 0.7, -0.6)
  }
  lab <- function(m) {
    z <- sample.int(K, m, TRUE, prev) - 1L
    zh <- vapply(z, function(k) sample.int(K, 1L, prob = Pi[, k + 1L]) - 1L,
                 0L)
    list(z = z, z_hat = zh)
  }
  d <- lab(n)
  x <- cbind(1, rnorm(n))
  eta <- gamma[d$z + 1L] - 0.3 + 0.5 * x[, 2L]
  y <- switch(family, poisson = rpois(n, exp(eta)),
              binomial = rbinom(n, 1L, plogis(eta)),
              gaussian = rnorm(n, eta, 0.7))
  v <- if (type == "internal") {
    idx <- sample.int(n, nv)
    list(z = d$z[idx], index = idx)
  } else lab(nv)
  list(y = y, z = d$z, z_hat = d$z_hat, x = x, K = K, family = family,
       Pi = Pi, prev = prev, v = v, psi = c(gamma[-1L], -0.3, 0.5))
}

# n x K class log densities
.ei_logf <- function(theta, d) {
  s <- d$K - 1L
  p <- s + ncol(d$x)
  eta <- outer(drop(d$x %*% theta[s + seq_len(ncol(d$x))]),
               c(0, theta[seq_len(s)]), "+")
  n <- length(d$y)
  matrix(switch(d$family,
    poisson  = dpois(d$y, exp(eta), log = TRUE),
    binomial = dbinom(d$y, 1L, plogis(eta), log = TRUE),
    gaussian = dnorm(d$y, eta, exp(theta[p + 1L]), log = TRUE)), n, d$K)
}
.ei_grad_rows <- function(f, theta, h = 1e-5) {
  vapply(seq_along(theta), function(k) {
    e <- replace(numeric(length(theta)), k, h)
    (f(theta + e) - f(theta - e)) / (2 * h)
  }, numeric(length(f(theta))))
}
.ei_jac <- function(f, theta, h = 1e-4) {
  vapply(seq_along(theta), function(k) {
    e <- replace(numeric(length(theta)), k, h)
    (f(theta + e) - f(theta - e)) / (2 * h)
  }, numeric(length(f(theta))))
}
.ei_theta <- function(fit, m, d) {
  th <- unname(coef(fit, method = m))
  if (d$family == "gaussian") th <- c(th, log(fit$nuisance[[m]]$sigma))
  th
}

test_that("with known probabilities ec and il coincide and match onestep", {
  for (fam in c("poisson", "binomial", "gaussian")) {
    d <- .ei_data(family = fam, seed = 410L)
    fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = fam,
                 method = c("ec", "il"), Pi = d$Pi, pi_z = d$prev[2])
    expect_equal(coef(fit, method = "ec"), coef(fit, method = "il"))
    expect_true(fit$convergence$ec$converged)
    skip_if_not_installed("RTMB")
    os <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = fam,
                method = "onestep", fix_omega = TRUE, p01 = d$Pi[2, 1],
                p10 = d$Pi[1, 2], pi_z = d$prev[2])
    expect_equal(unname(coef(fit, method = "ec")),
                 unname(coef(os, method = "onestep")), tolerance = 1e-4)
  }
})

test_that("known-probability sandwich matches a numerical reference", {
  for (cfg in list(c("gaussian", 2L), c("poisson", 3L))) {
    d <- .ei_data(family = cfg[1], K = as.integer(cfg[2]), seed = 420L)
    fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = d$family,
                 method = "ec", Pi = d$Pi,
                 pi_z = if (d$K == 2L) d$prev[2] else d$prev)
    th <- .ei_theta(fit, "ec", d)
    joint <- t(t(d$Pi) * d$prev)
    W <- joint / rowSums(joint)
    ll <- function(t) log(rowSums(exp(.ei_logf(t, d)) * W[d$z_hat + 1L, ]))
    U <- .ei_grad_rows(ll, th)
    expect_lt(max(abs(colSums(U))), 1e-5 * length(d$y))
    A <- .ei_jac(function(t) colSums(.ei_grad_rows(ll, t)), th)
    V <- solve(A) %*% crossprod(U) %*% t(solve(A))
    p <- length(coef(fit, method = "ec"))
    expect_equal(unname(vcov(fit, method = "ec")), V[seq_len(p), seq_len(p)],
                 tolerance = 1e-3)
  }
})

test_that("validated ec matches an independent stacked reference", {
  for (cfg in list(c("binomial", 2L, "external"), c("gaussian", 2L, "internal"),
                   c("poisson", 3L, "internal"))) {
    d <- .ei_data(family = cfg[1], K = as.integer(cfg[2]), type = cfg[3],
                  n = 700L, nv = 220L, seed = 430L)
    fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = d$family,
                 method = "ec", validation = d$v)
    K <- d$K; s <- K - 1L; n <- length(d$y)
    th <- .ei_theta(fit, "ec", d)
    eta <- unname(fit$mc_estimate$eta)
    q <- length(th)
    internal <- !is.null(d$v$index)
    zv <- d$v$z
    zhv <- if (internal) d$z_hat[d$v$index] else d$v$z_hat
    rows <- function(full) {
      t <- full[seq_len(q)]; e <- full[q + seq_len(length(eta))]
      Pi <- rbind(1 - colSums(matrix(e[seq_len(K * s)], s, K)),
                  matrix(e[seq_len(K * s)], s, K))
      pi <- c(1 - sum(e[K * s + seq_len(s)]), e[K * s + seq_len(s)])
      joint <- t(t(Pi) * pi)
      W <- joint / rowSums(joint)
      lmix <- function(tt) log(rowSums(exp(.ei_logf(tt, d)) *
                                         W[d$z_hat + 1L, , drop = FALSE]))
      U <- .ei_grad_rows(lmix, t)
      nuis <- cbind(do.call(cbind, lapply(0:s, function(l)
        sapply(seq_len(s), function(j) (zv == l) * ((zhv == j) - Pi[j + 1L, l + 1L])))),
        sapply(seq_len(s), function(l) (zv == l) - pi[l + 1L]))
      if (internal) {
        idx <- d$v$index
        lc <- function(tt) .ei_logf(tt, d)[cbind(idx, zv + 1L)]
        U[idx, ] <- .ei_grad_rows(lc, t)
        G <- cbind(U, matrix(0, n, length(eta)))
        G[idx, q + seq_along(eta)] <- nuis
        G
      } else {
        rbind(cbind(U, matrix(0, n, length(eta))),
              cbind(matrix(0, length(zv), q), nuis))
      }
    }
    full <- c(th, eta)
    G <- rows(full)
    expect_lt(max(abs(colSums(G))), 1e-4 * n)
    A <- .ei_jac(function(f) colSums(rows(f)), full)
    V <- solve(A) %*% crossprod(G) %*% t(solve(A))
    p <- length(coef(fit, method = "ec"))
    expect_equal(unname(vcov(fit, method = "ec")), V[seq_len(p), seq_len(p)],
                 tolerance = 2e-3, label = paste(cfg, collapse = " "))
  }
})

test_that("validated il maximises the joint likelihood; sandwich reference", {
  for (cfg in list(c("poisson", 2L, "external", "validation"),
                   c("binomial", 2L, "external", "em"),
                   c("gaussian", 2L, "internal", "validation"),
                   c("poisson", 3L, "internal", "validation"))) {
    d <- .ei_data(family = cfg[1], K = as.integer(cfg[2]), type = cfg[3],
                  n = 700L, nv = 220L, seed = 440L)
    fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = d$family,
                 method = "il", validation = d$v,
                 mc_control = control_mc(prevalence = cfg[4]))
    K <- d$K; s <- K - 1L; n <- length(d$y)
    internal <- !is.null(d$v$index)
    th <- .ei_theta(fit, "il", d)
    q <- length(th)
    Pi <- unname(fit$nuisance$il$Pi); pi <- unname(fit$nuisance$il$pi_z)
    full <- c(th, as.numeric(log(Pi[-1L, , drop = FALSE]) -
                               rep(log(Pi[1L, ]), each = s)),
              log(pi[-1L] / pi[1L]))
    zv <- d$v$z
    zhv <- if (internal) d$z_hat[d$v$index] else d$v$z_hat
    unit_ll <- function(f) {
      t <- f[seq_len(q)]
      a <- rbind(0, matrix(f[q + seq_len(K * s)], s, K))
      P <- exp(a) / rep(colSums(exp(a)), each = K)
      b <- c(0, f[q + K * s + seq_len(s)])
      pr <- exp(b) / sum(exp(b))
      lf <- .ei_logf(t, d)
      ll <- log(rowSums(exp(lf) * P[d$z_hat + 1L, , drop = FALSE] *
                          rep(pr, each = n)))
      if (internal) {
        idx <- d$v$index
        ll[idx] <- lf[cbind(idx, zv + 1L)] +
          log(P[cbind(zhv + 1L, zv + 1L)]) + log(pr[zv + 1L])
        return(ll)
      }
      c(ll, log(P[cbind(zhv + 1L, zv + 1L)]) +
          if (cfg[4] == "validation") log(pr[zv + 1L]) else 0)
    }
    G <- .ei_grad_rows(unit_ll, full)
    expect_lt(max(abs(colSums(G))), 1e-4 * n)
    A <- .ei_jac(function(f) colSums(.ei_grad_rows(unit_ll, f)), full)
    V <- solve(A) %*% crossprod(G) %*% t(solve(A))
    p <- length(coef(fit, method = "il"))
    expect_equal(unname(vcov(fit, method = "il")), V[seq_len(p), seq_len(p)],
                 tolerance = 2e-3, label = paste(cfg, collapse = " "))
  }
})

test_that("ec/il are at least as efficient as sub with known probabilities", {
  skip_on_cran()
  d <- .ei_data(n = 30000L, family = "poisson", seed = 450L)
  fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, method = c("sub", "ec"),
               Pi = d$Pi, pi_z = d$prev[2])
  ratio <- diag(vcov(fit, method = "ec")) / diag(vcov(fit, method = "sub"))
  expect_true(all(ratio < 1.01))
  expect_lt(max(abs(coef(fit, method = "ec") - d$psi)), 0.06)
})

test_that("input checks for ec and il", {
  d <- .ei_data(n = 200L, family = "binomial", seed = 460L)
  y <- d$y
  y[1] <- 0.5
  expect_error(mcglm(y, z_hat = d$z_hat, x = d$x, family = "binomial",
                     method = "ec", Pi = d$Pi, pi_z = 0.4), "0/1 response")
  expect_error(mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "binomial",
                     method = "il", c1 = 0.05, c2 = 0.01), "For il, supply Pi")
  expect_error(
    mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "binomial", method = "il",
          validation = d$v, mc_control = control_mc(prevalence = "inverse")),
    "maximum likelihood")
})

test_that("validated ec and il are unbiased and calibrated", {
  skip_on_cran()
  reps <- 100L
  for (type in c("external", "internal")) {
    out <- matrix(NA_real_, reps, 4L)
    for (r in seq_len(reps)) {
      d <- .ei_data(n = 1500L, nv = 300L, type = type, seed = 5000L + r)
      v <- d$v
      if (type == "internal") {
        # outcome-dependent audit: EC uses the design weights, IL needs none
        p_incl <- ifelse(d$y == 0, 0.08, 0.25)
        idx <- which(runif(length(d$y)) < p_incl)
        v <- validation_sample(d$z[idx], index = idx, weights = 1 / p_incl[idx])
      }
      fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, method = c("ec", "il"),
                   validation = v, mc_control = control_mc(on_ill = "none"))
      out[r, ] <- c(coef(fit, method = "ec")[1], fit$se$ec[1],
                    coef(fit, method = "il")[1], fit$se$il[1])
    }
    for (k in c(1L, 3L)) {
      expect_lt(abs(mean(out[, k]) - 0.8), 0.04, label = type)
      ratio <- mean(out[, k + 1L]) / sd(out[, k])
      expect_true(ratio > 0.8 && ratio < 1.25, label = paste(type, k, ratio))
    }
  }
})
