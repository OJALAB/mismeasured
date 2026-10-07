# Corrected score of Akazawa, Kinukawa and Nakamura / Yi et al. (2019, eq.
# 17) ("cs_akn") with an estimated misclassification matrix.

.akv_data <- function(n = 1500L, K = 2L, family = "poisson", seed = 601L,
                      nv = 300L, type = "external", slope_z = 0) {
  set.seed(seed)
  if (K == 2L) {
    Pi <- matrix(c(0.86, 0.14, 0.12, 0.88), 2L)
    gamma <- c(0, 0.8)
  } else {
    Pi <- matrix(c(0.85, 0.10, 0.05, 0.10, 0.80, 0.10, 0.05, 0.10, 0.85), 3L)
    gamma <- c(0, 0.7, -0.6)
  }
  draw_z <- function(m, x1) {
    if (K == 2L) return(rbinom(m, 1L, plogis(-0.4 + slope_z * x1)))
    sample.int(K, m, TRUE, c(0.5, 0.3, 0.2)) - 1L
  }
  proxy <- function(z) vapply(z, function(k)
    sample.int(K, 1L, prob = Pi[, k + 1L]) - 1L, 0L)
  x1 <- rnorm(n)
  z <- draw_z(n, x1)
  z_hat <- proxy(z)
  eta <- gamma[z + 1L] - 0.3 + 0.5 * x1
  y <- switch(family, poisson = rpois(n, exp(eta)),
              binomial = rbinom(n, 1L, plogis(eta)))
  v <- if (type == "internal") {
    idx <- sample.int(n, nv)
    list(z = z[idx], index = idx)
  } else {
    zv <- draw_z(nv, rnorm(nv))
    list(z = zv, z_hat = proxy(zv))
  }
  list(y = y, z = z, z_hat = z_hat, x = cbind(1, x1), K = K, Pi = Pi,
       family = family, v = v, psi = c(gamma[-1L], -0.3, 0.5))
}

test_that("the AKN surrogate is the column of Pi^{-1} at the proxy", {
  Pi <- matrix(c(0.85, 0.10, 0.05, 0.10, 0.80, 0.10, 0.05, 0.10, 0.85), 3L)
  z_hat <- c(0L, 1L, 2L, 2L, 0L)
  Q <- mismeasured:::.akn_build_Q(Pi, 3L)
  x_akn <- mismeasured:::.akn_unbiased_surrogate(z_hat, Q$Q_inv, Q$p0, 3L)
  expect_equal(x_akn, t(solve(Pi))[z_hat + 1L, -1L], tolerance = 1e-12)
})

test_that("validated cs_akn matches an independent stacked reference", {
  for (cfg in list(c("poisson", 2L, "external"), c("binomial", 2L, "internal"),
                   c("poisson", 3L, "internal"), c("binomial", 3L, "external"))) {
    d <- .akv_data(n = 900L, K = as.integer(cfg[2]), family = cfg[1],
                   type = cfg[3], nv = 250L, seed = 610L)
    fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = d$family,
                 method = "cs_akn", validation = d$v)
    K <- d$K; s <- K - 1L; n <- length(d$y)
    psi <- unname(coef(fit, method = "cs_akn"))
    p <- length(psi)
    internal <- !is.null(d$v$index)
    zv <- d$v$z
    zhv <- if (internal) d$z_hat[d$v$index] else d$v$z_hat
    Pi_hat <- unname(fit$nuisance$cs_akn$Pi)
    linkinv <- if (d$family == "poisson") exp else plogis
    rows <- function(full) {
      b <- full[seq_len(p)]
      Pi <- rbind(1 - colSums(matrix(full[p + seq_len(K * s)], s, K)),
                  matrix(full[p + seq_len(K * s)], s, K))
      wts <- t(solve(Pi))[d$z_hat + 1L, , drop = FALSE]  # n x K, rows sum to 1
      eta0 <- drop(d$x %*% b[s + seq_len(2L)])
      U <- matrix(0, n, p)
      for (l in seq_len(K)) {
        xi <- cbind(matrix(0, n, s), d$x)
        if (l > 1L) xi[, l - 1L] <- 1
        U <- U + wts[, l] * xi * (d$y - linkinv(eta0 + c(0, b[seq_len(s)])[l]))
      }
      nuis <- do.call(cbind, lapply(0:s, function(l)
        sapply(seq_len(s), function(j) (zv == l) * ((zhv == j) - Pi[j + 1L, l + 1L]))))
      if (internal) {
        idx <- d$v$index
        xi <- cbind(outer(zv, seq_len(s), "==") * 1, d$x[idx, , drop = FALSE])
        U[idx, ] <- xi * (d$y[idx] - linkinv(drop(xi %*% b)))
        G <- cbind(U, matrix(0, n, K * s))
        G[idx, p + seq_len(K * s)] <- nuis
        G
      } else {
        rbind(cbind(U, matrix(0, n, K * s)),
              cbind(matrix(0, length(zv), p), nuis))
      }
    }
    full <- c(psi, as.numeric(Pi_hat[-1L, , drop = FALSE]))
    G <- rows(full)
    expect_lt(max(abs(colSums(G))), 1e-6 * n)
    A <- vapply(seq_along(full), function(k) {
      h <- 1e-6; e <- replace(numeric(length(full)), k, h)
      (colSums(rows(full + e)) - colSums(rows(full - e))) / (2 * h)
    }, numeric(length(full)))
    V <- solve(A) %*% crossprod(G) %*% t(solve(A))
    expect_equal(unname(vcov(fit, method = "cs_akn")), V[seq_len(p), seq_len(p)],
                 tolerance = 1e-5, label = paste(cfg, collapse = " "))
  }
})

test_that("cs_akn stays consistent when Z depends on x; sub does not", {
  skip_on_cran()
  d <- .akv_data(n = 40000L, nv = 3000L, slope_z = 1.2, seed = 620L)
  fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, method = c("sub", "cs_akn"),
               validation = d$v, mc_control = control_mc(prevalence = "em"))
  expect_lt(max(abs(coef(fit, method = "cs_akn") - d$psi)), 0.07)
  expect_gt(abs(coef(fit, method = "sub")[1] - d$psi[1]), 0.08)
})

test_that("ill-conditioned Q is flagged for known and estimated Pi", {
  d <- .akv_data(n = 300L, seed = 630L)
  weak <- matrix(c(0.51, 0.49, 0.48, 0.52), 2L)
  expect_warning(mcglm(d$y, z_hat = d$z_hat, x = d$x, method = "cs_akn",
                       Pi = weak), "cs_akn: \\|1 - p01 - p10\\|")
  expect_error(mcglm(d$y, z_hat = d$z_hat, x = d$x, method = "cs_akn",
                     Pi = weak, mc_control = control_mc(on_ill = "error")),
               "unstable")
  # columns 2 and 3 nearly equal: kappa(Q) is about 104
  Pi3 <- matrix(c(0.8, 0.1, 0.1, 0.3, 0.4, 0.3, 0.3, 0.395, 0.305), 3L)
  expect_false(suppressWarnings(
    mismeasured:::.akn_check_Q(Pi3, 3L, control_mc())))
  expect_true(mismeasured:::.akn_check_Q(diag(3), 3L, control_mc()))
})

test_that("cs_akn validation follows the shared rules", {
  d <- .akv_data(n = 400L, seed = 640L, type = "internal", nv = 150L)
  expect_error(mcglm(d$y, z_hat = d$z_hat, x = d$x, method = "cs_akn",
                     validation = d$v, Pi = d$Pi), "remove Pi")
  p_incl <- ifelse(d$z_hat == 1L, 0.4, 0.2)
  idx <- which(runif(400L) < p_incl)
  fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, method = c("sub", "cs_akn"),
               validation = validation_sample(d$z[idx], index = idx,
                                              weights = 1 / p_incl[idx]))
  expect_identical(fit$nuisance$cs_akn$beta_equation, "weighted")
  expect_true(all(is.finite(fit$se$cs_akn)))
  expect_named(fit$nuisance$cs_akn$eta, c("Pi[1,0]", "Pi[1,1]"))
})

test_that("validated cs_akn is unbiased and calibrated under x-dependent prevalence", {
  skip_on_cran()
  reps <- 100L
  for (type in c("external", "internal")) {
    out <- matrix(NA_real_, reps, 2L)
    for (r in seq_len(reps)) {
      d <- .akv_data(n = 1500L, nv = 300L, type = type, slope_z = 1,
                     seed = 7000L + r)
      fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, method = "cs_akn",
                   validation = d$v, mc_control = control_mc(on_ill = "none"))
      out[r, ] <- c(coef(fit, method = "cs_akn")[1], fit$se$cs_akn[1])
    }
    expect_lt(abs(mean(out[, 1]) - 0.8), 0.05, label = type)
    ratio <- mean(out[, 2]) / sd(out[, 1])
    expect_true(ratio > 0.8 && ratio < 1.3, label = paste(type, ratio))
  }
})
