# Subtraction correction with misclassification probabilities estimated
# from an internal or external validation sample (Yi et al., 2019, Sec. 4).

.subv_data <- function(K = 2L, type = "external", family = "poisson",
                       n = 1500L, nv = 300L, seed = 202L) {
  set.seed(seed)
  if (K == 2L) {
    Pi <- matrix(c(0.85, 0.15, 0.10, 0.90), 2L)
    prev <- c(0.6, 0.4)
    gamma <- c(0, 0.8)
  } else {
    Pi <- matrix(c(0.85, 0.10, 0.05,
                   0.10, 0.80, 0.10,
                   0.05, 0.10, 0.85), 3L)
    prev <- c(0.5, 0.3, 0.2)
    gamma <- c(0, 0.7, -0.6)
  }
  labels <- function(m) {
    z <- sample.int(K, m, TRUE, prev) - 1L
    zh <- vapply(z, function(k) sample.int(K, 1L, prob = Pi[, k + 1L]) - 1L,
                 0L)
    list(z = z, z_hat = zh)
  }
  d <- labels(n)
  x <- cbind(1, rnorm(n))
  eta <- gamma[d$z + 1L] - 0.3 + 0.5 * x[, 2L]
  y <- switch(family,
              poisson  = rpois(n, exp(eta)),
              binomial = rbinom(n, 1L, plogis(eta)),
              gaussian = rnorm(n, eta, 0.7))
  if (type == "internal") {
    index <- sample.int(n, nv)
    v <- list(z = d$z[index], index = index)
    proxy_v <- d$z_hat[index]
  } else {
    v <- labels(nv)
    proxy_v <- v$z_hat
  }
  list(y = y, z = d$z, z_hat = d$z_hat, x = x, v = v, proxy_v = proxy_v,
       K = K, family = family, psi = c(gamma[-1L], -0.3, 0.5))
}

.subv_fit <- function(d, prevalence = "validation", on_ill = "warn", ...) {
  mcglm(d$y, z_hat = d$z_hat, x = d$x, family = d$family,
        method = c("naive", "sub"), validation = d$v,
        mc_control = control_mc(prevalence = prevalence, on_ill = on_ill), ...)
}

# Independent stacked estimating function: rows for every unit, columns
# (psi, eta), written without the package's helpers.
.subv_reference_vcov <- function(d, fit, prevalence) {
  K <- d$K; s <- K - 1L; n <- length(d$y)
  fam <- stats::family(fit)
  internal <- !is.null(d$v$index)
  psi <- unname(coef(fit, method = "sub"))
  eta <- unname(fit$nuisance$sub$eta)
  p <- length(psi); q <- length(eta)
  xi_hat <- cbind(outer(d$z_hat, seq_len(s), "==") * 1, d$x)

  rows <- function(theta) {
    b <- theta[seq_len(p)]; e <- theta[p + seq_len(q)]
    Pi <- matrix(e[seq_len(K * s)], s, K)
    Pi <- rbind(1 - colSums(Pi), Pi)
    pr <- e[K * s + seq_len(s)]
    pi <- c(1 - sum(pr), pr)
    joint <- t(t(Pi) * pi)
    W <- joint / rowSums(joint)
    mu <- fam$linkinv(outer(drop(d$x %*% b[-seq_len(s)]), c(0, b[seq_len(s)]),
                            "+"))
    U <- xi_hat * (d$y - rowSums(W[d$z_hat + 1L, , drop = FALSE] * mu))
    nuis <- function(z, zh) {
      out <- NULL
      for (l in 0:s) for (j in seq_len(s))
        out <- cbind(out, (z == l) * ((zh == j) - Pi[j + 1L, l + 1L]))
      out
    }
    if (internal) {
      idx <- d$v$index
      zt <- d$v$z
      xi_t <- cbind(outer(zt, seq_len(s), "==") * 1, d$x[idx, , drop = FALSE])
      U[idx, ] <- xi_t * (d$y[idx] - mu[cbind(idx, zt + 1L)])
      G <- cbind(U, matrix(0, n, q))
      G[idx, p + seq_len(K * s)] <- nuis(zt, d$z_hat[idx])
      if (prevalence == "validation")
        G[idx, p + K * s + seq_len(s)] <-
          outer(zt, seq_len(s), "==") - rep(pr, each = length(idx))
    } else {
      zv <- d$v$z; m <- length(zv)
      Gv <- cbind(matrix(0, m, p), nuis(zv, d$v$z_hat), matrix(0, m, s))
      if (prevalence == "validation")
        Gv[, p + K * s + seq_len(s)] <-
          outer(zv, seq_len(s), "==") - rep(pr, each = m)
      G <- rbind(cbind(U, matrix(0, n, q)), Gv)
    }
    # EM fixed point: P(Z = l | Z_hat_i) - pi_l on every main unit
    if (prevalence == "em")
      G[seq_len(n), p + K * s + seq_len(s)] <-
        W[d$z_hat + 1L, -1L, drop = FALSE] - rep(pr, each = n)
    G
  }
  theta <- c(psi, eta)
  A <- vapply(seq_along(theta), function(k) {
    h <- 1e-6
    e <- replace(numeric(length(theta)), k, h)
    (colSums(rows(theta + e)) - colSums(rows(theta - e))) / (2 * h)
  }, numeric(length(theta)))
  G <- rows(theta)
  list(score = colSums(G),
       vcov = solve(A) %*% crossprod(G) %*% t(solve(A)))
}

test_that("nuisance estimates are the validation-sample proportions", {
  for (type in c("external", "internal")) {
    d <- .subv_data(K = 3L, type = type, seed = 210L)
    fit <- .subv_fit(d)
    tab <- table(factor(d$proxy_v, 0:2), factor(d$v$z, 0:2))
    expect_equal(unname(fit$nuisance$sub$Pi), unname(unclass(prop.table(tab, 2))))
    expect_equal(unname(fit$nuisance$sub$pi_z),
                 as.numeric(prop.table(table(d$v$z))))
    expect_identical(fit$validation_design, type)
  }
  d <- .subv_data(seed = 211L)
  fit <- .subv_fit(d, prevalence = "em")
  expect_equal(unname(fit$nuisance$sub$pi_z),
               as.numeric(solve(fit$nuisance$sub$Pi,
                                prop.table(table(d$z_hat)))))
})

test_that("stacked sandwich matches an independent numerical reference", {
  for (K in 2:3) for (type in c("external", "internal"))
    for (src in c("validation", "em")) {
      d <- .subv_data(K = K, type = type, family = "binomial",
                      n = 900L, nv = 250L, seed = 220L + K)
      fit <- .subv_fit(d, prevalence = src)
      ref <- .subv_reference_vcov(d, fit, src)
      # The fitted (psi, eta) solve the stacked equations, including the
      # true-category score on internal validation rows.
      expect_lt(max(abs(ref$score)), 1e-7)
      p <- length(coef(fit, method = "sub"))
      expect_equal(unname(vcov(fit, method = "sub")),
                   ref$vcov[seq_len(p), seq_len(p)], tolerance = 1e-5,
                   label = paste(K, type, src))
      expect_equal(unname(fit$nuisance$sub$vcov),
                   ref$vcov[-seq_len(p), -seq_len(p)], tolerance = 1e-5)
    }
})

test_that("integer frequency weights equal replicated main-study rows", {
  d <- .subv_data(n = 400L, nv = 150L, seed = 230L)
  w <- rep(1:3, length.out = 400L)
  idx <- rep(seq_len(400L), w)
  fw <- .subv_fit(d, weights = w)
  fr <- mcglm(d$y[idx], z_hat = d$z_hat[idx], x = d$x[idx, ],
              family = "poisson", method = "sub", validation = d$v)
  expect_equal(coef(fw, method = "sub"), coef(fr, method = "sub"),
               tolerance = 1e-8)
  expect_equal(vcov(fw, method = "sub"), vcov(fr, method = "sub"),
               tolerance = 1e-7)
})

test_that("supplied probabilities conflict with a validation sample for sub", {
  d <- .subv_data(n = 300L, nv = 100L, seed = 240L)
  Pi <- matrix(c(0.85, 0.15, 0.10, 0.90), 2L)
  expect_error(.subv_fit(d, Pi = Pi), "remove Pi")
  expect_error(.subv_fit(d, pi_z = 0.4), "remove pi_z")
  df <- data.frame(y = d$y, z = d$z_hat, x1 = d$x[, 2L])
  expect_error(mcglm(y ~ mc(z, Pi) + x1, data = df, method = "sub",
                     validation = d$v), "mc\\(\\)")
  ok <- mcglm(y ~ mc(z) + x1, data = df, method = "sub", validation = d$v)
  expect_equal(unname(coef(ok, method = "sub")),
               unname(coef(.subv_fit(d), method = "sub")))
})

test_that("cs and bca/bcm take the validation-sample probabilities", {
  d <- .subv_data(n = 500L, nv = 200L, type = "internal", seed = 250L)
  tab <- prop.table(table(factor(d$proxy_v, 0:1), factor(d$v$z, 0:1)))
  pi_v <- colSums(tab)
  supplied <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "poisson",
                    method = c("bca", "cs"), validation = d$v,
                    Pi = sweep(unclass(tab), 2, pi_v, "/"), pi_z = pi_v[2])
  derived <- mcglm(d$y, z_hat = d$z_hat, x = d$x, family = "poisson",
                   method = c("bca", "cs", "sub"), validation = d$v)
  for (m in c("bca", "cs")) {
    expect_equal(coef(derived, method = m), coef(supplied, method = m))
    expect_equal(vcov(derived, method = m), vcov(supplied, method = m))
  }
})

test_that("validation is rejected for methods that cannot use it", {
  d <- .subv_data(n = 300L, nv = 100L, seed = 260L)
  expect_error(mcglm(d$y, z_hat = d$z_hat, x = d$x, method = "bcm",
                     validation = d$v), "requires method")
  expect_error(mcglm(d$y, z_hat = d$z_hat, x = d$x,
                     method = c("sub", "onestep"), validation = d$v),
               "not supported for method\\(s\\) onestep")
  bad <- d
  bad$v$z <- bad$v$z * 0L
  expect_error(.subv_fit(bad), "Every true category")
})

test_that("main-study prevalence: EM boundary warns, inversion errors", {
  d <- .subv_data(n = 400L, nv = 200L, seed = 270L)
  # Validation proxies almost never flip, main proxies are rare: Pi^{-1} p
  # then leaves (0, 1).
  d$v$z_hat <- d$v$z
  d$v$z_hat[1:2] <- 1L - d$v$z_hat[1:2]
  d$z_hat[] <- 0L
  d$z_hat[1:3] <- 1L
  expect_error(suppressWarnings(.subv_fit(d, prevalence = "inverse")),
               "outside \\(0, 1\\)")
  msgs <- character()
  fit <- withCallingHandlers(
    .subv_fit(d, prevalence = "em"),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("boundary", msgs)))
  expect_true(any(grepl("singular Jacobian", msgs)))
  expect_true(all(is.na(fit$se$sub)))
})

test_that("validated sub is consistent and its SEs are calibrated", {
  skip_on_cran()
  for (type in c("external", "internal")) {
    est <- se <- matrix(NA_real_, 150L, 3L)
    for (r in seq_len(150L)) {
      d <- .subv_data(type = type, n = 1500L, nv = 300L, seed = 1000L + r)
      fit <- .subv_fit(d, on_ill = "none")
      est[r, ] <- coef(fit, method = "sub")
      se[r, ] <- fit$se$sub
    }
    expect_lt(max(abs(colMeans(est) - d$psi)), 0.04, label = type)
    ratio <- colMeans(se) / apply(est, 2, sd)
    expect_true(all(ratio > 0.8 & ratio < 1.25), label = type)
  }
})
