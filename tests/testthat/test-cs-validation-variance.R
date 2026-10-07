# Independent equation-level references for the estimated-probability CS
# sandwich. In particular, these do not call the new variance helpers.
.csv_fixture <- function(K = 2L, type = "internal", family = "poisson",
                         n = 700L, nv = 300L, seed = 723L) {
  set.seed(seed)
  Pi <- if (K == 2L) matrix(c(.85, .15, .1, .9), 2L) else
    matrix(c(.85, .1, .05, .1, .85, .05, .05, .1, .85), 3L)
  prevalence <- if (K == 2L) c(.6, .4) else c(.5, .3, .2)
  labels <- function(n) {
    z <- sample.int(K, n, TRUE, prevalence) - 1L
    proxy <- vapply(z, function(k) sample.int(K, 1L, prob = Pi[, k + 1L]) - 1L, 0L)
    list(z = z, z_hat = proxy)
  }
  d <- labels(n)
  x <- cbind(1, rnorm(n))
  eta <- c(0, .65, -.4)[d$z + 1L] - .4 + .5 * x[, 2L]
  y <- switch(family, poisson = rpois(n, exp(eta)),
              binomial = rbinom(n, 1, plogis(eta)), gaussian = rnorm(n, eta))
  if (type == "internal") {
    index <- sample.int(n, nv)
    v <- list(z = d$z[index], index = index)
    proxy_v <- d$z_hat[index]
  } else {
    v <- labels(nv)
    proxy_v <- v$z_hat
  }
  B <- matrix(tabulate(proxy_v + 1L + K * v$z, nbins = K * K), K) / nv
  pi <- colSums(B)
  args <- list(formula = y, z_hat = d$z_hat, x = x, family = family,
               Pi = sweep(B, 2, pi, "/"), pi_z = if (K == 2L) pi[2L] else pi,
               method = c("naive", "bca", "bcm", "cs"))
  list(args = args, v = v, B = B, proxy_v = proxy_v)
}

.csv_jacobian <- function(f, a, h = 1e-5) {
  vapply(seq_along(a), function(j) {
    plus <- minus <- a
    plus[j] <- plus[j] + h
    minus[j] <- minus[j] - h
    (f(plus) - f(minus)) / (2 * h)
  }, numeric(length(f(a))))
}

.csv_reference <- function(fit, d) {
  psi <- unname(coef(fit, method = "cs"))
  K <- fit$K
  q <- fit$x
  n <- fit$n
  nv <- length(d$v$z)
  mu <- fit$family$linkinv
  # The multicategory equation directly in joint probabilities, also at K=2.
  score <- function(beta, b) {
    gamma <- c(0, beta[seq_len(K - 1L)])
    eta <- drop(q %*% beta[-seq_len(K - 1L)])
    mus <- vapply(gamma, function(g) mu(eta + g), numeric(n))
    ef <- fit$xi_hat * (fit$y - mu(drop(fit$xi_hat %*% beta)))
    for (l in seq_len(K)) for (j in seq_len(K)) {
      dummy <- numeric(K - 1L)
      if (j > 1L) dummy[j - 1L] <- 1
      xi <- cbind(matrix(rep(dummy, each = n), n), q)
      ef <- ef + b[j + K * (l - 1L)] * xi * (mus[, j] - mus[, l])
    }
    ef
  }
  b <- as.numeric(d$B)
  phi <- score(psi, b)
  J <- .csv_jacobian(function(beta) colMeans(score(beta, b)), psi)
  D <- .csv_jacobian(function(b) colMeans(score(psi, b)), b)
  cell <- d$proxy_v + 1L + K * d$v$z
  s <- diag(K * K)[cell, , drop = FALSE] - rep(b, each = nv)
  Sigma <- crossprod(s) / nv
  meat <- crossprod(phi) / n + n / nv * D %*% Sigma %*% t(D)
  if (!is.null(d$v$index)) {
    pv <- phi[d$v$index, , drop = FALSE]
    Gamma <- crossprod(sweep(pv, 2, colMeans(pv)), s) / nv
    meat <- meat + Gamma %*% t(D) + D %*% t(Gamma)
  }
  Ji <- solve(J)
  list(V = Ji %*% meat %*% t(Ji) / n,
       old = Ji %*% (crossprod(phi) / n) %*% t(Ji) / n,
       phi = phi, J = J, D = D, s = s)
}

test_that("validation sandwich matches independent derivatives for all GLM families", {
  for (K in c(2L, 3L)) for (family in c("poisson", "binomial", "gaussian")) {
    for (type in c("external", "internal")) {
      d <- .csv_fixture(K, type, family)
      old <- do.call(mcglm, d$args)
      fit <- do.call(mcglm, c(d$args, list(validation = d$v)))
      ref <- .csv_reference(fit, d)
      expect_equal(unname(vcov(fit, method = "cs")), ref$V, tolerance = 1e-7)
      expect_equal(unname(vcov(old, method = "cs")), ref$old, tolerance = 1e-7)
      expect_identical(fit$coefficients, old$coefficients)
      for (method in c("naive", "bca", "bcm"))
        expect_identical(vcov(fit, method = method), vcov(old, method = method))
    }
  }
})

test_that("binary eta influence formula equals joint-cell implementation", {
  d <- .csv_fixture()
  fit <- do.call(mcglm, c(d$args, list(validation = d$v)))
  psi <- unname(coef(fit, method = "cs"))
  pi <- d$args$pi_z
  p01 <- d$args$Pi[2, 1]
  p10 <- d$args$Pi[1, 2]
  z <- d$v$z
  zh <- d$proxy_v
  nv <- length(z)
  # Equations s-binary-validation and D-binary-explicit, verbatim algebra.
  s <- cbind(z - pi, (z == 0) * (zh - p01) / (1 - pi),
             (z == 1) * ((1 - zh) - p10) / pi)
  eta <- drop(fit$x %*% psi[-1])
  delta <- fit$family$linkinv(eta + psi[1]) - fit$family$linkinv(eta)
  dq <- colMeans(delta * fit$x)
  D <- rbind(c(-p01 * mean(delta), (1 - pi) * mean(delta), 0),
             cbind(-(p01 + p10) * dq, (1 - pi) * dq, -pi * dq))
  ref <- .csv_reference(fit, d)
  pv <- ref$phi[d$v$index, , drop = FALSE]
  Gamma <- crossprod(sweep(pv, 2, colMeans(pv)), s) / nv
  meat <- crossprod(ref$phi) / fit$n + fit$n / nv * D %*% (crossprod(s) / nv) %*% t(D) +
    Gamma %*% t(D) + D %*% t(Gamma)
  Ji <- solve(ref$J)
  expect_equal(unname(vcov(fit, method = "cs")), Ji %*% meat %*% t(Ji) / fit$n,
               tolerance = 1e-7)
})

test_that("internal covariance agrees with a stacked estimating-equation sandwich", {
  for (nv in c(300L, 700L)) {
    d <- .csv_fixture(nv = nv)
    fit <- do.call(mcglm, c(d$args, list(validation = d$v)))
    ref <- .csv_reference(fit, d)
    nuisance <- matrix(0, fit$n, ncol(ref$s))
    nuisance[d$v$index, ] <- fit$n / nv * ref$s
    ef <- cbind(ref$phi, nuisance)
    # Upper block of inverse stacked Jacobian [[J,D],[0,-I]].
    Ji <- solve(ref$J)
    top <- cbind(Ji, Ji %*% ref$D)
    stacked <- top %*% crossprod(ef) %*% t(top) / fit$n^2
    expect_equal(unname(vcov(fit, method = "cs")), stacked, tolerance = 1e-7)
    ext <- list(z = d$v$z, z_hat = d$proxy_v)
    independent <- do.call(mcglm, c(d$args, list(validation = ext)))
    expect_gt(max(abs(vcov(fit, method = "cs") - vcov(independent, method = "cs"))), 1e-5)
  }
})

test_that("validation SE propagates to summary and confidence intervals", {
  d <- .csv_fixture()
  fit <- do.call(mcglm, c(d$args, list(validation = d$v)))
  se <- sqrt(diag(vcov(fit, method = "cs")))
  expect_equal(fit$se$cs, se)
  ci <- confint(fit, method = "cs")
  expect_equal(unname(ci[, 2] - ci[, 1]), unname(2 * qnorm(.975) * se))
  expect_equal(unname(summary(fit)$coefficients$cs[, "Std. Error"]), unname(se))
})

test_that("validation metadata rejects inconsistent probabilities and unsupported designs", {
  d <- .csv_fixture()
  fit <- function(v = d$v, args = d$args) do.call(mcglm, c(args, list(validation = v)))
  expect_error(fit(list(z = d$v$z, index = rep(1, length(d$v$z)))), "distinct")
  expect_error(fit(list(z = d$v$z, index = d$v$index + 1000L)), "row numbers")
  expect_error(fit(list(z = d$v$z, index = d$v$index, z_hat = 1 - d$proxy_v)), "does not match")
  # factor labels "0"/"1" match the integer-coded model by label
  expect_equal(vcov(fit(list(z = factor(d$v$z), index = d$v$index)),
                    method = "cs"),
               vcov(fit(), method = "cs"))
  expect_error(fit(list(z = c(NA, d$v$z[-1]), index = d$v$index)), "codes")
  expect_error(fit(list(z = rep(0, length(d$v$z)), index = d$v$index)), "Every true category")
  expect_error(fit(list(z = d$v$z)), "z_hat")
  wrong <- d$args
  wrong$pi_z <- .8
  expect_error(fit(args = wrong), "empirical estimates")
  wrong <- d$args
  wrong$weights <- rep(2, length(wrong$formula))
  expect_error(fit(args = wrong), "unweighted")
  wrong$weights[] <- 1
  expect_equal(vcov(fit(args = wrong), method = "cs"), vcov(fit(), method = "cs"))
  wrong <- d$args
  wrong$method <- "bcm"
  expect_error(fit(args = wrong), "requires method")
})

test_that("formula, binary c1/c2, and multicategory numerical Jacobian interfaces work", {
  d <- .csv_fixture()
  a <- d$args
  df <- data.frame(y = a$formula, z = a$z_hat, x = a$x[, 2])
  fitted <- mcglm(y ~ mc(z) + x, data = df, family = a$family, method = "cs",
                  Pi = a$Pi, pi_z = a$pi_z, validation = d$v)
  direct <- do.call(mcglm, c(a, list(validation = d$v)))
  expect_equal(unname(vcov(fitted, method = "cs")), unname(vcov(direct, method = "cs")))
  a$Pi <- a$pi_z <- NULL
  a$c1 <- d$B[2, 1]
  a$c2 <- d$B[2, 1] - d$B[1, 2]
  constants <- do.call(mcglm, c(a, list(validation = d$v)))
  expect_equal(vcov(constants, method = "cs"), vcov(direct, method = "cs"))
  d3 <- .csv_fixture(3L)
  analytical <- do.call(mcglm, c(d3$args, list(validation = d3$v)))
  numerical <- do.call(mcglm, c(d3$args, list(validation = d3$v, jacobian = "numerical")))
  expect_equal(vcov(analytical, method = "cs"), vcov(numerical, method = "cs"), tolerance = 1e-5)
})

test_that("external validation uses 1/n_validation scaling even when n_validation > n", {
  d <- .csv_fixture(type = "external", nv = 900L)
  old <- do.call(mcglm, d$args)
  fit <- do.call(mcglm, c(d$args, list(validation = d$v)))
  doubled <- lapply(d$v, rep, times = 2L)
  fit2 <- do.call(mcglm, c(d$args, list(validation = doubled)))
  base <- vcov(old, method = "cs")
  expect_equal(vcov(fit2, method = "cs") - base,
               (vcov(fit, method = "cs") - base) / 2, tolerance = 1e-12)
})

test_that("zero empirical misclassification adds no validation uncertainty", {
  for (K in c(2L, 3L)) {
    d <- .csv_fixture(K, type = "external")
    d$v$z_hat <- d$v$z
    d$args$Pi <- diag(K)
    prevalence <- tabulate(d$v$z + 1L, nbins = K) / length(d$v$z)
    d$args$pi_z <- if (K == 2L) prevalence[2] else prevalence
    old <- do.call(mcglm, d$args)
    fit <- do.call(mcglm, c(d$args, list(validation = d$v)))
    expect_equal(vcov(fit, method = "cs"), vcov(old, method = "cs"), tolerance = 1e-12)
  }
})
