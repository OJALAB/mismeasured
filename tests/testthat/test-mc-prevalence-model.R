# Covariate-dependent prevalence P(Z | x) (control_mc(prevalence_model = )).

.pm_data <- function(n = 2000L, nv = 400L, seed = 801L, a = -0.5, b = 1.2,
                     type = "internal") {
  set.seed(seed)
  df <- data.frame(x1 = rnorm(n))
  z <- rbinom(n, 1L, plogis(a + b * df$x1))
  df$z <- ifelse(z == 1L, rbinom(n, 1L, 0.85), rbinom(n, 1L, 0.12))
  df$y <- rpois(n, exp(-0.3 + 0.8 * z + 0.5 * df$x1))
  if (type == "internal") {
    idx <- sample.int(n, nv)
    val <- validation_sample(z[idx], index = idx)
  } else {
    dv <- data.frame(x1 = rnorm(nv))
    zv <- rbinom(nv, 1L, plogis(a + b * dv$x1))
    zhv <- ifelse(zv == 1L, rbinom(nv, 1L, 0.85), rbinom(nv, 1L, 0.12))
    val <- validation_sample(zv, zhv, data = dv)
  }
  list(df = df, z = z, val = val)
}
.pm <- control_mc(prevalence_model = ~ x1)

test_that("multinom wrapper and prevalence probabilities", {
  set.seed(1)
  x <- rnorm(500)
  z <- rbinom(500, 1, plogis(-0.3 + x))
  w <- runif(500, 0.5, 2)
  a <- mismeasured:::.mc_multinom(z, cbind(1, x), w, 2L)
  ref <- coef(suppressWarnings(glm(z ~ x, binomial, weights = w)))
  expect_equal(as.numeric(a), unname(ref), tolerance = 1e-5)
  P <- mismeasured:::.mc_prev_probs(a, cbind(1, x))
  expect_equal(P[, 2], unname(plogis(drop(cbind(1, x) %*% ref))),
               tolerance = 1e-5)
  expect_equal(rowSums(P), rep(1, 500))
})

test_that("estimate_mc() with a prevalence model: validation and em", {
  d <- .pm_data()
  est <- estimate_mc(d$val, z_hat = d$df$z, data = d$df, control = .pm)
  vb <- est$validation
  ref <- suppressWarnings(glm(vb$z ~ d$df$x1[vb$index], family = binomial))
  alpha <- est$map(est$eta)$alpha
  expect_equal(as.numeric(alpha), unname(coef(ref)), tolerance = 1e-6)
  # sandwich of the logistic fit for alpha
  X <- cbind(1, d$df$x1[vb$index])
  pr <- plogis(drop(X %*% coef(ref)))
  A <- crossprod(X * (pr * (1 - pr)), X)
  B <- crossprod(X * (vb$z - pr))
  V <- solve(A) %*% B %*% solve(A)
  expect_equal(unname(vcov(est)[3:4, 3:4]), V, tolerance = 1e-5)
  expect_named(est$eta, c("Pi[1,0]", "Pi[1,1]", "alpha[1,(Intercept)]",
                          "alpha[1,x1]"))
  expect_output(print(est), "Prevalence model ~x1")

  # em: score of sum log sum_l Pi[zhat, l] pi_l(x) is zero
  e <- .pm_data(type = "external", seed = 802L)
  ext <- validation_sample(e$val$z, e$val$z_hat)
  est_em <- estimate_mc(ext, z_hat = e$df$z, data = e$df,
                        control = control_mc(prevalence = "em",
                                             prevalence_model = ~ x1))
  a <- est_em$map(est_em$eta)$alpha
  X <- cbind(1, e$df$x1)
  P <- cbind(1 - plogis(drop(X %*% t(a))), plogis(drop(X %*% t(a))))
  joint <- unname(est_em$Pi)[e$df$z + 1L, ] * P
  post <- joint / rowSums(joint)
  expect_lt(max(abs(colSums((post[, 2] - P[, 2]) * X))), 1e-6)
})

test_that("sub with a prevalence model matches an independent stacked reference", {
  for (type in c("internal", "external")) {
    d <- .pm_data(n = 800L, nv = 250L, type = type, seed = 810L)
    fit <- mcglm(y ~ mc(z) + x1, data = d$df, method = "sub",
                 validation = d$val, mc_control = .pm)
    psi <- unname(coef(fit, method = "sub"))
    est <- fit$mc_estimate
    vb <- est$validation
    n <- nrow(d$df)
    X <- cbind(1, d$df$x1)
    Xv <- if (type == "internal") X[vb$index, ] else cbind(1, d$val$data$x1)
    zhv <- vb$proxy
    rows <- function(full) {
      b <- full[1:3]
      Pi <- rbind(1 - full[4:5], full[4:5])
      a <- full[6:7]
      P1 <- plogis(drop(X %*% a))
      prior <- cbind(1 - P1, P1)
      joint <- Pi[d$df$z + 1L, ] * prior
      Wi <- joint / rowSums(joint)
      mu <- exp(outer(drop(X %*% b[2:3]), c(0, b[1]), "+"))
      U <- cbind(d$df$z, X) * (d$df$y - rowSums(Wi * mu))
      nuis <- cbind((vb$z == 0) * (zhv - Pi[2, 1]), (vb$z == 1) * (zhv - Pi[2, 2]),
                    (vb$z - plogis(drop(Xv %*% a))) * Xv)
      if (type == "internal") {
        idx <- vb$index
        xi <- cbind(vb$z, X[idx, ])
        U[idx, ] <- xi * (d$df$y[idx] - exp(drop(xi %*% b)))
        G <- cbind(U, matrix(0, n, 4))
        G[idx, 4:7] <- nuis
        G
      } else {
        rbind(cbind(U, matrix(0, n, 4)), cbind(matrix(0, vb$n, 3), nuis))
      }
    }
    full <- c(psi, unname(est$eta))
    G <- rows(full)
    expect_lt(max(abs(colSums(G))), 1e-6 * n)
    A <- vapply(seq_along(full), function(k) {
      h <- 1e-6; e <- replace(numeric(7), k, h)
      (colSums(rows(full + e)) - colSums(rows(full - e))) / (2 * h)
    }, numeric(7))
    V <- solve(A) %*% crossprod(G) %*% t(solve(A))
    expect_equal(unname(vcov(fit, method = "sub")), V[1:3, 1:3],
                 tolerance = 1e-5, label = type)
  }
})

test_that("il with a prevalence model maximises the joint likelihood", {
  d <- .pm_data(n = 700L, nv = 220L, seed = 820L)
  fit <- mcglm(y ~ mc(z) + x1, data = d$df, method = "il",
               validation = d$val, mc_control = .pm)
  vb <- fit$mc_estimate$validation
  X <- cbind(1, d$df$x1)
  idx <- vb$index
  Pi <- unname(fit$nuisance$il$Pi)
  full <- c(unname(coef(fit, method = "il")),
            log(Pi[2, ] / Pi[1, ]), as.numeric(fit$nuisance$il$alpha))
  unit_ll <- function(f) {
    b <- f[1:3]
    P <- rbind(1 / (1 + exp(f[4:5])), exp(f[4:5]) / (1 + exp(f[4:5])))
    p1 <- plogis(drop(X %*% f[6:7]))
    prior <- cbind(1 - p1, p1)
    lf <- matrix(dpois(d$df$y, exp(outer(drop(X %*% b[2:3]), c(0, b[1]), "+")),
                       log = TRUE), ncol = 2)
    ll <- log(rowSums(exp(lf) * P[d$df$z + 1L, ] * prior))
    ll[idx] <- lf[cbind(idx, vb$z + 1L)] + log(P[cbind(d$df$z[idx] + 1L, vb$z + 1L)]) +
      log(prior[cbind(idx, vb$z + 1L)])
    ll
  }
  grad <- function(f, h = 1e-5) vapply(seq_along(f), function(k) {
    e <- replace(numeric(length(f)), k, h)
    (unit_ll(f + e) - unit_ll(f - e)) / (2 * h)
  }, numeric(nrow(d$df)))
  G <- grad(full)
  expect_lt(max(abs(colSums(G))), 1e-4 * nrow(d$df))
  A <- vapply(seq_along(full), function(k) {
    h <- 1e-4; e <- replace(numeric(7), k, h)
    (colSums(grad(full + e)) - colSums(grad(full - e))) / (2 * h)
  }, numeric(7))
  V <- solve(A) %*% crossprod(G) %*% t(solve(A))
  expect_equal(unname(vcov(fit, method = "il")), V[1:3, 1:3], tolerance = 2e-3)
})

test_that("a prevalence model removes the bias when Z depends on x", {
  skip_on_cran()
  d <- .pm_data(n = 20000L, nv = 1500L, seed = 830L)
  expect_warning(
    f0 <- mcglm(y ~ mc(z) + x1, data = d$df, method = c("sub", "il"),
                validation = d$val),
    "associated with the covariates")
  expect_lt(f0$independence_test$p.value, 1e-6)
  f <- mcglm(y ~ mc(z) + x1, data = d$df, method = c("sub", "ec", "il"),
             validation = d$val, mc_control = .pm)
  truth <- c(0.8, -0.3, 0.5)
  expect_gt(abs(coef(f0, method = "sub")[3] - 0.5), 0.06)
  for (m in c("sub", "ec", "il"))
    expect_lt(max(abs(coef(f, method = m) - truth)), 0.05, label = m)
  expect_null(f$independence_test)
})

test_that("prevalence models: factor covariates and K = 3", {
  set.seed(840)
  n <- 1500
  df <- data.frame(x1 = rnorm(n), region = factor(sample(c("n", "s", "w"), n, TRUE)))
  eta <- cbind(0, 0.3 + 0.8 * df$x1, -0.4 + 0.7 * (df$region == "s"))
  pr <- exp(eta) / rowSums(exp(eta))
  z <- apply(pr, 1, function(p) sample(0:2, 1, prob = p))
  Pi <- matrix(c(0.85, 0.10, 0.05, 0.10, 0.80, 0.10, 0.05, 0.10, 0.85), 3L)
  df$occ <- factor(c("a", "b", "c")[vapply(z, function(k)
    sample.int(3L, 1L, prob = Pi[, k + 1L]), 0L)])
  df$y <- rpois(n, exp(c(0, 0.6, -0.4)[z + 1L] - 0.3 + 0.5 * df$x1))
  idx <- sample.int(n, 500)
  fit <- mcglm(y ~ mc(occ) + x1 + region, data = df,
               method = c("sub", "ec", "il", "cs_akn"),
               validation = validation_sample(c("a", "b", "c")[z[idx] + 1L],
                                              index = idx),
               mc_control = control_mc(prevalence_model = ~ x1 + region))
  for (m in c("sub", "ec", "il", "cs_akn"))
    expect_true(all(is.finite(fit$se[[m]])), label = m)
  expect_identical(dim(fit$nuisance$il$alpha), c(2L, 4L))
  expect_identical(names(fit$mc_estimate$pi), c("a", "b", "c"))
})

test_that("prevalence model input errors", {
  d <- .pm_data(n = 300L, nv = 100L, seed = 850L)
  expect_error(control_mc(prevalence = "inverse", prevalence_model = ~ x1),
               "no covariate-dependent")
  expect_error(control_mc(prevalence_model = y ~ x1), "one-sided")
  expect_error(mcglm(y ~ mc(z) + x1, data = d$df, method = c("sub", "cs"),
                     validation = d$val, mc_control = .pm),
               "constant prevalence")
  expect_error(mcglm(d$df$y, z_hat = d$df$z, x = cbind(1, d$df$x1),
                     method = "sub", validation = d$val, mc_control = .pm),
               "needs the data")
  e <- .pm_data(n = 300L, nv = 100L, seed = 851L, type = "external")
  no_data <- validation_sample(e$val$z, e$val$z_hat)
  expect_error(mcglm(y ~ mc(z) + x1, data = e$df, method = "sub",
                     validation = no_data, mc_control = .pm), "needs the data")
  expect_error(validation_sample(c(0, 1), index = 1:2,
                                 data = data.frame(x1 = 1:2)), "external")
})

test_that("prevalence-model fits are unbiased and calibrated", {
  skip_on_cran()
  reps <- 60L
  for (type in c("internal", "external")) {
    out <- matrix(NA_real_, reps, 6L)
    ctrl <- if (type == "internal") .pm else
      control_mc(prevalence = "em", prevalence_model = ~ x1)
    for (r in seq_len(reps)) {
      d <- .pm_data(n = 2000L, nv = 400L, seed = 9000L + r, type = type)
      val <- if (type == "internal") d$val else
        validation_sample(d$val$z, d$val$z_hat)
      fit <- mcglm(y ~ mc(z) + x1, data = d$df, method = c("sub", "ec", "il"),
                   validation = val, mc_control = ctrl)
      out[r, ] <- c(sapply(c("sub", "ec", "il"), function(m)
        c(coef(fit, method = m)[1], fit$se[[m]][1])))
    }
    for (k in c(1L, 3L, 5L)) {
      expect_lt(abs(mean(out[, k]) - 0.8), 0.04, label = paste(type, k))
      ratio <- mean(out[, k + 1L]) / sd(out[, k])
      expect_true(ratio > 0.75 && ratio < 1.3, label = paste(type, k, ratio))
    }
  }
})

test_that("the independence warning follows control_mc(on_ill)", {
  d <- .pm_data(n = 3000L, nv = 600L, seed = 860L, b = 2)
  expect_error(mcglm(y ~ mc(z) + x1, data = d$df, method = "sub",
                     validation = d$val,
                     mc_control = control_mc(on_ill = "error")),
               "sub assumes a constant prevalence")
  fit <- mcglm(y ~ mc(z) + x1, data = d$df, method = c("sub", "cs_akn"),
               validation = d$val, mc_control = control_mc(on_ill = "none"))
  expect_lt(fit$independence_test$p.value, 0.01)
  # cs_akn alone does not depend on the prevalence: no test
  fit2 <- mcglm(y ~ mc(z) + x1, data = d$df, method = "cs_akn",
                validation = d$val)
  expect_null(fit2$independence_test)
})
