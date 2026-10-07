# Shrinkage estimators of Pi, regularisation, bootstrap and posterior
# variances.

.sh_sparse <- function(K = 6L, n = 6000L, nv = 150L, seed = 1301L,
                       type = "internal") {
  set.seed(seed)
  Pi <- matrix(0.03, K, K)
  diag(Pi) <- 1 - 0.03 * (K - 1L)
  z <- sample.int(K, n, TRUE) - 1L
  zh <- vapply(z, function(k) sample.int(K, 1L, prob = Pi[, k + 1L]) - 1L, 0L)
  x1 <- rnorm(n)
  gam <- c(0, 0.5, -0.3, 0.2, 0.4, -0.5, 0.3, -0.2)[seq_len(K)]
  y <- rpois(n, exp(gam[z + 1L] - 0.3 + 0.5 * x1))
  if (type == "internal") {
    idx <- sample.int(n, nv)
    val <- validation_sample(z[idx], index = idx)
  } else {
    zv <- sample.int(K, nv, TRUE) - 1L
    val <- validation_sample(zv, vapply(zv, function(k)
      sample.int(K, 1L, prob = Pi[, k + 1L]) - 1L, 0L))
  }
  list(y = y, z = z, z_hat = zh, x = cbind(1, x1), Pi = Pi, val = val,
       psi = c(gam[-1L], -0.3, 0.5))
}

test_that("eb limits, targets and refusals", {
  d <- .sh_sparse()
  est0 <- estimate_mc(d$val, z_hat = d$z_hat,
                      control = control_mc(on_ill = "none"))
  zero <- estimate_mc(d$val, z_hat = d$z_hat,
                      control = control_mc(estimator = "eb", concentration = 0,
                                           on_ill = "none"))
  expect_equal(zero$Pi, est0$Pi)
  big <- estimate_mc(d$val, z_hat = d$z_hat,
                     control = control_mc(estimator = "eb",
                                          concentration = 1e9,
                                          on_ill = "none"))
  expect_equal(unname(big$Pi), big$shrinkage$T, tolerance = 1e-5)
  eb <- estimate_mc(d$val, z_hat = d$z_hat,
                    control = control_mc(estimator = "eb", on_ill = "none"))
  expect_equal(unname(colSums(eb$Pi)), rep(1, 6))
  expect_equal(sum(eb$Pi == 0), 0)
  expect_lt(max(abs(eb$Pi - d$Pi)), max(abs(est0$Pi - d$Pi)))
  # pooled target: one diagonal value, equal off-diagonal values
  T <- eb$shrinkage$T
  expect_equal(length(unique(round(diag(T), 12))), 1L)
  expect_lt(max(abs(eb$score(eb$eta))), 1e-8)
  # groups: within-group errors differ from between-group errors
  gr <- estimate_mc(d$val, z_hat = d$z_hat,
                    control = control_mc(estimator = "eb", target = "groups",
                                         groups = rep(1:2, each = 3),
                                         on_ill = "none"))
  expect_equal(unname(colSums(gr$shrinkage$T)), rep(1, 6))
  ll <- estimate_mc(d$val, z_hat = d$z_hat,
                    control = control_mc(estimator = "eb",
                                         target = "loglinear",
                                         on_ill = "none"))
  expect_equal(unname(colSums(ll$shrinkage$T)), rep(1, 6))
  # refusals: uniform-like pooled target, bad matrices, missing pieces
  flat <- validation_sample(rep(0:1, each = 50), rep(0:1, 50))
  expect_error(estimate_mc(flat, control = control_mc(estimator = "eb",
                                                      on_ill = "none")),
               "ill-conditioned")
  expect_error(control_mc(estimator = "eb", target = "matrix"),
               "target_matrix")
  expect_error(estimate_mc(d$val, z_hat = d$z_hat,
                           control = control_mc(estimator = "eb",
                                                target = "matrix",
                                                target_matrix = diag(3))),
               "K x K")
  expect_error(estimate_mc(d$val, z_hat = d$z_hat,
                           control = control_mc(estimator = "eb",
                                                target = "groups")),
               "groups")
  expect_warning(estimate_mc(d$val, z_hat = d$z_hat,
                             control = control_mc(estimator = "dirichlet",
                                                  on_ill = "none")),
                 "uniform matrix")
})

test_that("the marginal-likelihood concentration is of the right size", {
  set.seed(1310)
  K <- 8L
  c0 <- 40
  T <- matrix(0.2 / (K - 1), K, K)
  diag(T) <- 0.8
  C <- vapply(seq_len(K), function(l) {
    p <- rgamma(K, c0 * T[, l]); p <- p / sum(p)
    tabulate(sample.int(K, 80L, TRUE, p), nbins = K)
  }, numeric(K))
  cc <- mismeasured:::.mc_shrink_concentration(C, T)
  expect_gt(cc, c0 / 3)
  expect_lt(cc, c0 * 3)
})

test_that("shrunk sub sandwich matches a penalised stacked reference", {
  d <- .sh_sparse(K = 3L, n = 900L, nv = 120L, type = "external", seed = 1320L)
  ctrl <- control_mc(estimator = "eb", concentration = 25, on_ill = "none")
  fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, method = "sub",
               validation = d$val, mc_control = ctrl)
  est <- fit$mc_estimate
  T <- est$shrinkage$T
  K <- 3L; n <- length(d$y)
  zv <- est$validation$z; zhv <- est$validation$proxy
  psi <- unname(coef(fit, method = "sub"))
  xi_hat <- cbind(outer(d$z_hat, 1:2, "==") * 1, d$x)
  data_rows <- function(f) {
    b <- f[1:4]
    Pi <- matrix(f[5:10], 2, 3); Pi <- rbind(1 - colSums(Pi), Pi)
    pr <- c(1 - sum(f[11:12]), f[11:12])
    W <- t(t(Pi) * pr); W <- W / rowSums(W)
    mu <- exp(outer(drop(d$x %*% b[3:4]), c(0, b[1:2]), "+"))
    U <- xi_hat * (d$y - rowSums(W[d$z_hat + 1L, ] * mu))
    nuis <- cbind(do.call(cbind, lapply(0:2, function(l)
      sapply(1:2, function(j) (zv == l) * ((zhv == j) - Pi[j + 1L, l + 1L])))),
      sapply(1:2, function(l) (zv == l) - pr[l + 1L]))
    rbind(cbind(U, matrix(0, n, 8)), cbind(matrix(0, length(zv), 4), nuis))
  }
  total <- function(f) {
    Pi <- matrix(f[5:10], 2, 3); Pi <- rbind(1 - colSums(Pi), Pi)
    out <- colSums(data_rows(f))
    out[5:10] <- out[5:10] + as.numeric(25 * (T[-1, ] - Pi[-1, ]))
    out
  }
  full <- c(psi, unname(est$eta))
  expect_lt(max(abs(total(full))), 1e-6 * n)
  A <- vapply(seq_along(full), function(k) {
    h <- 1e-6; e <- replace(numeric(12), k, h)
    (total(full + e) - total(full - e)) / (2 * h)
  }, numeric(12))
  V <- solve(A) %*% crossprod(data_rows(full)) %*% t(solve(A))
  expect_equal(unname(vcov(fit, method = "sub")), V[1:4, 1:4], tolerance = 1e-5)
})

test_that("il uses the Dirichlet prior and converges on sparse audits", {
  skip_on_cran()
  d <- .sh_sparse()
  fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, method = c("sub", "il"),
               validation = d$val,
               mc_control = control_mc(estimator = "eb", on_ill = "none"))
  expect_true(fit$convergence$il$converged)
  expect_lt(max(abs(coef(fit, method = "il") - d$psi)), 0.2)
  expect_true(all(is.finite(fit$se$il)))
})

test_that("on_ill = 'regularize' switches to eb", {
  weak <- validation_sample(rep(0:1, each = 100),
                            c(rep(0:1, c(51, 49)), rep(0:1, c(49, 51))))
  msgs <- character()
  est <- withCallingHandlers(
    estimate_mc(weak, control = control_mc(on_ill = "regularize")),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("re-estimating with estimator = 'eb'", msgs)))
  # an uninformative classifier stays flagged after shrinkage
  expect_true(any(grepl("smallest singular value", msgs[-1])))
  expect_identical(est$shrinkage$estimator, "eb")
})

test_that("bootstrap and posterior variances", {
  d <- .sh_sparse(K = 2L, n = 1200L, nv = 250L, type = "external", seed = 1330L)
  df <- data.frame(y = d$y, z = d$z_hat, x1 = d$x[, 2])
  boot <- function(seed) mcglm(y ~ mc(z) + x1, data = df, method = "sub",
                               validation = d$val,
                               mc_control = control_mc(variance = "bootstrap",
                                                       B = 30, seed = seed))
  set.seed(99); u1 <- runif(1)
  set.seed(99); b1 <- boot(7); u2 <- runif(1)
  expect_identical(u1, u2)                 # global RNG state restored
  expect_identical(b1$se$sub, boot(7)$se$sub)
  expect_identical(b1$replicates$type, "bootstrap")
  expect_equal(dim(b1$replicates$coefficients$sub), c(30L, 3L))

  post <- mcglm(y ~ mc(z) + x1, data = df, method = c("sub", "bca", "il"),
                validation = d$val,
                mc_control = control_mc(variance = "posterior", B = 30,
                                        seed = 1))
  expect_setequal(post$replicates$delta_kept, c("bca", "il"))
  dl <- mcglm(y ~ mc(z) + x1, data = df, method = c("sub", "bca", "il"),
              validation = d$val)
  expect_identical(post$se$il, dl$se$il)

  int <- .sh_sparse(K = 2L, n = 600L, nv = 100L, seed = 1331L)
  dfi <- data.frame(y = int$y, z = int$z_hat, x1 = int$x[, 2])
  expect_error(mcglm(y ~ mc(z) + x1, data = dfi, method = "sub",
                     validation = int$val,
                     mc_control = control_mc(variance = "posterior", B = 5)),
               "external audits")
  expect_error(mcglm(y ~ mc(z) + x1, data = df, method = "sub",
                     validation = estimate_mc(d$val),
                     mc_control = control_mc(variance = "bootstrap", B = 5)),
               "validation sample itself")
})

test_that("resampling keeps internal audit membership", {
  val <- validation_sample(c("a", "b", "a"), index = c(2, 5, 7),
                           weights = c(1, 2, 3))
  b <- mismeasured:::.mc_resample_validation(val, c(5, 5, 1, 7, 3))
  expect_identical(b$index, c(1L, 2L, 4L))
  expect_identical(b$z, c("b", "b", "a"))
  expect_identical(b$weights, c(2, 2, 3))
})

test_that("bootstrap and posterior standard errors agree with the delta method", {
  skip_on_cran()
  d <- .sh_sparse(K = 2L, n = 2000L, nv = 300L, type = "external", seed = 1340L)
  df <- data.frame(y = d$y, z = d$z_hat, x1 = d$x[, 2])
  fit <- function(v) mcglm(y ~ mc(z) + x1, data = df, method = c("sub", "cs"),
                           validation = d$val,
                           mc_control = control_mc(variance = v, B = 150,
                                                   seed = 3))
  dl <- fit("delta")
  for (v in c("bootstrap", "posterior")) {
    f <- fit(v)
    for (m in c("sub", "cs")) {
      ratio <- f$se[[m]] / dl$se[[m]]
      expect_true(all(ratio > 0.75 & ratio < 1.33), label = paste(v, m))
    }
  }
})

test_that("validation designs combine with every estimator and variance", {
  skip_on_cran()
  set.seed(1350)
  n <- 900
  lev <- c("a", "b", "c")
  Pi <- matrix(c(0.85, 0.10, 0.05, 0.10, 0.80, 0.10, 0.05, 0.10, 0.85), 3L)
  df <- data.frame(x1 = rnorm(n), region = factor(sample(c("n", "s"), n, TRUE)))
  z <- sample(0:2, n, TRUE, c(0.5, 0.3, 0.2))
  df$occ <- factor(lev[vapply(z, function(k) sample(0:2, 1, prob = Pi[, k + 1]),
                              0L) + 1L], levels = lev)
  df$y <- rpois(n, exp(c(0, 0.6, -0.4)[z + 1] - 0.3 + 0.5 * df$x1))
  zv <- sample(0:2, 250, TRUE, c(0.5, 0.3, 0.2))
  zhv <- vapply(zv, function(k) sample(0:2, 1, prob = Pi[, k + 1]), 0L)
  dv <- data.frame(x1 = rnorm(250), region = factor(sample(c("n", "s"), 250, TRUE)))
  p_incl <- ifelse(df$occ == "a", 0.15, 0.35)
  iw <- which(runif(n) < p_incl)
  is <- unlist(lapply(lev, function(l) sample(which(df$occ == l), 60)))
  vals <- list(
    validation_sample(lev[z[iw] + 1], index = iw, weights = 1 / p_incl[iw],
                      estimator = "ht"),
    validation_sample(lev[z[is] + 1], index = is, strata = "z_hat"),
    validation_sample(lev[zv + 1], lev[zhv + 1], strata = "z_hat",
                      N_strata = c(a = 2000, b = 1500, c = 1000), data = dv))
  ctrls <- list(control_mc(on_ill = "none"),
                control_mc(estimator = "eb", target = "groups",
                           groups = c(1, 1, 2), on_ill = "none"),
                control_mc(prevalence = "em", on_ill = "none"),
                control_mc(prevalence_model = ~ x1 + region, on_ill = "none"),
                control_mc(variance = "bootstrap", B = 10, seed = 1,
                           on_ill = "none"))
  for (v in vals) for (ct in ctrls) {
    m <- if (is.null(ct$prevalence_model))
      c("bca", "bcm", "cs", "sub", "ec", "il", "cs_akn") else
      c("sub", "ec", "il", "cs_akn")
    f <- suppressWarnings(mcglm(y ~ mc(occ) + x1 + region, data = df,
                                method = m, validation = v, mc_control = ct))
    for (mm in m) expect_true(all(is.finite(f$se[[mm]])), label = mm)
  }
})
