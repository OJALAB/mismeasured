# estimate_mc(), control_mc(), EM prevalence and diagnose_mc().

.mce_labels <- function(m, Pi, prev) {
  K <- nrow(Pi)
  z <- sample.int(K, m, TRUE, prev) - 1L
  zh <- vapply(z, function(k) sample.int(K, 1L, prob = Pi[, k + 1L]) - 1L, 0L)
  list(z = z, z_hat = zh)
}
.mce_Pi2 <- matrix(c(0.85, 0.15, 0.10, 0.90), 2L)
.mce_Pi3 <- matrix(c(0.85, 0.10, 0.05,
                     0.10, 0.80, 0.10,
                     0.05, 0.10, 0.85), 3L)

test_that("EM prevalence equals the inversion when interior", {
  Pi <- .mce_Pi3
  pi <- c(0.5, 0.3, 0.2)
  em <- mismeasured:::.mc_prevalence_em(drop(Pi %*% pi), Pi)
  expect_equal(em$pi, pi, tolerance = 1e-10)
  expect_false(em$boundary)
  # unchanged interface of the mcglm helper
  z_hat <- rep(0:1, round(1000 * drop(.mce_Pi2 %*% c(0.6, 0.4))))
  expect_equal(mismeasured:::.mcglm_estimate_pi_z(z_hat, .mce_Pi2), 0.4,
               tolerance = 1e-12)
})

test_that("EM prevalence returns the boundary MLE instead of clamping", {
  Pi <- .mce_Pi2
  p <- c(0.95, 0.05)              # below the false-positive rate 0.15
  expect_true(solve(Pi, p)[2] < 0)
  em <- mismeasured:::.mc_prevalence_em(p, Pi)
  expect_true(em$boundary)
  expect_equal(sum(em$pi), 1)
  expect_true(all(em$pi >= 0))
  loglik <- function(pi) sum(p * log(drop(Pi %*% pi)))
  expect_gte(loglik(em$pi), loglik(c(0.99, 0.01)))
  z_hat <- rep(0:1, c(95, 5))
  expect_warning(pz <- mismeasured:::.mcglm_estimate_pi_z(z_hat, Pi),
                 "boundary")
  expect_lt(pz, 1e-4)
})

test_that("unweighted estimates are proportions with binomial variances", {
  set.seed(1)
  v <- .mce_labels(400, .mce_Pi2, c(0.6, 0.4))
  est <- estimate_mc(validation_sample(v$z, v$z_hat))
  tab <- table(factor(v$z_hat, 0:1), factor(v$z, 0:1))
  expect_equal(unname(est$Pi), unname(unclass(prop.table(tab, 2))))
  pi1 <- mean(v$z)
  expect_equal(est$pi, c(1 - pi1, pi1))
  n_l <- colSums(tab)
  se <- sqrt(diag(vcov(est)))
  expect_equal(unname(se[1:2]),
               unname(sqrt(est$Pi[2, ] * (1 - est$Pi[2, ]) / n_l)),
               tolerance = 1e-8)
  expect_equal(unname(se[3]), sqrt(pi1 * (1 - pi1) / 400), tolerance = 1e-8)
  joint <- sweep(est$Pi, 2, est$pi, "*")
  expect_equal(est$W, joint / rowSums(joint))
  expect_s3_class(summary(est), "summary.mc_estimate")
  expect_output(print(est), "Pi = P\\(z_hat")
  expect_output(print(summary(est)), "No problems detected")
})

test_that("design-weighted Hajek and HT estimators", {
  set.seed(2)
  n <- 3000
  main <- .mce_labels(n, .mce_Pi2, c(0.6, 0.4))
  # Poisson audit of the main study, oversampling predicted positives
  p_incl <- ifelse(main$z_hat == 1, 0.2, 0.05)
  idx <- which(runif(n) < p_incl)
  d <- 1 / p_incl[idx]
  z_v <- main$z[idx]
  vh <- estimate_mc(validation_sample(z_v, index = idx, weights = d),
                    z_hat = main$z_hat)
  vt <- estimate_mc(validation_sample(z_v, index = idx, weights = d,
                                      estimator = "ht"),
                    z_hat = main$z_hat)
  # column ratios are the same; prevalence differs
  expect_equal(vh$Pi, vt$Pi)
  expect_equal(vh$pi[2], sum(d * z_v) / sum(d))
  expect_equal(vt$pi[2], sum(d * z_v) / n)
  # the unweighted column ratios are biased by the oversampling of z_hat = 1
  unw <- estimate_mc(validation_sample(z_v, index = idx),
                     z_hat = main$z_hat)
  expect_gt(abs(unw$Pi[2, 1] - 0.15), abs(vh$Pi[2, 1] - 0.15))
  expect_true(all(is.finite(vt$vcov)))
})

test_that("z_hat-stratified audit gives within-stratum predictive values", {
  set.seed(3)
  n <- 2000
  main <- .mce_labels(n, .mce_Pi3, c(0.5, 0.3, 0.2))
  idx <- unlist(lapply(0:2, function(j) sample(which(main$z_hat == j), 40)))
  est <- estimate_mc(validation_sample(main$z[idx], index = idx,
                                      strata = "z_hat"),
                     z_hat = main$z_hat)
  within <- prop.table(table(factor(main$z_hat[idx], 0:2),
                             factor(main$z[idx], 0:2)), 1)
  expect_equal(unname(est$W), unname(unclass(within)), tolerance = 1e-12)
  # external z_hat strata need the frame sizes
  v <- validation_sample(main$z[idx], main$z_hat[idx], strata = "z_hat")
  expect_error(estimate_mc(v), "N_strata")
  v$N_strata <- c(`0` = 900, `1` = 600, `2` = 500)
  ext <- estimate_mc(v)
  expect_true(all(is.finite(ext$vcov)))
})

test_that("prevalence from the main study: em and inverse", {
  set.seed(4)
  main <- .mce_labels(2000, .mce_Pi2, c(0.6, 0.4))
  v <- .mce_labels(300, .mce_Pi2, c(0.6, 0.4))
  val <- validation_sample(v$z, v$z_hat)
  expect_error(estimate_mc(val, control = control_mc(prevalence = "em")),
               "needs the main-study proxies")
  em <- estimate_mc(val, z_hat = main$z_hat,
                    control = control_mc(prevalence = "em"))
  inv <- estimate_mc(val, z_hat = main$z_hat,
                     control = control_mc(prevalence = "inverse"))
  expect_equal(em$pi, inv$pi, tolerance = 1e-8)
  p_main <- mean(main$z_hat)
  expect_equal(drop(em$Pi %*% em$pi)[2], p_main, tolerance = 1e-8)
  # the two parameterisations give the same delta-method variance for pi
  g <- mismeasured:::.mc_num_jacobian(function(e) inv$map(e)$pi[2],
                                      unname(inv$eta))
  expect_equal(drop(g %*% vcov(inv) %*% t(g)), vcov(em)["pi[1]", "pi[1]"],
               tolerance = 1e-6)
})

test_that("diagnose_mc() flags ill-conditioned matrices", {
  # Pi-hat = [.51 .49; .49 .51]: singular values 1 and 0.02
  z <- rep(0:1, each = 100)
  z_hat <- c(rep(0:1, c(51, 49)), rep(0:1, c(49, 51)))
  val <- validation_sample(z, z_hat)
  expect_warning(est <- estimate_mc(val), "singular value")
  expect_true(length(est$diagnostics$problems) > 0)
  expect_error(estimate_mc(val, control = control_mc(on_ill = "error")),
               "singular value")
  expect_silent(estimate_mc(val, control = control_mc(on_ill = "none")))
  small <- validation_sample(c(0, 0, 0, 1, 1, 1), c(0, 0, 1, 1, 1, 0))
  expect_warning(estimate_mc(small), "fewer than 10")
  expect_output(print(diagnose_mc(estimate_mc(
    small, control = control_mc(on_ill = "none")))), "Problems")
})

test_that("delta-method SEs of design-weighted estimates are calibrated", {
  skip_on_cran()
  set.seed(6)
  n <- 3000
  reps <- 300
  out <- matrix(NA_real_, reps, 6)
  for (r in seq_len(reps)) {
    main <- .mce_labels(n, .mce_Pi2, c(0.6, 0.4))
    p_incl <- ifelse(main$z_hat == 1, 0.15, 0.04)
    idx <- which(runif(n) < p_incl)
    est <- estimate_mc(validation_sample(main$z[idx], index = idx,
                                         weights = 1 / p_incl[idx]),
                       z_hat = main$z_hat,
                       control = control_mc(on_ill = "none"))
    out[r, ] <- c(est$eta, sqrt(diag(vcov(est))))
  }
  expect_lt(abs(mean(out[, 1]) - 0.15), 0.01)
  expect_lt(abs(mean(out[, 3]) - 0.4), 0.01)
  ratio <- colMeans(out[, 4:6]) / apply(out[, 1:3], 2, sd)
  expect_true(all(ratio > 0.85 & ratio < 1.2))
})
