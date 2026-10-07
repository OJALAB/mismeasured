# simex() with estimate_mc() objects in mc().

.se_data <- function(seed = 1601L, n = 1200L, nv = 250L) {
  set.seed(seed)
  z <- rbinom(n, 1L, 0.4)
  zh <- ifelse(z == 1L, rbinom(n, 1L, 0.88), rbinom(n, 1L, 0.1))
  x1 <- rnorm(n)
  y <- 1 + 0.8 * z + 0.5 * x1 + rnorm(n)
  zv <- rbinom(nv, 1L, 0.4)
  zhv <- ifelse(zv == 1L, rbinom(nv, 1L, 0.88), rbinom(nv, 1L, 0.1))
  list(df = data.frame(y, z = factor(zh), x1), z = z, zv = zv, zhv = zhv,
       val = validation_sample(zv, zhv))
}

test_that("mc(z, estimate) uses the estimated matrix", {
  d <- .se_data()
  est <- estimate_mc(d$val)
  known <- simex(y ~ mc(z, est$Pi) + x1, data = d$df, method = "standard",
                 B = 20)
  fit <- simex(y ~ mc(z, est) + x1, data = d$df, method = "standard", B = 20,
               mc_variance = "conditional")
  expect_identical(coef(fit), coef(known))
  expect_identical(vcov(fit), vcov(known))
  expect_identical(fit$mc.variance, "conditional")
  expect_s3_class(fit$mc.estimate$z, "mc_estimate")
  # the improved method takes the estimate's prevalence
  imp <- simex(y ~ mc(z, est) + x1, data = d$df, mc_variance = "conditional")
  expect_equal(imp$pi.vec, unname(est$pi))
})

test_that("draws propagate the uncertainty of the estimate", {
  d <- .se_data()
  est <- estimate_mc(d$val)
  f1 <- simex(y ~ mc(z, est) + x1, data = d$df, mc_draws = 20)
  f2 <- simex(y ~ mc(z, est) + x1, data = d$df, mc_draws = 20)
  expect_identical(vcov(f1), vcov(f2))
  expect_gt(vcov(f1)[1, 1], f1$vcov.conditional[1, 1])
  expect_equal(nrow(f1$mc.draws$coefficients), 20L)
  expect_gt(sd(f1$mc.draws$coefficients[, 1]), 0)
  expect_match(f1$vcov.assumption, "Rubin's rule")
  expect_error(simex(y ~ mc(z, est) + x1, data = d$df, mc_draws = 1),
               "at least 2")
})

test_that("labels, data order and main-study proxies are checked", {
  d <- .se_data()
  lev <- c("no", "yes")
  est <- estimate_mc(validation_sample(lev[d$zv + 1], lev[d$zhv + 1]))
  df <- d$df
  df$z <- factor(lev[as.integer(as.character(d$df$z)) + 1], levels = rev(lev))
  fit <- simex(y ~ mc(z, est) + x1, data = df, method = "standard", B = 20,
               mc_variance = "conditional")
  df2 <- df
  df2$z <- factor(as.character(df$z), levels = lev)
  fit2 <- simex(y ~ mc(z, est) + x1, data = df2, method = "standard", B = 20,
                mc_variance = "conditional")
  expect_identical(coef(fit), coef(fit2))
  df3 <- df
  df3$z <- factor(c("no", "maybe")[as.integer(df$z == "yes") + 1])
  expect_error(simex(y ~ mc(z, est) + x1, data = df3), "'maybe'")
  idx <- sample(nrow(d$df), 200)
  est_int <- estimate_mc(validation_sample(d$z[idx], index = idx),
                         z_hat = as.integer(as.character(d$df$z)))
  other <- d$df
  other$z <- factor(rev(as.character(other$z)), levels = c("0", "1"))
  expect_error(simex(y ~ mc(z, est_int) + x1, data = other),
               "other main-study proxies")
})

test_that("matrices without valid powers are repaired with a warning", {
  bad <- matrix(c(0.60, 0.30, 0.10, 0.40, 0.30, 0.30, 0.05, 0.05, 0.90), 3L)
  expect_false(check.mc.matrix(list(bad)))
  expect_warning(fixed <- mismeasured:::.simex_valid_power(bad, "z"),
                 "build.mc.matrix")
  expect_true(check.mc.matrix(list(fixed)))
})

test_that("response misclassification accepts an estimate", {
  set.seed(1610)
  n <- 1000
  x1 <- rnorm(n)
  ytrue <- rbinom(n, 1, plogis(-0.3 + 0.8 * x1))
  yobs <- ifelse(ytrue == 1, rbinom(n, 1, 0.9), rbinom(n, 1, 0.1))
  yv <- rbinom(200, 1, 0.45)
  yvh <- ifelse(yv == 1, rbinom(200, 1, 0.9), rbinom(200, 1, 0.1))
  est <- estimate_mc(validation_sample(yv, yvh))
  df <- data.frame(y = factor(yobs), x1 = x1)
  fit <- simex(mc(y, est) ~ x1, data = df, family = binomial(), B = 20,
               mc_draws = 10)
  expect_true(all(is.finite(sqrt(diag(vcov(fit))))))
  expect_identical(names(fit$mc.estimate), "y")
})

test_that("prevalence from proxy frequencies uses EM", {
  Pi <- matrix(c(0.9, 0.1, 0.2, 0.8), 2L)
  z_hat <- rep(0:1, round(1000 * drop(Pi %*% c(0.6, 0.4))))
  expect_equal(mismeasured:::.estimate_pi_vec(z_hat, Pi), c(0.6, 0.4),
               tolerance = 1e-8)
  expect_warning(mismeasured:::.estimate_pi_vec(rep(0:1, c(95, 5)), Pi),
                 "boundary")
})

test_that("propagated simex standard errors are calibrated", {
  skip_on_cran()
  one <- function(r) {
    d <- .se_data(seed = 8000L + r, n = 2000L, nv = 200L)
    est <- estimate_mc(d$val, control = control_mc(on_ill = "none"))
    f <- suppressWarnings(simex(y ~ mc(z, est) + x1, data = d$df, mc_draws = 30))
    c(coef(f)[1], sqrt(vcov(f)[1, 1]), sqrt(f$vcov.conditional[1, 1]))
  }
  R <- t(vapply(seq_len(100L), one, numeric(3)))
  expect_lt(abs(mean(R[, 1]) - 0.8), 0.04)
  ratio <- mean(R[, 2]) / sd(R[, 1])
  expect_true(ratio > 0.8 && ratio < 1.25, label = paste("draws", ratio))
  expect_lt(mean(R[, 3]) / sd(R[, 1]), ratio)
})
