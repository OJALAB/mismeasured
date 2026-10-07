# validation_sample() class, list coercion and binding to the main study.

test_that("validation_sample() records the design", {
  set.seed(1)
  z <- rbinom(60, 1, 0.4)
  zh <- ifelse(z == 1, rbinom(60, 1, 0.85), rbinom(60, 1, 0.1))
  v <- validation_sample(z, zh)
  expect_s3_class(v, "mc_validation")
  expect_identical(v$type, "external")
  expect_identical(v$estimator, "hajek")
  expect_output(print(v), "external.*simple random sample.*Hajek")

  vi <- validation_sample(z, index = 1:60, weights = rep(3, 60),
                          strata = "z_hat", estimator = "ht")
  expect_identical(vi$type, "internal")
  expect_identical(vi$strata_type, "z_hat")
  expect_output(print(vi), "design weights, stratified by z_hat")
})

test_that("validation_sample() rejects malformed input", {
  z <- c(0, 1, 1, 0)
  expect_error(validation_sample(c(NA, 1, 1, 0), z), "codes")
  expect_error(validation_sample(c(0.5, 1, 1, 0), z), "whole numbers")
  expect_error(validation_sample(list(0, 1), z), "labels")
  expect_error(validation_sample(z), "z_hat \\(external\\) or index")
  expect_error(validation_sample(z, index = c(1, 1, 2, 3)), "distinct")
  expect_error(validation_sample(z, z, weights = c(1, 1, 0, 1)), "positive")
  expect_error(validation_sample(z, z, strata = 1:3), "strata")
  expect_error(validation_sample(z, z, estimator = "ht"), "frame size N")
  expect_error(validation_sample(z, z, N_strata = c(10, 20)), "named")
})

test_that("as_validation_sample() coerces the list interface", {
  v <- as_validation_sample(list(z = c(0, 1, 1), z_hat = c(0, 1, 0)))
  expect_s3_class(v, "mc_validation")
  expect_identical(as_validation_sample(v), v)
  expect_error(as_validation_sample(list(z = 1, foo = 2)), "must be a list")
})

test_that("binding resolves internal proxies and design weights", {
  bind <- mismeasured:::.mc_bind_validation
  z_hat <- c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
  # internal SRS: constant weight n / m
  b <- bind(validation_sample(c(0, 1, 1, 0), index = c(1, 5, 6, 2)),
            z_hat, K = 2L)
  expect_equal(b$proxy, c(0L, 1L, 1L, 0L))
  expect_equal(b$d, rep(10 / 4, 4))
  expect_equal(b$N, 10)
  # internal, stratified by z_hat: N_j / n_j from main-study proxy counts
  b <- bind(validation_sample(c(0, 1, 1, 0), index = c(1, 5, 6, 2),
                              strata = "z_hat"), z_hat, K = 2L)
  expect_equal(b$d, c(4 / 2, 6 / 2, 6 / 2, 4 / 2))
  # external with user strata needs N_strata
  v <- validation_sample(c(0, 1, 1, 0), c(0, 1, 0, 0), strata = c("a", "a", "b", "b"))
  expect_error(bind(v, NULL, 2L), "N_strata")
  v$N_strata <- c(a = 50, b = 10)
  expect_equal(bind(v, NULL, 2L)$d, c(25, 25, 5, 5))
  # proxy mismatch, empty true category, codes beyond K
  expect_error(bind(validation_sample(c(0, 1), z_hat = c(1, 1), index = 1:2),
                    z_hat, 2L), "does not match")
  expect_error(bind(validation_sample(c(0, 0), z_hat = c(0, 1)), NULL, 2L),
               "Every true category")
  expect_error(bind(validation_sample(c(0, 2), z_hat = c(0, 1)), NULL, 2L),
               "'2' that are not categories")
})

test_that("mcglm() accepts a validation_sample() identical to the list form", {
  set.seed(5)
  n <- 400
  z <- rbinom(n, 1, 0.4)
  zh <- ifelse(z == 1, rbinom(n, 1, 0.85), rbinom(n, 1, 0.1))
  x <- cbind(1, rnorm(n))
  y <- rpois(n, exp(-0.3 + 0.8 * z + 0.5 * x[, 2]))
  idx <- sample.int(n, 120)
  f_list <- mcglm(y, z_hat = zh, x = x, method = "sub",
                  validation = list(z = z[idx], index = idx))
  f_obj <- mcglm(y, z_hat = zh, x = x, method = "sub",
                 validation = validation_sample(z[idx], index = idx))
  expect_identical(coef(f_obj, method = "sub"), coef(f_list, method = "sub"))
  expect_identical(vcov(f_obj, method = "sub"), vcov(f_list, method = "sub"))
})

test_that("audit categories are matched by label, not by position", {
  set.seed(9)
  n <- 600
  lev <- c("clerk", "manager", "sales")
  Pi <- matrix(c(0.85, 0.10, 0.05, 0.10, 0.80, 0.10, 0.05, 0.10, 0.85), 3L)
  z <- sample.int(3L, n, TRUE, c(0.5, 0.3, 0.2)) - 1L
  zh <- vapply(z, function(k) sample.int(3L, 1L, prob = Pi[, k + 1L]) - 1L, 0L)
  x1 <- rnorm(n)
  y <- rpois(n, exp(c(0, 0.6, -0.4)[z + 1L] - 0.3 + 0.5 * x1))
  df <- data.frame(y = y, occ = factor(lev[zh + 1L], levels = lev), x1 = x1)
  idx <- sample.int(n, 200L)

  ref <- mcglm(y, z_hat = zh, x = cbind(1, x1), method = "sub",
               validation = list(z = z[idx], index = idx))
  # character labels; and a factor whose levels are in another order
  f_chr <- mcglm(y ~ mc(occ) + x1, data = df, method = "sub",
                 validation = validation_sample(lev[z[idx] + 1L],
                                                index = idx))
  f_rev <- mcglm(y ~ mc(occ) + x1, data = df, method = "sub",
                 validation = validation_sample(
                   factor(lev[z[idx] + 1L], levels = rev(lev)), index = idx))
  expect_equal(unname(coef(f_chr, method = "sub")),
               unname(coef(ref, method = "sub")))
  expect_equal(coef(f_rev, method = "sub"), coef(f_chr, method = "sub"))
  expect_identical(rownames(f_chr$nuisance$sub$Pi), lev)
  expect_identical(names(f_chr$nuisance$sub$pi_z), lev)

  # numeric codes against labelled categories: error with a hint
  expect_error(mcglm(y ~ mc(occ) + x1, data = df, method = "sub",
                     validation = validation_sample(z[idx], index = idx)),
               "'0'.*not categories of the model \\(clerk, manager, sales\\).*levels")
  # an unknown label
  bad <- lev[z[idx] + 1L]
  bad[1] <- "cleark"
  expect_error(mcglm(y ~ mc(occ) + x1, data = df, method = "sub",
                     validation = validation_sample(bad, index = idx)),
               "'cleark'")
  # ISCO-like numeric labels 1..3: numbers match the labels, no shift
  df_num <- df
  df_num$occ <- factor(zh + 1L, levels = 1:3)
  f_num <- mcglm(y ~ mc(occ) + x1, data = df_num, method = "sub",
                 validation = validation_sample(z[idx] + 1L, index = idx))
  expect_equal(unname(coef(f_num, method = "sub")),
               unname(coef(ref, method = "sub")))
})

test_that("estimate_mc() categories must match the model", {
  set.seed(10)
  lev <- c("a", "b")
  zv <- rbinom(300, 1, 0.4)
  zhv <- ifelse(zv == 1, rbinom(300, 1, 0.85), rbinom(300, 1, 0.1))
  n <- 500
  z <- rbinom(n, 1, 0.4)
  zh <- ifelse(z == 1, rbinom(n, 1, 0.85), rbinom(n, 1, 0.1))
  df <- data.frame(y = rpois(n, exp(0.5 * z)),
                   occ = factor(lev[zh + 1L], levels = rev(lev)))
  est <- estimate_mc(validation_sample(lev[zv + 1L], lev[zhv + 1L]))
  expect_identical(est$levels, c("a", "b"))
  expect_error(mcglm(y ~ mc(occ), data = df, method = "sub",
                     validation = est), "differ from the model's \\(b, a\\)")
  est2 <- estimate_mc(validation_sample(lev[zv + 1L], lev[zhv + 1L]),
                      levels = levels(df$occ))
  expect_identical(rownames(est2$Pi), c("b", "a"))
  fit <- mcglm(y ~ mc(occ), data = df, method = "sub", validation = est2)
  expect_equal(unname(est2$Pi), unname(est$Pi[2:1, 2:1]))
  expect_true(all(is.finite(fit$se$sub)))
  # main-study proxies as a factor fix the labels too
  est3 <- estimate_mc(validation_sample(lev[zv + 1L], lev[zhv + 1L]),
                      z_hat = df$occ, control = control_mc(prevalence = "em"))
  expect_identical(est3$levels, c("b", "a"))
})
