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
  expect_error(validation_sample(factor(z), z), "codes")
  expect_error(validation_sample(c(NA, 1, 1, 0), z), "codes")
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
               "codes")
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
