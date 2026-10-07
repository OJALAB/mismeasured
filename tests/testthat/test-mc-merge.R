# suggest_merge() and merge_levels().

.mg_data <- function(seed = 1501L, nv = 600L) {
  set.seed(seed)
  K <- 4L
  Pi <- matrix(c(0.90, 0.04, 0.03, 0.03,
                 0.04, 0.47, 0.45, 0.04,
                 0.04, 0.45, 0.47, 0.04,
                 0.03, 0.03, 0.04, 0.90), K, K)
  lev <- c("a", "b", "c", "d")
  draw <- function(m) {
    z <- sample.int(K, m, TRUE) - 1L
    zh <- vapply(z, function(k) sample.int(K, 1L, prob = Pi[, k + 1L]) - 1L, 0L)
    list(z = lev[z + 1L], z_hat = lev[zh + 1L], code = z)
  }
  list(draw = draw, Pi = Pi, lev = lev, v = draw(nv))
}

test_that("collapsing averages columns with the prevalences and sums rows", {
  Pi <- matrix(c(0.8, 0.1, 0.1, 0.2, 0.7, 0.1, 0.1, 0.3, 0.6), 3L)
  pi <- c(0.5, 0.3, 0.2)
  col <- mismeasured:::.mc_collapse(Pi, pi, list(1L, 2:3))
  expect_equal(col$pi, c(0.5, 0.5))
  expect_equal(col$Pi[, 1], c(0.8, 0.2))
  expect_equal(col$Pi[2, 2], (0.3 * 0.8 + 0.2 * 0.9) / 0.5)
  expect_equal(colSums(col$Pi), c(1, 1))
})

test_that("the planted confused pair is merged first", {
  d <- .mg_data()
  est <- estimate_mc(validation_sample(d$v$z, d$v$z_hat),
                     control = control_mc(on_ill = "none"))
  sm <- suggest_merge(est)
  expect_s3_class(sm, "mc_merge")
  expect_false(sm$steps$pass[1])
  expect_identical(sm$steps$merged[2], "b & c")
  expect_identical(sm$steps$new_label[2], "b+c")
  expect_true(sm$steps$pass[2])
  expect_gt(sm$steps$sigma_min[2], sm$steps$sigma_min[1])
  expect_gt(sm$steps$confusion_share[2], 0.35)   # each goes to the other 45% of the time
  expect_identical(unname(sm$mapping), c("a", "b+c", "b+c", "d"))
  expect_gt(sm$steps$p_homogeneity[2], 0.05)
  expect_output(print(sm), "Suggestions only")
  sk <- suggest_merge(est, criterion = "kappa")
  expect_identical(sk$steps$merged[2], "b & c")
})

test_that("groups restrict the merges", {
  d <- .mg_data()
  est <- estimate_mc(validation_sample(d$v$z, d$v$z_hat),
                     control = control_mc(on_ill = "none"))
  sm <- suggest_merge(est, groups = c(1, 1, 2, 2), max_merges = 2)
  merged <- sm$steps$merged[-1]
  expect_false(any(grepl("b & c", merged, fixed = TRUE)))
  expect_true(all(merged %in% c("a & b", "c & d")))
  expect_error(suggest_merge(est, groups = 1:2), "group for each")
})

test_that("merge_levels recodes data and validation samples consistently", {
  f <- merge_levels(factor(c("a", "c", "b", "d"), levels = c("a", "b", "c", "d")),
                    c(b = "b+c", c = "b+c"))
  expect_identical(levels(f), c("a", "b+c", "d"))
  expect_identical(as.character(f), c("a", "b+c", "b+c", "d"))

  d <- .mg_data()
  est <- estimate_mc(validation_sample(d$v$z, d$v$z_hat),
                     control = control_mc(on_ill = "none"))
  sm <- suggest_merge(est)
  vs <- validation_sample(d$v$z, d$v$z_hat, strata = "z_hat",
                          N_strata = c(a = 100, b = 50, c = 70, d = 80))
  vm <- merge_levels(vs, sm)
  expect_identical(sort(unique(vm$z)), c("a", "b+c", "d"))
  expect_equal(vm$N_strata[["b+c"]], 120)
  # step 0 is the identity
  expect_identical(merge_levels(vs, sm, step = 0)$z, vs$z)

  # refit on merged categories, internal audit
  n <- 2000
  main <- d$draw(n)
  x1 <- rnorm(n)
  y <- rpois(n, exp(c(0, 0.5, 0.5, -0.4)[main$code + 1L] - 0.3 + 0.5 * x1))
  df <- data.frame(y = y, occ = factor(main$z_hat, levels = d$lev), x1 = x1)
  idx <- sample(n, 400)
  val <- validation_sample(main$z[idx], index = idx)
  est_main <- estimate_mc(val, z_hat = df$occ,
                          control = control_mc(on_ill = "none"))
  sm_main <- suggest_merge(est_main)
  df$occ_merged <- merge_levels(df$occ, sm_main)
  fit <- mcglm(y ~ mc(occ_merged) + x1, data = df, method = c("sub", "cs_akn"),
               validation = merge_levels(val, sm_main))
  expect_identical(fit$z_levels, c("a", "b+c", "d"))
  expect_true(all(is.finite(fit$se$cs_akn)))
  expect_error(merge_levels(df$occ, list(a = "x")), "named character")
})

test_that("indistinguishable columns are merged even when diagnostics pass", {
  d <- .mg_data(seed = 1510L, nv = 600L)
  est <- estimate_mc(validation_sample(d$v$z, d$v$z_hat),
                     control = control_mc(on_ill = "none"))
  sm <- suggest_merge(est)
  if (sm$steps$pass[1]) {
    # sampling noise separated the columns: the homogeneity test still merges
    expect_identical(sm$steps$merged[2], "b & c")
    expect_gt(sm$steps$p_homogeneity[2], 0.05)
  }
  # alpha = 1 stops as soon as the diagnostics pass
  expect_true(nrow(suggest_merge(est, alpha = 1)$steps) == 1L ||
                !sm$steps$pass[1])
  expect_equal(mismeasured:::.mc_homogeneity_p(c(50, 0), c(0, 50)),
               0, tolerance = 1e-10)
  expect_gt(mismeasured:::.mc_homogeneity_p(c(25, 25), c(24, 26)), 0.5)
})
