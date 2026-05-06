test_that("misclass() produces misclassified data", {
  set.seed(42)
  x1 <- factor(rbinom(200, 1, 0.5))
  p1 <- matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2)
  colnames(p1) <- levels(x1)
  x <- data.frame(x1 = x1)
  x.mc <- misclass(x, list(x1 = p1), k = 1)

  expect_s3_class(x.mc, "data.frame")
  expect_equal(ncol(x.mc), 1)
  expect_equal(levels(x.mc$x1), levels(x1))
  # Should not be identical (misclassification occurred)
  expect_false(identical(x$x1, x.mc$x1))
})

test_that("check.mc.matrix identifies valid/invalid matrices", {
  P1 <- matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2)
  P2 <- matrix(c(0.4, 0.6, 0.6, 0.4), nrow = 2)
  result <- check.mc.matrix(list(P1, P2))
  expect_equal(result, c(TRUE, FALSE))
})

test_that("build.mc.matrix produces valid matrix", {
  Pi <- matrix(c(0.989, 0.01, 0.001, 0.17, 0.829, 0.001,
                 0.001, 0.18, 0.819), nrow = 3, byrow = FALSE)
  result <- build.mc.matrix(Pi)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(3, 3))
  # Column sums should be ~1
  expect_true(all(abs(colSums(result) - 1) < 0.01))
})

# --- Additional coverage tests ---

test_that("misclass() input validation", {
  expect_error(misclass(data.frame(x = 1:3), "not_a_list"),
               "must be a list")
  expect_error(misclass(matrix(1, 3, 3), list(x = matrix(1, 2, 2))),
               "must be a data.frame")
  P <- matrix(c(0.9, 0.1, 0.2, 0.8), 2, 2)
  expect_error(
    misclass(data.frame(a = factor(0:1)),
             list(b = P)),
    "do not match"
  )
  expect_error(misclass(data.frame(a = factor(0:1)),
                        list(a = P), k = -1),
               "non-negative")
})

test_that("misclass() handles K=3 factors", {
  set.seed(123)
  x <- factor(sample(letters[1:3], 300, replace = TRUE))
  P <- matrix(0.05, 3, 3); diag(P) <- 0.9
  colnames(P) <- rownames(P) <- levels(x)
  out <- misclass(data.frame(x = x), list(x = P), k = 1)
  expect_equal(levels(out$x), levels(x))
  # With high diagonal, most observations should remain in the same class.
  expect_gt(mean(out$x == x), 0.7)
})

test_that("build.mc.matrix supports method='log' and method='jlt'", {
  Pi <- matrix(c(0.989, 0.01, 0.001, 0.17, 0.829, 0.001,
                 0.001, 0.18, 0.819), nrow = 3, byrow = FALSE)
  res_log <- build.mc.matrix(Pi, method = "log")
  expect_true(is.matrix(res_log))
  expect_equal(dim(res_log), c(3, 3))

  res_jlt <- build.mc.matrix(Pi, method = "jlt")
  expect_true(is.matrix(res_jlt))
  expect_equal(dim(res_jlt), c(3, 3))
})

test_that("build.mc.matrix(method='jlt') errors when diagonal contains 0", {
  Pi_bad <- matrix(c(0.0, 1.0, 0.5, 0.5), nrow = 2, byrow = FALSE)
  expect_error(build.mc.matrix(Pi_bad, method = "jlt"),
               "0 on the diagonal")
})

test_that("build.mc.matrix with diag.cor=TRUE produces valid output", {
  Pi <- matrix(c(0.989, 0.01, 0.001, 0.17, 0.829, 0.001,
                 0.001, 0.18, 0.819), nrow = 3, byrow = FALSE)
  result <- build.mc.matrix(Pi, diag.cor = TRUE)
  expect_equal(dim(result), c(3, 3))
  expect_true(all(abs(colSums(result) - 1) < 0.01))
})

test_that("check.mc.matrix handles K=3 matrices", {
  Pi_valid <- matrix(c(0.989, 0.01, 0.001, 0.17, 0.829, 0.001,
                       0.001, 0.18, 0.819), nrow = 3, byrow = FALSE)
  Pi_fixed <- build.mc.matrix(Pi_valid)
  res <- check.mc.matrix(list(Pi_valid, Pi_fixed))
  expect_length(res, 2L)
  expect_type(res, "logical")
})

# --- diag.block ---

test_that("diag.block(list) builds block-diagonal matrix from a list", {
  a <- matrix(1, 2, 2)
  b <- matrix(2, 2, 3)
  out <- diag.block(list(a, b))
  expect_equal(dim(out), c(4, 5))
  expect_equal(out[1:2, 1:2], a)
  expect_equal(out[3:4, 3:5], b)
  expect_true(all(out[1:2, 3:5] == 0))
  expect_true(all(out[3:4, 1:2] == 0))
})

test_that("diag.block(matrix, n) replicates via Kronecker", {
  a <- matrix(c(1, 2, 3, 4), 2, 2)
  out <- diag.block(a, 3L)
  expect_equal(dim(out), c(6, 6))
  for (i in 0:2) {
    expect_equal(out[(2 * i + 1):(2 * i + 2), (2 * i + 1):(2 * i + 2)], a)
  }
})

test_that("diag.block accepts mixed list of matrices and vectors", {
  out <- diag.block(list(matrix(1, 2, 2), c(1, 2, 3)))
  # NROW(c(1,2,3)) = 3, NCOL = 1
  expect_equal(dim(out), c(5, 3))
})
