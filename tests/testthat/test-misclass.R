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
