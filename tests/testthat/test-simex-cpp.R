# --------------------------------------------------------------------------
# Compiled SIMEX / MC-SIMEX backend guardrails.
# --------------------------------------------------------------------------

test_that("standard MC-SIMEX C++ rejects invalid transition matrices loudly", {
  y <- c(1, 2, 1, 2)
  z_hat <- c(0L, 1L, 0L, 1L)
  x <- matrix(1, nrow = 4, ncol = 1)
  wt <- rep(1, 4)

  expect_error(
    mismeasured:::mcsimex_sim_cpp(
      y, z_hat, x, matrix(c(0.9, 0.1, 0.2, 0.7), 2, 2),
      2L, 1L, 1, 1L, wt, 1L
    ),
    "columns must sum"
  )
  expect_error(
    mismeasured:::mcsimex_sim_cpp(
      y, z_hat, x, matrix(c(0.9, 0.1, -0.2, 1.2), 2, 2),
      2L, 1L, 1, 1L, wt, 1L
    ),
    "negative"
  )
  expect_error(
    mismeasured:::mcsimex_sim_cpp(
      y, z_hat, x, matrix(c(0.9, 0.1, NaN, 0.9), 2, 2),
      2L, 1L, 1, 1L, wt, 1L
    ),
    "non-finite"
  )
  expect_error(
    mismeasured:::mcsimex_sim_cpp(y, z_hat, x, diag(2), 2L, 99L, 1, 1L, wt, 1L),
    "Unknown dist_code"
  )
})
