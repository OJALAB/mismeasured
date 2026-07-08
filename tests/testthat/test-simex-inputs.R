# Input validation: silent-misbehavior guards (improvement plan Phase 2).

.make_mc_df <- function(n = 500, seed = 1) {
  set.seed(seed)
  z <- rbinom(n, 1, 0.4)
  y <- 1 + 2 * z + rnorm(n)
  z_star <- ifelse(z == 1, 1 - rbinom(n, 1, 0.12), rbinom(n, 1, 0.1))
  data.frame(y = y, z = factor(z_star))
}
test_that("offset() terms are rejected loudly (plan 2.1)", {
  Pi <- matrix(c(0.9, 0.1, 0.12, 0.88), 2, 2)
  df <- .make_mc_df()
  df$expo <- runif(nrow(df), 1, 2)

  expect_error(
    simex(y ~ mc(z, Pi) + offset(log(expo)), data = df),
    "offset"
  )
  expect_error(
    mcglm(y ~ mc(z, Pi) + offset(log(expo)), data = df,
          family = gaussian()),
    "offset"
  )
})
