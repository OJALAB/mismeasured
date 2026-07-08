# SE calibration for the improved MC-SIMEX (improvement plan 1.4).
# The mean reported SE must track the Monte-Carlo SD of the estimator;
# a ratio outside [0.8, 1.25] means the variance formula is mis-calibrated
# (e.g. SEs shrinking with B, or using the naive vcov at lambda = 0).

.simex_se_calibration <- function(n_rep, gen, B = 20) {
  res <- t(vapply(seq_len(n_rep), function(r) {
    df <- gen(r)
    fit <- simex(y ~ mc(z, attr(df, "Pi")), data = df,
                 family = attr(df, "family"), B = B, seed = r)
    c(coef(fit), sqrt(diag(fit$vcov)))
  }, numeric(4)))
  p <- 2L
  colMeans(res[, p + seq_len(p)]) / apply(res[, seq_len(p)], 2, sd)
}

test_that("improved MC-SIMEX SEs are calibrated (gaussian)", {
  skip_on_cran()

  Pi <- matrix(c(0.9, 0.1, 0.12, 0.88), 2, 2)
  gen <- function(r) {
    set.seed(r)
    n <- 1000
    z <- rbinom(n, 1, 0.4)
    y <- 1 + 2 * z + rnorm(n)
    z_star <- ifelse(z == 1, 1 - rbinom(n, 1, 0.12), rbinom(n, 1, 0.1))
    df <- data.frame(y = y, z = factor(z_star))
    attr(df, "Pi") <- Pi
    attr(df, "family") <- gaussian()
    df
  }

  ratios <- .simex_se_calibration(200, gen)
  expect_true(all(ratios > 0.8 & ratios < 1.25),
              info = paste("SE/MC-SD ratios:",
                           paste(round(ratios, 3), collapse = ", ")))
})

test_that("improved MC-SIMEX SEs are calibrated (binomial)", {
  skip_on_cran()

  Pi <- matrix(c(0.9, 0.1, 0.12, 0.88), 2, 2)
  gen <- function(r) {
    set.seed(r)
    n <- 1000
    z <- rbinom(n, 1, 0.4)
    y <- rbinom(n, 1, plogis(0.5 + 1 * z))
    z_star <- ifelse(z == 1, 1 - rbinom(n, 1, 0.12), rbinom(n, 1, 0.1))
    df <- data.frame(y = y, z = factor(z_star))
    attr(df, "Pi") <- Pi
    attr(df, "family") <- binomial()
    df
  }

  ratios <- .simex_se_calibration(200, gen)
  expect_true(all(ratios > 0.8 & ratios < 1.25),
              info = paste("SE/MC-SD ratios:",
                           paste(round(ratios, 3), collapse = ", ")))
})
