# Input validation: guards against silent misbehavior on valid inputs.

.make_mc_df <- function(n = 500, seed = 1) {
  set.seed(seed)
  z <- rbinom(n, 1, 0.4)
  y <- 1 + 2 * z + rnorm(n)
  z_star <- ifelse(z == 1, 1 - rbinom(n, 1, 0.12), rbinom(n, 1, 0.1))
  data.frame(y = y, z = factor(z_star))
}
test_that("offset() terms are rejected", {
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

test_that("grouped binomial cbind(succ, fail) uses trial counts", {
  Pi <- matrix(c(0.9, 0.1, 0.12, 0.88), 2, 2)
  set.seed(3)
  n_grp <- 60
  trials <- sample(5:30, n_grp, replace = TRUE)
  z <- rbinom(n_grp, 1, 0.4)
  succ <- rbinom(n_grp, trials, plogis(0.5 + z))
  z_star <- ifelse(z == 1, 1 - rbinom(n_grp, 1, 0.12), rbinom(n_grp, 1, 0.1))

  grouped <- data.frame(succ = succ, fail = trials - succ,
                        z = factor(z_star))
  expanded <- data.frame(
    y = unlist(mapply(function(s, t) rep(c(1, 0), c(s, t - s)), succ, trials,
                      SIMPLIFY = FALSE)),
    z = factor(rep(z_star, trials))
  )

  # lambda = 0 makes the improved fit deterministic (Pi^0 = I), so the
  # grouped and expanded representations must agree exactly.
  fit_g <- simex(cbind(succ, fail) ~ mc(z, Pi), data = grouped,
                 family = binomial(), method = "improved",
                 lambda = 0, B = 1, jackknife = FALSE)
  fit_e <- simex(y ~ mc(z, Pi), data = expanded,
                 family = binomial(), method = "improved",
                 lambda = 0, B = 1, jackknife = FALSE)
  expect_equal(coef(fit_g), coef(fit_e), tolerance = 1e-6)

  # General path (resampling + variance) runs and is sane.
  fit <- simex(cbind(succ, fail) ~ mc(z, Pi), data = grouped,
               family = binomial(), method = "improved", B = 20, seed = 42)
  expect_true(all(is.finite(coef(fit))))
  expect_true(all(is.finite(fit$vcov)))
  # Ignoring trial counts (weight 1 per group) would inflate SEs by
  # roughly sqrt(mean(trials)); pin against the expanded-data fit.
  fit_e2 <- simex(y ~ mc(z, Pi), data = expanded,
                  family = binomial(), method = "improved", B = 20,
                  seed = 42)
  expect_equal(sqrt(diag(fit$vcov)), sqrt(diag(fit_e2$vcov)),
               tolerance = 0.25)
})

test_that("mcglm matrix interface validates z_hat coding", {
  set.seed(5)
  n <- 200
  z <- rbinom(n, 1, 0.4)
  x <- cbind(1, rnorm(n))
  y <- 1 + 2 * z + rnorm(n)
  Pi <- matrix(c(0.9, 0.1, 0.12, 0.88), 2, 2)

  # 1-based coding must be rejected, not silently distorted
  expect_error(
    mcglm(y, z_hat = z + 1L, x = x, family = gaussian(), Pi = Pi,
          method = "cs"),
    "0, ..., K-1", fixed = TRUE
  )
  expect_error(
    mcglm(y, z_hat = factor(z), x = x, family = gaussian(), Pi = Pi,
          method = "cs"),
    "factor"
  )
  expect_error(
    mcglm(y, z_hat = z + 0.5, x = x, family = gaussian(), Pi = Pi,
          method = "cs"),
    "integer"
  )

  # K comes from Pi when a category is unobserved
  Pi3 <- matrix(c(0.9, 0.05, 0.05,
                  0.05, 0.9, 0.05,
                  0.05, 0.05, 0.9), 3, 3)
  z3 <- ifelse(z == 1, 2L, 0L)  # category 1 unobserved
  # The all-zero dummy makes the naive vcov singular; mcglm warns about
  # that (correctly) -- the pin here is that K comes from Pi, not from
  # length(unique(z_hat)) which would give 2 and misindex categories.
  suppressWarnings(
    fit <- mcglm(y, z_hat = z3, x = x, family = gaussian(), Pi = Pi3,
                 method = "naive")
  )
  expect_identical(fit$K, 3L)
})

test_that("non-canonical links are rejected", {
  Pi <- matrix(c(0.9, 0.1, 0.12, 0.88), 2, 2)
  df <- .make_mc_df()
  df$yb <- rbinom(nrow(df), 1, 0.5)

  expect_error(
    simex(yb ~ mc(z, Pi), data = df, family = binomial("probit")),
    "canonical"
  )
  expect_error(
    mcglm(yb ~ mc(z, Pi), data = df, family = binomial("probit"),
          method = "cs"),
    "canonical"
  )
  expect_error(
    simex(y ~ mc(z, Pi), data = df, family = gaussian("log")),
    "canonical"
  )

  # canonical links keep working
  fit <- simex(yb ~ mc(z, Pi), data = df, family = binomial(), B = 5,
               seed = 1)
  expect_true(all(is.finite(coef(fit))))
})

test_that("rows with NAs are dropped consistently up front", {
  Pi <- matrix(c(0.9, 0.1, 0.12, 0.88), 2, 2)
  df <- .make_mc_df(n = 400)
  df$x <- rnorm(nrow(df))
  df_na <- df
  df_na$x[c(5, 50)] <- NA
  df_na$y[100] <- NA
  df_na$z[200] <- NA

  expect_warning(
    fit_na <- simex(y ~ mc(z, Pi) + x, data = df_na, B = 5, seed = 9),
    "dropped"
  )
  fit_cc <- simex(y ~ mc(z, Pi) + x,
                  data = df_na[complete.cases(df_na), ], B = 5, seed = 9)
  expect_equal(coef(fit_na), coef(fit_cc))
  expect_identical(fit_na$n, nrow(df_na) - 4L)

  # me() path too
  df_me <- df
  df_me$xs <- df_me$x + rnorm(nrow(df_me), sd = 0.3)
  df_me$xs[10] <- NA
  expect_warning(
    fit_me <- simex(y ~ me(xs, 0.3), data = df_me, B = 10, seed = 9),
    "dropped"
  )
  expect_true(all(is.finite(coef(fit_me))))

  # mismatched weights after NA drop must error, not recycle
  expect_error(
    suppressWarnings(simex(y ~ mc(z, Pi) + x, data = df_na, B = 5,
                           weights = rep(1, 10))),
    "weights"
  )
})

test_that("me()/mc() variables in non-main-effect terms are rejected", {
  Pi <- matrix(c(0.9, 0.1, 0.12, 0.88), 2, 2)
  df <- .make_mc_df(n = 300)
  df$x <- rnorm(nrow(df))
  df$xs <- df$x + rnorm(nrow(df), sd = 0.3)

  expect_error(
    simex(y ~ mc(z, Pi) + x:z, data = df),
    "bare main effect"
  )
  expect_error(
    simex(y ~ me(xs, 0.3) + xs:x, data = df),
    "bare main effect"
  )
  expect_error(
    simex(y ~ me(xs, 0.3) + I(xs^2), data = df),
    "bare main effect"
  )
  expect_error(
    simex(y ~ me(xs, 0.3) + log(abs(xs)), data = df),
    "bare main effect"
  )

  # interactions among clean covariates stay allowed
  df$w <- rnorm(nrow(df))
  fit <- simex(y ~ mc(z, Pi) + x * w, data = df, B = 5, seed = 1)
  expect_true(all(is.finite(coef(fit))))
})

test_that("default method is improved for gaussian, standard otherwise", {
  Pi <- matrix(c(0.9, 0.1, 0.12, 0.88), 2, 2)
  df <- .make_mc_df(n = 300)
  df$yb <- rbinom(nrow(df), 1, plogis(df$y - 1))

  expect_identical(simex(y ~ mc(z, Pi), data = df, B = 1)$method,
                   "improved")
  expect_identical(
    simex(yb ~ mc(z, Pi), data = df, family = binomial(), B = 10,
          jackknife = FALSE)$method,
    "standard"
  )
})
