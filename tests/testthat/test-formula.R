test_that("me() creates me_term object", {
  term <- me(x, 0.5)
  expect_s3_class(term, "me_term")
  expect_equal(term$variable, "x")
})

test_that("mc() creates mc_term object", {
  Pi <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
  term <- mc(z, Pi)
  expect_s3_class(term, "mc_term")
  expect_equal(term$variable, "z")
})

test_that("parse_simex_formula strips me() and returns clean formula", {
  df <- data.frame(y = 1:10, x = rnorm(10), w = rnorm(10))
  parsed <- parse_simex_formula(y ~ me(x, 0.5) + w, df, globalenv())
  expect_equal(deparse(parsed$clean_formula), "y ~ x + w")
  expect_equal(parsed$error_type, "me")
  expect_length(parsed$me_terms, 1)
  expect_equal(parsed$me_terms[[1]]$variable, "x")
  expect_equal(parsed$me_terms[[1]]$sd, 0.5)
})

test_that("parse_simex_formula strips mc() and returns clean formula", {
  Pi <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
  df <- data.frame(y = 1:10, z = factor(rep(0:1, 5)), x = rnorm(10))
  parsed <- parse_simex_formula(y ~ mc(z, Pi) + x, df, environment())
  expect_equal(deparse(parsed$clean_formula), "y ~ z + x")
  expect_equal(parsed$error_type, "mc")
  expect_length(parsed$mc_terms, 1)
})

test_that("multiple me() terms collected correctly", {
  df <- data.frame(y = 1:10, x1 = rnorm(10), x2 = rnorm(10))
  parsed <- parse_simex_formula(y ~ me(x1, 0.3) + me(x2, 0.5), df, globalenv())
  expect_length(parsed$me_terms, 2)
  expect_equal(parsed$me_terms[[1]]$variable, "x1")
  expect_equal(parsed$me_terms[[2]]$variable, "x2")
})

test_that("no me/mc terms gives error", {
  df <- data.frame(y = 1:10, x = rnorm(10))
  expect_error(parse_simex_formula(y ~ x, df, globalenv()),
               "no me\\(\\) or mc\\(\\)")
})

test_that("me() on factor gives error in simex()", {
  df <- data.frame(y = 1:10, z = factor(rep(0:1, 5)))
  expect_error(simex(y ~ me(z, 0.5), data = df),
               "factor")
})

test_that("mc() with wrong K gives error", {
  Pi3 <- matrix(1/3, 3, 3)
  df <- data.frame(y = 1:10, z = factor(rep(0:1, 5)))
  expect_error(parse_simex_formula(y ~ mc(z, Pi3), df, environment()),
               "3 rows.*2 levels")
})

test_that("mc() with non-stochastic matrix gives error", {
  Pi_bad <- matrix(c(0.9, 0.3, 0.1, 0.7), 2, 2)
  df <- data.frame(y = 1:10, z = factor(rep(0:1, 5)))
  expect_error(parse_simex_formula(y ~ mc(z, Pi_bad), df, environment()),
               "sum to 1")
})

test_that("mc() coerces non-factor with warning", {
  Pi <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
  df <- data.frame(y = 1:10, z = rep(0:1, 5), x = rnorm(10))
  expect_warning(parse_simex_formula(y ~ mc(z, Pi) + x, df, environment()),
                 "not a factor")
})

test_that("heteroscedastic me() with column name works", {
  df <- data.frame(y = rnorm(20), x = rnorm(20), sd_x = abs(rnorm(20, 0.5, 0.1)))
  parsed <- parse_simex_formula(y ~ me(x, sd_x), df, globalenv())
  expect_length(parsed$me_terms[[1]]$sd, 20)
})

test_that("parse_simex_formula rejects unsupported or ambiguous error terms", {
  Pi <- diag(2)
  df <- data.frame(
    y = rnorm(8),
    x = rnorm(8),
    z = factor(rep(0:1, 4)),
    bad_sd = letters[1:8]
  )

  expect_error(
    parse_simex_formula(me(y, 0.1) ~ x, df, environment()),
    "Response measurement error"
  )
  expect_error(
    parse_simex_formula(y ~ me(x, 0.1):z, df, environment()),
    "cannot appear inside interactions"
  )
  expect_error(
    parse_simex_formula(y ~ I(mc(z, Pi)), df, environment()),
    "cannot appear inside I"
  )
  expect_error(
    parse_simex_formula(y ~ me(x, bad_sd), df, environment()),
    "must be numeric"
  )
  expect_error(
    parse_simex_formula(y ~ me(x, -1), df, environment()),
    "non-negative"
  )
  expect_error(
    parse_simex_formula(y ~ mc(z, matrix(1, nrow = 2, ncol = 3)), df,
                        environment()),
    "must be square"
  )
})

test_that("response mc() formula path validates factors and matrix columns", {
  Pi <- diag(2)
  df_num <- data.frame(y = rep(0:1, 4), x = rnorm(8))

  expect_warning(
    parsed <- parse_simex_formula(mc(y, Pi) ~ x, df_num, environment()),
    "not a factor"
  )
  expect_equal(parsed$error_type, "mc")
  expect_equal(parsed$response_mc$variable, "y")

  df_fac <- data.frame(y = factor(rep(0:1, 4)), x = rnorm(8))
  Pi_bad <- matrix(c(0.9, 0.1, 0.2, 0.7), 2, 2)
  expect_error(
    parse_simex_formula(mc(y, Pi_bad) ~ x, df_fac, environment()),
    "response matrix must sum to 1"
  )
})
