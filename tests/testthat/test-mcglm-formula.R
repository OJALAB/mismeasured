# --------------------------------------------------------------------------
# Tests for the mcglm() formula interface: y ~ mc(z, Pi) + ...
# --------------------------------------------------------------------------

# Helper: generate binary misclassification data
gen_binary <- function(n = 2000, seed = 42) {
  set.seed(seed)
  p01 <- 0.10; p10 <- 0.15
  Pi <- matrix(c(1 - p01, p01, p10, 1 - p10), 2, 2)
  z <- rbinom(n, 1, 0.4)
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  eta <- 0.8 * z + -0.5 + 0.7 * x1
  y <- rpois(n, exp(eta))
  df <- data.frame(y = y, z = factor(z_hat), x1 = x1, x2 = x2)
  list(df = df, Pi = Pi, p01 = p01, p10 = p10, psi0 = c(0.8, -0.5, 0.7),
       z_hat_int = as.integer(z_hat), n = n)
}

# Helper: generate multicategory misclassification data
gen_multi <- function(n = 2000, K = 3, seed = 55) {
  set.seed(seed)
  pi_z <- c(0.5, 0.3, 0.2)
  z <- sample(0:(K - 1), n, replace = TRUE, prob = pi_z)
  Pi <- matrix(c(0.8, 0.1, 0.1,
                 0.1, 0.8, 0.1,
                 0.1, 0.1, 0.8), nrow = K, byrow = TRUE)
  z_hat <- integer(n)
  for (i in seq_len(n)) z_hat[i] <- sample(0:(K - 1), 1, prob = Pi[, z[i] + 1])
  x1 <- rnorm(n)
  gamma <- c(0, 1.0, -0.9)
  alpha <- c(0.8, -0.7)
  eta <- alpha[1] + alpha[2] * x1 + gamma[z + 1]
  y <- rpois(n, exp(eta))
  df <- data.frame(y = y, z = factor(z_hat), x1 = x1)
  list(df = df, Pi = Pi, pi_z = pi_z, K = K, n = n)
}


# ===== MINIMAL FORMULA CALL (only Pi, no pi_z) =====

test_that("minimal formula: only Pi, pi_z auto-estimated (binary)", {
  d <- gen_binary(); Pi <- d$Pi
  # No pi_z, no p01, no p10 — all derived from Pi + data
  fit <- mcglm(y ~ mc(z, Pi) + x1, data = d$df, family = "poisson")
  expect_s3_class(fit, "mcglm")
  expect_equal(fit$K, 2L)
  expect_equal(fit$n, d$n)
  expect_equal(fit$p, 3)
  expect_named(fit$coefficients, c("naive", "bca", "bcm", "cs"))
  for (nm in names(fit$coefficients)) {
    expect_true(all(is.finite(fit$coefficients[[nm]])),
                info = paste("method:", nm))
  }
})

test_that("minimal formula: only Pi, pi_z auto-estimated (multicategory)", {
  d <- gen_multi(); Pi <- d$Pi
  fit <- mcglm(y ~ mc(z, Pi) + x1, data = d$df, family = "poisson")
  expect_s3_class(fit, "mcglm")
  expect_equal(fit$K, d$K)
  expect_equal(fit$p, d$K - 1 + 2)
})


# ===== BASIC FORMULA INTERFACE =====

test_that("formula interface extracts Pi from mc() term", {
  d <- gen_binary(); Pi <- d$Pi
  fit <- mcglm(y ~ mc(z, Pi) + x1, data = d$df, family = "poisson")
  expect_s3_class(fit, "mcglm")
  expect_equal(fit$K, 2L)
})

test_that("formula interface auto-derives p01/p10/pi_z from Pi", {
  d <- gen_binary(); Pi <- d$Pi
  fit <- mcglm(y ~ mc(z, Pi) + x1, data = d$df, family = "poisson",
               method = c("naive", "bca"))
  expect_s3_class(fit, "mcglm")
  expect_true(all(is.finite(fit$coefficients$bca)))
})

test_that("explicit pi_z overrides auto-estimated value", {
  d <- gen_binary(); Pi <- d$Pi
  fit1 <- mcglm(y ~ mc(z, Pi) + x1, data = d$df, family = "poisson",
                 pi_z = 0.4, method = c("naive", "bca"))
  fit2 <- mcglm(y ~ mc(z, Pi) + x1, data = d$df, family = "poisson",
                 method = c("naive", "bca"))
  # Both should work; results may differ slightly

  expect_true(all(is.finite(fit1$coefficients$bca)))
  expect_true(all(is.finite(fit2$coefficients$bca)))
})


# ===== FORMULA WITH VARIOUS X TRANSFORMATIONS =====

test_that("formula with I(x^2) works", {
  d <- gen_binary(); Pi <- d$Pi
  fit <- mcglm(y ~ mc(z, Pi) + x1 + I(x1^2), data = d$df,
               family = "poisson")
  expect_s3_class(fit, "mcglm")
  expect_equal(fit$p, 4)  # gamma + intercept + x1 + I(x1^2)
})

test_that("formula with log() transformation works", {
  d <- gen_binary(); Pi <- d$Pi
  d$df$x_pos <- abs(d$df$x1) + 0.1
  fit <- mcglm(y ~ mc(z, Pi) + log(x_pos), data = d$df,
               family = "poisson", method = "naive")
  expect_s3_class(fit, "mcglm")
  expect_equal(fit$p, 3)
})

test_that("formula with multiple covariates works", {
  d <- gen_binary(); Pi <- d$Pi
  fit <- mcglm(y ~ mc(z, Pi) + x1 + x2, data = d$df,
               family = "poisson", method = "naive")
  expect_s3_class(fit, "mcglm")
  expect_equal(fit$p, 4)  # gamma + intercept + x1 + x2
})

test_that("formula with interaction x1:x2 works", {
  d <- gen_binary(); Pi <- d$Pi
  fit <- mcglm(y ~ mc(z, Pi) + x1 + x2 + x1:x2, data = d$df,
               family = "poisson", method = "naive")
  expect_s3_class(fit, "mcglm")
  expect_equal(fit$p, 5)  # gamma + intercept + x1 + x2 + x1:x2
})

test_that("formula with poly() works", {
  d <- gen_binary(); Pi <- d$Pi
  fit <- mcglm(y ~ mc(z, Pi) + poly(x1, 2), data = d$df,
               family = "poisson", method = "naive")
  expect_s3_class(fit, "mcglm")
  expect_equal(fit$p, 4)  # gamma + intercept + poly1 + poly2
})

test_that("formula with no extra covariates (intercept only)", {
  d <- gen_binary(); Pi <- d$Pi
  fit <- mcglm(y ~ mc(z, Pi), data = d$df,
               family = "poisson", method = "naive")
  expect_s3_class(fit, "mcglm")
  expect_equal(fit$p, 2)  # gamma + intercept
})


# ===== FORMULA vs MATRIX INTERFACE CONSISTENCY =====

test_that("formula and matrix interfaces give same naive results", {
  d <- gen_binary(); Pi <- d$Pi
  fit_formula <- mcglm(y ~ mc(z, Pi) + x1, data = d$df,
                        family = "poisson", method = "naive")

  x_mat <- cbind(1, d$df$x1)
  fit_matrix <- mcglm(d$df$y, z_hat = d$z_hat_int, x = x_mat,
                       family = "poisson", method = "naive")

  expect_equal(unname(fit_formula$coefficients$naive),
               unname(fit_matrix$coefficients$naive), tolerance = 1e-8)
})

test_that("formula and matrix give same corrected results (binary)", {
  d <- gen_binary(); Pi <- d$Pi
  fit_f <- mcglm(y ~ mc(z, Pi) + x1, data = d$df,
                  family = "poisson", pi_z = 0.4,
                  method = c("naive", "bca", "cs"))

  x_mat <- cbind(1, d$df$x1)
  fit_m <- mcglm(d$df$y, z_hat = d$z_hat_int, x = x_mat,
                  family = "poisson", p01 = d$p01, p10 = d$p10, pi_z = 0.4,
                  method = c("naive", "bca", "cs"))

  expect_equal(unname(fit_f$coefficients$bca),
               unname(fit_m$coefficients$bca), tolerance = 1e-8)
  expect_equal(unname(fit_f$coefficients$cs),
               unname(fit_m$coefficients$cs), tolerance = 1e-8)
})


# ===== BINARY MISCLASSIFICATION - ALL FAMILIES =====

test_that("formula: binary Poisson BCA/BCM/CS all finite", {
  d <- gen_binary(); Pi <- d$Pi
  fit <- mcglm(y ~ mc(z, Pi) + x1, data = d$df, family = "poisson")
  for (nm in names(fit$coefficients)) {
    expect_true(all(is.finite(fit$coefficients[[nm]])),
                info = paste("method:", nm))
  }
})

test_that("formula: binary Binomial works", {
  set.seed(99)
  n <- 1000; Pi <- matrix(c(0.92, 0.08, 0.12, 0.88), 2, 2)
  z <- rbinom(n, 1, 0.4)
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, 0.08)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.12)
  x1 <- rnorm(n)
  y <- rbinom(n, 1, plogis(0.8 * z - 0.5 + 0.7 * x1))
  df <- data.frame(y = y, z = factor(z_hat), x1 = x1)

  fit <- mcglm(y ~ mc(z, Pi) + x1, data = df, family = "binomial")
  expect_s3_class(fit, "mcglm")
  for (nm in names(fit$coefficients)) {
    expect_true(all(is.finite(fit$coefficients[[nm]])),
                info = paste("method:", nm))
  }
})

test_that("formula: binary Gaussian works", {
  set.seed(33)
  n <- 1000; Pi <- matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2)
  z <- rbinom(n, 1, 0.4)
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, 0.1)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.1)
  x1 <- rnorm(n)
  y <- 0.8 * z + -0.5 + 0.7 * x1 + rnorm(n, 0, 0.5)
  df <- data.frame(y = y, z = factor(z_hat), x1 = x1)

  fit <- mcglm(y ~ mc(z, Pi) + x1, data = df, family = "gaussian")
  expect_s3_class(fit, "mcglm")
  for (nm in names(fit$coefficients)) {
    expect_true(all(is.finite(fit$coefficients[[nm]])),
                info = paste("method:", nm))
  }
})


# ===== MULTICATEGORY - FORMULA =====

test_that("formula: multicategory K=3 all methods (auto pi_z)", {
  d <- gen_multi(); Pi <- d$Pi
  fit <- mcglm(y ~ mc(z, Pi) + x1, data = d$df, family = "poisson",
               method = c("naive", "bca", "bcm", "cs"))
  expect_equal(fit$K, 3)
  expect_named(fit$coefficients, c("naive", "bca", "bcm", "cs"))
  for (nm in names(fit$coefficients)) {
    expect_true(all(is.finite(fit$coefficients[[nm]])),
                info = paste("method:", nm))
  }
})

test_that("formula: multicategory K=4", {
  set.seed(77)
  n <- 2000; K <- 4
  pi_z <- rep(0.25, K)
  z <- sample(0:(K - 1), n, replace = TRUE, prob = pi_z)
  Pi <- matrix(0.05, K, K); diag(Pi) <- 0.85
  z_hat <- integer(n)
  for (i in 1:n) z_hat[i] <- sample(0:(K - 1), 1, prob = Pi[, z[i] + 1])
  x1 <- rnorm(n)
  y <- rpois(n, exp(c(0, 0.5, -0.3, 0.8)[z + 1] + 0.5 + 0.3 * x1))
  df <- data.frame(y = y, z = factor(z_hat), x1 = x1)

  fit <- mcglm(y ~ mc(z, Pi) + x1, data = df, family = "poisson",
               method = c("naive", "bca"))
  expect_equal(fit$K, K)
  expect_equal(fit$p, K - 1 + 2)
})

test_that("formula: multicategory with I(x^2)", {
  d <- gen_multi(); Pi <- d$Pi
  fit <- mcglm(y ~ mc(z, Pi) + x1 + I(x1^2), data = d$df,
               family = "poisson", method = "naive")
  expect_equal(fit$p, d$K - 1 + 3)  # 2 gamma + intercept + x1 + I(x1^2)
})


# ===== BIAS REDUCTION VIA FORMULA =====

test_that("formula: corrections reduce bias (binary Poisson)", {
  d <- gen_binary(n = 5000); Pi <- d$Pi
  fit <- mcglm(y ~ mc(z, Pi) + x1, data = d$df, family = "poisson")

  bias_naive <- abs(fit$coefficients$naive[1] - d$psi0[1])
  bias_bca   <- abs(fit$coefficients$bca[1] - d$psi0[1])
  bias_cs    <- abs(fit$coefficients$cs[1] - d$psi0[1])
  expect_lt(bias_bca, bias_naive)
  expect_lt(bias_cs, bias_naive)
})


# ===== ERROR HANDLING =====

test_that("formula: error when no mc() term", {
  d <- gen_binary(); Pi <- d$Pi
  expect_error(
    mcglm(y ~ x1 + x2, data = d$df, family = "poisson"),
    "mc\\(\\)"
  )
})

test_that("formula: error when variable not in data", {
  d <- gen_binary(); Pi <- d$Pi
  expect_error(
    mcglm(y ~ mc(missing_var, Pi) + x1, data = d$df, family = "poisson"),
    "not found"
  )
})

test_that("formula: error when Pi columns don't sum to 1", {
  d <- gen_binary(); Pi <- d$Pi
  bad_Pi <- matrix(c(0.9, 0.2, 0.1, 0.8), 2, 2)
  expect_error(
    mcglm(y ~ mc(z, bad_Pi) + x1, data = d$df, family = "poisson"),
    "sum to 1"
  )
})

test_that("formula: error when data is NULL", {
  expect_error(
    mcglm(y ~ mc(z, Pi) + x1, family = "poisson"),
    "data"
  )
})
