# mcglm() with validation_sample() / estimate_mc() objects: entry points,
# design weights, Hajek vs HT, beta equations and variance options.

.sd_data <- function(n = 1500L, seed = 301L) {
  set.seed(seed)
  z <- rbinom(n, 1L, 0.4)
  z_hat <- ifelse(z == 1L, rbinom(n, 1L, 0.85), rbinom(n, 1L, 0.12))
  x1 <- rnorm(n)
  y <- rpois(n, exp(-0.3 + 0.8 * z + 0.5 * x1))
  list(y = y, z = z, z_hat = z_hat, x1 = x1, x = cbind(1, x1),
       df = data.frame(y = y, z = z_hat, x1 = x1))
}

test_that("list, validation_sample, estimate_mc and mc(z, est) agree", {
  d <- .sd_data()
  idx <- sample.int(length(d$y), 300L)
  val <- validation_sample(d$z[idx], index = idx)
  est <- estimate_mc(val, z_hat = d$z_hat)
  fits <- list(
    mcglm(y ~ mc(z) + x1, data = d$df, method = "sub",
          validation = list(z = d$z[idx], index = idx)),
    mcglm(y ~ mc(z) + x1, data = d$df, method = "sub", validation = val),
    mcglm(y ~ mc(z) + x1, data = d$df, method = "sub", validation = est),
    mcglm(y ~ mc(z, est) + x1, data = d$df, method = "sub"))
  for (f in fits[-1]) {
    expect_equal(coef(f, method = "sub"), coef(fits[[1]], method = "sub"))
    expect_equal(vcov(f, method = "sub"), vcov(fits[[1]], method = "sub"))
  }
  expect_s3_class(fits[[4]]$mc_estimate, "mc_estimate")
  expect_identical(fits[[4]]$nuisance$sub$beta_equation, "yi")
})

test_that("estimate objects are checked against the fitted data", {
  d <- .sd_data()
  idx <- 1:200
  est <- estimate_mc(validation_sample(d$z[idx], index = idx),
                     z_hat = d$z_hat)
  other <- d$df
  other$z[1:5] <- 1L - other$z[1:5]
  expect_error(mcglm(y ~ mc(z, est) + x1, data = other, method = "sub"),
               "different main-study proxies")
  expect_error(mcglm(y ~ mc(z, est) + x1, data = d$df, method = "sub",
                     validation = list(z = d$z[idx], index = idx)),
               "not both")
  expect_error(mcglm(y ~ mc(z, est) + x1, data = d$df, method = "sub",
                     pi_z = 0.4), "remove pi_z")
  expect_error(mcglm(y ~ mc(z, est) + x1, data = d$df, method = "sub", K = 3),
               "differs")
  # an external estimate with validation prevalence is portable
  ext <- estimate_mc(validation_sample(d$z[1:300], d$z_hat[1:300]))
  expect_silent(mcglm(y ~ mc(z, ext) + x1, data = other, method = "sub"))
  expect_error(simex(y ~ mc(z, ext) + x1, data = d$df, family = poisson()),
               "does not accept estimate_mc")
})

test_that("Hajek and HT differ only through the prevalence", {
  d <- .sd_data(seed = 302L)
  n <- length(d$y)
  p_incl <- ifelse(d$z_hat == 1L, 0.25, 0.06)
  idx <- which(runif(n) < p_incl)
  mk <- function(estimator)
    mcglm(d$y, z_hat = d$z_hat, x = d$x, method = "sub",
          validation = validation_sample(d$z[idx], index = idx,
                                         weights = 1 / p_incl[idx],
                                         estimator = estimator))
  fh <- mk("hajek")
  ft <- mk("ht")
  expect_equal(fh$nuisance$sub$Pi, ft$nuisance$sub$Pi)
  expect_false(isTRUE(all.equal(fh$nuisance$sub$pi_z, ft$nuisance$sub$pi_z)))
  expect_identical(fh$nuisance$sub$beta_equation, "weighted")
  expect_true(all(is.finite(c(fh$se$sub, ft$se$sub))))
})

test_that("beta_equation = 'weighted' solves the design-weighted equation", {
  d <- .sd_data(n = 800L, seed = 303L)
  p_incl <- ifelse(d$y == 0, 0.08, 0.3)
  idx <- which(runif(800L) < p_incl)
  dw <- 1 / p_incl[idx]
  fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, method = "sub",
               validation = validation_sample(d$z[idx], index = idx,
                                              weights = dw))
  psi <- unname(coef(fit, method = "sub"))
  W <- fit$nuisance$sub$W
  mu <- exp(outer(drop(d$x %*% psi[2:3]), c(0, psi[1]), "+"))
  xi_hat <- cbind(d$z_hat, d$x)
  U <- xi_hat * (d$y - rowSums(W[d$z_hat + 1L, ] * mu))
  xi <- cbind(d$z[idx], d$x[idx, ])
  S <- xi * (d$y[idx] - exp(drop(xi %*% psi)))
  expect_lt(max(abs(colSums(U) + colSums(dw * (S - U[idx, ])))), 1e-6)
})

test_that("conditional variance ignores the nuisance uncertainty", {
  d <- .sd_data(seed = 304L)
  v <- validation_sample(d$z[1:150], d$z_hat[1:150])
  fd <- mcglm(d$y, z_hat = d$z_hat, x = d$x, method = "sub", validation = v)
  fc <- mcglm(d$y, z_hat = d$z_hat, x = d$x, method = "sub", validation = v,
              mc_control = control_mc(variance = "conditional"))
  expect_equal(coef(fd, method = "sub"), coef(fc, method = "sub"))
  expect_lt(fc$se$sub[1], fd$se$sub[1])
  # known probabilities give the same conditional sandwich
  fk <- mcglm(d$y, z_hat = d$z_hat, x = d$x, method = "sub",
              Pi = fd$nuisance$sub$Pi, pi_z = fd$nuisance$sub$pi_z[2])
  expect_equal(unname(vcov(fc, method = "sub")),
               unname(vcov(fk, method = "sub")), tolerance = 1e-8)
})

test_that("stratified external audit and cs restrictions", {
  d <- .sd_data(seed = 305L)
  ext <- .sd_data(n = 3000L, seed = 306L)
  idx <- unlist(lapply(0:1, function(j) sample(which(ext$z_hat == j), 80L)))
  v <- validation_sample(ext$z[idx], ext$z_hat[idx], strata = "z_hat",
                         N_strata = c(`0` = sum(ext$z_hat == 0),
                                      `1` = sum(ext$z_hat == 1)))
  fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, method = "sub", validation = v)
  expect_true(all(is.finite(fit$se$sub)))
  expect_error(mcglm(d$y, z_hat = d$z_hat, x = d$x, method = c("cs", "sub"),
                     validation = v), "unweighted, unstratified")
})

test_that("design-weighted SUB is unbiased and calibrated under Y-dependent audits", {
  skip_on_cran()
  reps <- 150L
  out <- matrix(NA_real_, reps, 6L)
  for (r in seq_len(reps)) {
    d <- .sd_data(n = 3000L, seed = 2000L + r)
    p_incl <- ifelse(d$y == 0, 0.03, 0.25)
    idx <- which(runif(3000L) < p_incl)
    fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, method = "sub",
                 validation = validation_sample(d$z[idx], index = idx,
                                                weights = 1 / p_incl[idx]),
                 mc_control = control_mc(on_ill = "none"))
    out[r, ] <- c(coef(fit, method = "sub"), fit$se$sub)
  }
  expect_lt(max(abs(colMeans(out[, 1:3]) - c(0.8, -0.3, 0.5))), 0.03)
  ratio <- colMeans(out[, 4:6]) / apply(out[, 1:3], 2, sd)
  expect_true(all(ratio > 0.8 & ratio < 1.25))
})

test_that("external z_hat-stratified audits give calibrated SEs", {
  skip_on_cran()
  reps <- 150L
  out <- matrix(NA_real_, reps, 2L)
  for (r in seq_len(reps)) {
    d <- .sd_data(n = 3000L, seed = 3000L + r)
    frame <- .sd_data(n = 10000L, seed = 5000L + r)
    idx <- unlist(lapply(0:1, function(j)
      sample(which(frame$z_hat == j), 150L)))
    v <- validation_sample(frame$z[idx], frame$z_hat[idx], strata = "z_hat",
                           N_strata = c(`0` = sum(frame$z_hat == 0),
                                        `1` = sum(frame$z_hat == 1)))
    fit <- mcglm(d$y, z_hat = d$z_hat, x = d$x, method = "sub",
                 validation = v, mc_control = control_mc(on_ill = "none"))
    out[r, ] <- c(coef(fit, method = "sub")[1], fit$se$sub[1])
  }
  expect_lt(abs(mean(out[, 1]) - 0.8), 0.03)
  ratio <- mean(out[, 2]) / sd(out[, 1])
  expect_true(ratio > 0.8 && ratio < 1.25)
})
