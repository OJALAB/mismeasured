# nonprob_mc(): point estimation, variance, nonprob-class compatibility

library(survey)

make_setup <- function(seed = 1, N = 50000, n_p = 1000, K = 2) {
  set.seed(seed)
  x1 <- rnorm(N)
  if (K == 2) {
    Pi <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
    z  <- rbinom(N, 1, 0.4)
    z_hat <- ifelse(z == 1, 1 - rbinom(N, 1, 0.15), rbinom(N, 1, 0.10))
    eta <- 0.8 * z - 0.5 + 0.7 * x1
  } else {
    Pi <- matrix(c(0.90, 0.06, 0.04,
                   0.08, 0.85, 0.07,
                   0.05, 0.05, 0.90), 3, 3)
    z  <- sample(0:2, N, TRUE, c(0.3, 0.4, 0.3))
    z_hat <- vapply(z, function(zz) sample(0:2, 1, prob = Pi[, zz + 1]), 0L)
    eta <- c(0, 0.8, -0.5)[z + 1] - 0.2 + 0.4 * x1
  }
  y <- rpois(N, exp(eta))
  s_np  <- rbinom(N, 1, plogis(-2.5 + 0.5 * x1)) == 1
  df_np <- data.frame(y = y, z = z_hat, x1 = x1)[s_np, ]
  s_p   <- sample.int(N, n_p)
  des   <- svydesign(ids = ~1, weights = rep(N / n_p, n_p),
                     data = data.frame(z = z[s_p], x1 = x1[s_p]))
  list(df_np = df_np, des = des, Pi = Pi, ybar = mean(y),
       pi_z = if (K == 2) 0.4 else c(0.3, 0.4, 0.3))
}

test_that("nonprob_mc reproduces the manual public-API pipeline (binary)", {
  s <- make_setup()
  Pi <- s$Pi
  fit <- nonprob_mc(data = s$df_np, outcome = y ~ mc(z, Pi) + x1,
                    svydesign = s$des, family_outcome = "poisson",
                    mc_method = "cs", Pi = Pi)

  # manual pipeline (as in the mismeasured POC)
  m  <- mismeasured::mcglm(y ~ mc(z, Pi) + x1, data = s$df_np,
                           family = "poisson", method = c("naive", "cs"),
                           Pi = Pi)
  pr <- predict(m, newdata = s$des$variables, type = "response",
                method = "cs")
  des_upd <- update(s$des, .p = pr)
  mn <- svymean(~.p, des_upd)
  expect_equal(fit$output$mean, as.numeric(mn), tolerance = 1e-10)
  expect_equal(fit$SE$prob^2, as.numeric(attr(mn, "var")), tolerance = 1e-10)

  B <- sandwich::bread(m, method = "cs")
  S <- sandwich::meat(m, method = "cs")
  w <- weights(s$des)
  xi <- cbind(s$des$variables$z, 1, s$des$variables$x1)
  G <- colSums(w * poisson()$mu.eta(as.numeric(xi %*% coef(m))) * xi) / sum(w)
  v_np <- as.numeric(t(G) %*% B %*% S %*% t(B) %*% G) / nobs(m)
  expect_equal(fit$SE$nonprob^2, v_np, tolerance = 1e-10)

  # sanity: corrected estimate is close to the truth; z-coefficient is
  # de-attenuated relative to the naive fit
  expect_lt(abs(fit$output$mean - s$ybar), 4 * fit$output$SE)
  co <- fit$outcome[[1]]$model_fitted$coefficients
  expect_gt(co$cs[1], co$naive[1])
})

test_that("nonprobsvy S3 methods work on the returned object", {
  s <- make_setup()
  fit <- nonprob_mc(data = s$df_np, outcome = y ~ mc(z, Pi) + x1,
                    svydesign = s$des, family_outcome = "poisson",
                    Pi = s$Pi)
  expect_s3_class(fit, "nonprob")
  expect_output(print(fit), "mass imputation")
  expect_output(print(fit), "mcglm \\(poisson, cs\\)")
  expect_output(print(summary(fit)), "outcome residuals")
  ci <- confint(fit)
  expect_true(ci$lower_bound < fit$output$mean &&
              fit$output$mean < ci$upper_bound)
  expect_equal(unname(nobs(fit)), c(1000, nrow(s$df_np)))
  expect_equal(nonprobsvy::pop_size(fit), sum(weights(s$des)))
})

test_that("multicategory (K = 3) works end to end", {
  s <- make_setup(seed = 4, K = 3)
  fit <- nonprob_mc(data = s$df_np, outcome = y ~ mc(z, Pi) + x1,
                    svydesign = s$des, family_outcome = "poisson",
                    mc_method = "cs", Pi = s$Pi, pi_z = s$pi_z)
  expect_true(is.finite(fit$output$mean) && is.finite(fit$output$SE))
  expect_lt(abs(fit$output$mean - s$ybar), 4 * fit$output$SE)
  expect_length(fit$outcome[[1]]$coefficients, 4L)  # gamma1, gamma2, alpha0, alpha1
})

test_that("mc_method selects the estimator; se = FALSE skips variance", {
  s <- make_setup()
  fit_n <- nonprob_mc(data = s$df_np, outcome = y ~ mc(z, Pi) + x1,
                      svydesign = s$des, family_outcome = "poisson",
                      mc_method = "naive", Pi = s$Pi, se = FALSE)
  fit_c <- nonprob_mc(data = s$df_np, outcome = y ~ mc(z, Pi) + x1,
                      svydesign = s$des, family_outcome = "poisson",
                      mc_method = "cs", Pi = s$Pi, se = FALSE)
  expect_false(fit_n$output$mean == fit_c$output$mean)
  expect_true(is.na(fit_n$output$SE))
})

test_that("informative errors on bad input", {
  s <- make_setup()
  expect_error(nonprob_mc(data = s$df_np, outcome = y ~ z + x1,
                          svydesign = s$des, family_outcome = "poisson"),
               "mc\\(\\) term")
  des_noz <- svydesign(ids = ~1, weights = rep(1, 50),
                       data = data.frame(x1 = rnorm(50)))
  expect_error(nonprob_mc(data = s$df_np, outcome = y ~ mc(z, Pi) + x1,
                          svydesign = des_noz, family_outcome = "poisson",
                          Pi = s$Pi),
               "TRUE category")
  expect_error(nonprob_mc(data = s$df_np, outcome = y ~ mc(z, Pi) + x1,
                          svydesign = "not a design",
                          family_outcome = "poisson"),
               "svydesign")
})
