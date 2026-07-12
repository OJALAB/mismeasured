# estfun()/bread() sandwich integration, coeftest, and mc_parse_formula()

skip_if_not_installed("sandwich")

make_fit_bin <- function(n = 800, seed = 7,
                         method = c("naive", "bca", "bcm", "cs")) {
  set.seed(seed)
  x_mat <- cbind(1, rnorm(n))
  z     <- rbinom(n, 1, 0.4)
  y     <- rpois(n, exp(0.8 * z - 0.5 + 0.7 * x_mat[, 2]))
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, 0.10)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)
  mcglm(y, z_hat = z_hat, x = x_mat, family = "poisson", method = method,
        p01 = 0.10, p10 = 0.15, pi_z = 0.4)
}

# manual sandwich assembly: B S B' / n (sandwich::sandwich() itself omits
# the transpose and does not forward `method` to bread -- see ?estfun.mcglm)
manual_sandwich <- function(fit, m) {
  B <- sandwich::bread(fit, method = m)
  S <- sandwich::meat(fit, method = m)
  B %*% S %*% t(B) / nobs(fit)
}

test_that("estfun: cs estimating equations are solved, naive score is zero", {
  fit <- make_fit_bin()
  expect_equal(colMeans(sandwich::estfun(fit, method = "cs")),
               rep(0, 3), ignore_attr = TRUE, tolerance = 1e-6)
  expect_equal(colMeans(sandwich::estfun(fit, method = "naive")),
               rep(0, 3), ignore_attr = TRUE, tolerance = 1e-6)
})

test_that("bread/meat assembly reproduces vcov() for every method", {
  fit <- make_fit_bin()
  for (m in c("naive", "bca", "bcm", "cs")) {
    expect_equal(unname(manual_sandwich(fit, m)),
                 unname(vcov(fit, method = m)),
                 tolerance = 1e-10, label = m)
  }
})

test_that("sandwich::sandwich() matches vcov() for symmetric-bread fits", {
  # single-method fit so that sandwich()'s default bread/meat agree
  fit <- make_fit_bin(method = "naive")
  expect_equal(unname(sandwich::sandwich(fit)),
               unname(vcov(fit, method = "naive")), tolerance = 1e-10)
})

test_that("estfun/bread work for multicategory fits", {
  set.seed(11)
  n  <- 1200
  Pi <- matrix(c(0.90, 0.06, 0.04,
                 0.08, 0.85, 0.07,
                 0.05, 0.05, 0.90), 3, 3)
  pi_z <- c(0.3, 0.4, 0.3)
  z  <- sample(0:2, n, TRUE, pi_z)
  z_hat <- vapply(z, function(zz) sample(0:2, 1, prob = Pi[, zz + 1]), 0L)
  x_mat <- cbind(1, rnorm(n))
  y  <- rpois(n, exp(0.2 + 0.4 * x_mat[, 2] + c(0, 0.8, -0.5)[z + 1]))
  fit <- mcglm(y, z_hat = z_hat, x = x_mat, family = "poisson",
               method = c("naive", "cs"), Pi = Pi, pi_z = pi_z)

  expect_equal(colMeans(sandwich::estfun(fit, method = "cs")),
               rep(0, 4), ignore_attr = TRUE, tolerance = 1e-6)
  for (m in c("naive", "cs")) {
    expect_equal(unname(manual_sandwich(fit, m)),
                 unname(vcov(fit, method = m)),
                 tolerance = 1e-10, label = m)
  }
})

test_that("lmtest::coeftest works on mcglm fits", {
  skip_if_not_installed("lmtest")
  fit <- make_fit_bin()
  ct <- lmtest::coeftest(fit)   # coef() + vcov() of the default method
  expect_equal(unname(ct[, "Estimate"]),
               unname(coef(fit, method = "cs")))
  expect_equal(unname(ct[, "Std. Error"]),
               unname(sqrt(diag(vcov(fit, method = "cs")))))
  ct_naive <- lmtest::coeftest(
    fit, vcov. = function(x, ...) vcov(x, method = "naive"))
  expect_equal(unname(ct_naive[, "Std. Error"]),
               unname(sqrt(diag(vcov(fit, method = "naive")))))
})

test_that("coeftest, mc_parse_formula, predict(newdata), marginaleffects work for K = 3", {
  set.seed(21)
  n  <- 900
  Pi <- matrix(c(0.90, 0.06, 0.04,
                 0.08, 0.85, 0.07,
                 0.05, 0.05, 0.90), 3, 3)
  pi_z <- c(0.3, 0.4, 0.3)
  z  <- sample(0:2, n, TRUE, pi_z)
  z_hat <- vapply(z, function(zz) sample(0:2, 1, prob = Pi[, zz + 1]), 0L)
  x1 <- rnorm(n)
  y  <- rpois(n, exp(0.2 + 0.4 * x1 + c(0, 0.8, -0.5)[z + 1]))
  df <- data.frame(y = y, z = z_hat, x1 = x1)

  parts <- mc_parse_formula(y ~ mc(z, Pi) + x1, df)
  expect_equal(parts$K, 3L)
  expect_equal(parts$Pi, Pi)

  fit <- mcglm(y ~ mc(z, Pi) + x1, data = df, family = "poisson",
               method = c("naive", "cs"), pi_z = pi_z)

  # coeftest: z-tests from coef()/vcov() of the default (cs) method
  skip_if_not_installed("lmtest")
  ct <- lmtest::coeftest(fit)
  expect_equal(unname(ct[, "Estimate"]), unname(coef(fit, method = "cs")))
  expect_equal(nrow(ct), 4L)   # gamma1, gamma2, alpha0, alpha1

  # predict(newdata): category contrasts at x1 = 0 recover the gammas
  nd <- data.frame(z = 0:2, x1 = 0)
  pr <- predict(fit, newdata = nd, type = "link")
  expect_equal(pr[2] - pr[1], unname(coef(fit)["gamma1"]))
  expect_equal(pr[3] - pr[1], unname(coef(fit)["gamma2"]))
  expect_equal(predict(fit, newdata = df, type = "response"),
               predict(fit, type = "response"))

  # marginaleffects: average category comparisons vs manual predictions
  # (numeric custom contrasts take exactly two values -> one pair per call)
  skip_if_not_installed("marginaleffects")
  mu <- function(zv) mean(predict(fit, type = "response",
                                  newdata = transform(df, z = zv)))
  for (k in 1:2) {
    cmp <- marginaleffects::avg_comparisons(fit,
                                            variables = list(z = c(0, k)),
                                            newdata = df)
    expect_equal(cmp$estimate, mu(k) - mu(0), tolerance = 1e-8)
    expect_true(is.finite(cmp$std.error))
  }
})

test_that("mc_parse_formula() returns the mcglm() inputs", {
  set.seed(3)
  Pi <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
  df <- data.frame(y = rpois(100, 1), z = rbinom(100, 1, 0.4),
                   x1 = rnorm(100))
  parts <- mc_parse_formula(y ~ mc(z, Pi) + x1, df)
  expect_named(parts, c("y", "z_hat", "x", "Pi", "K"))
  expect_equal(parts$K, 2L)
  expect_equal(parts$z_hat, df$z)
  expect_equal(ncol(parts$x), 2L)
  expect_equal(parts$Pi, Pi)

  fit_f <- mcglm(y ~ mc(z, Pi) + x1, data = df, family = "poisson",
                 method = "cs")
  fit_m <- mcglm(parts$y, z_hat = parts$z_hat, x = parts$x,
                 family = "poisson", method = "cs", Pi = parts$Pi)
  expect_equal(fit_f$coefficients$cs, fit_m$coefficients$cs)
})

test_that("marginaleffects works on formula-interface mcglm fits", {
  skip_if_not_installed("marginaleffects")
  set.seed(9)
  n  <- 400
  Pi <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
  z  <- rbinom(n, 1, 0.4)
  x1 <- rnorm(n)
  y  <- rpois(n, exp(0.8 * z - 0.5 + 0.7 * x1))
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, 0.10)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)
  df  <- data.frame(y = y, z = z_hat, x1 = x1)
  fit <- mcglm(y ~ mc(z, Pi) + x1, data = df, family = "poisson",
               method = c("naive", "cs"))

  cmp <- marginaleffects::avg_comparisons(fit, variables = list(z = 0:1),
                                          newdata = df)
  expect_s3_class(cmp, "comparisons")
  expect_true(is.finite(cmp$estimate) && is.finite(cmp$std.error))
  # the z-contrast must reflect the corrected (default = cs) coefficients
  mu  <- function(zv) mean(predict(fit, type = "response",
                                   newdata = transform(df, z = zv)))
  expect_equal(cmp$estimate, mu(1) - mu(0), tolerance = 1e-8)

  sl <- marginaleffects::avg_slopes(fit, variables = "x1", newdata = df)
  expect_true(is.finite(sl$estimate) && is.finite(sl$std.error))

  # the plain user-facing call, no newdata: the original data must be
  # recoverable from the stored call
  sl0 <- marginaleffects::avg_slopes(fit)
  expect_true(all(is.finite(sl0$estimate)) && all(is.finite(sl0$std.error)))
  expect_setequal(sl0$term, c("z", "x1"))
})

test_that("predict(newdata) works for formula and matrix fits", {
  set.seed(5)
  Pi <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
  df <- data.frame(y = rpois(300, 1), z = rbinom(300, 1, 0.4),
                   x1 = rnorm(300))
  fit <- mcglm(y ~ mc(z, Pi) + x1, data = df, family = "poisson",
               method = "cs")

  # in-sample newdata reproduces in-sample predictions
  expect_equal(predict(fit, newdata = df, type = "response"),
               predict(fit, type = "response"))
  # newdata without the response column is accepted
  nd <- data.frame(z = c(0L, 1L), x1 = c(0, 0))
  pr <- predict(fit, newdata = nd, type = "link")
  expect_equal(pr[2] - pr[1], unname(coef(fit, method = "cs")[1]))

  # matrix-interface fits take list(z_hat, x)
  fit_m <- mcglm(df$y, z_hat = df$z, x = cbind(1, df$x1),
                 family = "poisson", method = "cs", Pi = Pi)
  pr_m <- predict(fit_m, newdata = list(z_hat = nd$z, x = cbind(1, nd$x1)),
                  type = "link")
  expect_equal(pr_m, pr, tolerance = 1e-8)
})
