# =========================================================================
# MC-SIMEX internal helpers
# =========================================================================
# These are used by .simex_discrete() in simex.R.
# No exported functions in this file.

#' Build design matrix [dummies(z, K-1), x] in R (mirrors C++ build_design)
#' @keywords internal
.build_xi_hat <- function(z_hat, x_mat, K) {
  n <- length(z_hat)
  s <- K - 1L
  dummies <- matrix(0, n, s)
  for (k in seq_len(s)) {
    dummies[, k] <- as.numeric(z_hat == k)
  }
  cbind(dummies, x_mat)
}

#' Generate parameter names for the SIMEX model
#' @keywords internal
.make_param_names <- function(z_levels, x_mat, K) {
  if (K == 2L) {
    dummy_names <- z_levels[2]
  } else {
    dummy_names <- z_levels[-1]
  }
  c(dummy_names, colnames(x_mat))
}


# =========================================================================
# Improved MC-SIMEX helpers (Sevilimedu & Yu, 2026)
# =========================================================================

#' Estimate marginal probability P(X=1) from observed data and Pi
#' @keywords internal
.estimate_pi_x <- function(z_hat, Pi) {
  p_w1 <- mean(z_hat == 1L)
  pi_00 <- Pi[1, 1]
  pi_11 <- Pi[2, 2]
  delta <- pi_00 + pi_11 - 1
  if (abs(delta) < 1e-10)
    stop("Misclassification matrix is degenerate (pi_00 + pi_11 = 1)")
  pi_x <- (p_w1 - 1 + pi_00) / delta
  max(1e-6, min(1 - 1e-6, pi_x))
}

#' Compute matrix power via eigendecomposition (R version)
#' @keywords internal
.mat_power_r <- function(Pi, power) {
  if (power == 0) return(diag(nrow(Pi)))
  if (power == 1) return(Pi)
  ev <- eigen(Pi)
  V <- ev$vectors
  d <- ev$values
  d_pow <- abs(d)^power * sign(d)
  zapsmall(V %*% diag(d_pow) %*% solve(V))
}

#' Compute correction factor c_lambda for the improved MC-SIMEX
#' @keywords internal
.compute_c_lambda <- function(Pi, pi_x, lambda) {
  sapply(lambda, function(lam) {
    Pi_lam <- .mat_power_r(Pi, 1 + lam)
    pi_00_lam <- Pi_lam[1, 1]
    pi_11_lam <- Pi_lam[2, 2]
    delta_lam <- pi_00_lam + pi_11_lam - 1
    if (abs(delta_lam) < 1e-10) return(NA_real_)
    denominator <- pi_x * (1 - pi_x) * delta_lam
    numerator <- (pi_00_lam - pi_x * delta_lam) *
                 (1 - pi_00_lam + delta_lam * pi_x)
    if (abs(denominator) < 1e-10) return(NA_real_)
    numerator / denominator
  })
}

#' Compute intercept correction factor for the improved MC-SIMEX
#' @keywords internal
.compute_intercept_correction <- function(Pi, pi_x, lambda) {
  sapply(lambda, function(lam) {
    Pi_lam <- .mat_power_r(Pi, 1 + lam)
    pi_00_lam <- Pi_lam[1, 1]
    pi_11_lam <- Pi_lam[2, 2]
    delta_lam <- pi_00_lam + pi_11_lam - 1
    denom <- pi_00_lam - delta_lam * pi_x
    if (abs(denom) < 1e-10) return(0)
    (1 - pi_11_lam) * pi_x / denom
  })
}

#' Find optimal lambda that minimizes |c_lambda|
#' @keywords internal
.find_optimal_lambda <- function(Pi, pi_x, grid = seq(0.1, 3, by = 0.01)) {
  c_values <- .compute_c_lambda(Pi, pi_x, grid)
  valid <- !is.na(c_values) & is.finite(c_values)
  if (!any(valid)) stop("Cannot find valid lambda for the improved method")
  grid[valid][which.min(abs(c_values[valid]))]
}

#' Variance estimation for the improved MC-SIMEX
#' @keywords internal
.variance_improved <- function(theta_list, corrected_coefs, c_lam_vec,
                               n_lambda, B, p, naive_fit, xi_hat, y, wt, fam) {
  N <- sum(wt)
  n <- length(y)
  if (n_lambda * B == 1) {
    psi_mean <- corrected_coefs
    eta_m <- as.numeric(xi_hat %*% psi_mean)
    w_m <- fam$mu_dot(eta_m)
    eps_m <- y - fam$mu(eta_m)
    A_m <- crossprod(xi_hat * (wt * w_m), xi_hat) / N
    C_m <- crossprod(xi_hat * (wt * eps_m), xi_hat * eps_m) / N
    A_m_inv <- tryCatch(solve(A_m), error = function(e) MASS::ginv(A_m))
    V_model <- A_m_inv %*% C_m %*% A_m_inv / N
    V_scaled <- V_model * c_lam_vec[1]^2
    return(V_scaled)
  }
  all_corrected <- do.call(rbind, theta_list)
  cov(all_corrected) / (n_lambda * B)
}
