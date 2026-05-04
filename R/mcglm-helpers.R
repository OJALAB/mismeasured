# ---------------------------------------------------------------------------
# Internal helper functions for computing drift, Fisher info, Jacobian
# Used by mcglm() estimators
# ---------------------------------------------------------------------------

#' Estimate pi_z (true prevalence) from observed proxy and Pi
#'
#' Uses Bayesian inversion: if pi_obs = Pi %*% pi_z, then
#' pi_z = solve(Pi) %*% pi_obs. Clamps to [0.01, 0.99] for stability.
#' @param z_hat Integer vector of observed proxy values (0-based).
#' @param Pi K x K misclassification matrix.
#' @return Numeric vector of length K (estimated true prevalences).
#' @keywords internal
.mcglm_estimate_pi_z <- function(z_hat, Pi) {
  K <- nrow(Pi)
  # Observed proportions
  tab <- tabulate(z_hat + 1L, nbins = K)
  pi_obs <- tab / sum(tab)
  # Invert: pi_z = Pi^{-1} %*% pi_obs
  pi_z <- as.numeric(solve(Pi) %*% pi_obs)
  # Clamp to valid range

  pi_z <- pmax(pi_z, 0.01)
  pi_z <- pi_z / sum(pi_z)
  if (K == 2L) return(pi_z[2])  # scalar for binary case
  pi_z
}

#' Build xi_hat design matrix for mcglm
#'
#' For binary Z: xi_hat = cbind(z_hat, x).
#' For multicategory Z: xi_hat = cbind(d_hat, x) where d_hat is the
#' dummy-encoding of z_hat with baseline = 0.
#' @param z_hat Integer proxy covariate vector.
#' @param x Covariate matrix (n x r).
#' @param K Number of Z categories.
#' @return n x p matrix.
#' @keywords internal
.mcglm_build_xi_hat <- function(z_hat, x, K) {
  if (K == 2L) return(cbind(z_hat, x))
  n <- length(z_hat)
  s <- K - 1L
  d_hat <- matrix(0, n, s)
  for (k in seq_len(s)) d_hat[, k] <- as.numeric(z_hat == k)
  cbind(d_hat, x)
}

# ---- Binary misclassification helpers ----

#' Compute the toggle gap delta_i(psi) for binary misclassification
#' @keywords internal
.mcglm_compute_delta <- function(psi, x, mu_fun) {
  gamma <- psi[1]
  alpha <- psi[-1]
  eta0 <- as.numeric(x %*% alpha)
  mu_fun(gamma + eta0) - mu_fun(eta0)
}

#' Compute per-observation drift m_hat(psi) for binary misclassification
#' @keywords internal
.mcglm_compute_m_bin <- function(psi, x, mu_fun, c1, c2) {
  n <- nrow(x)
  p <- length(psi)
  delta <- .mcglm_compute_delta(psi, x, mu_fun)
  m <- matrix(0, n, p)
  m[, 1]  <- -c1 * delta
  m[, -1] <- -c2 * delta * x
  m
}

#' Compute m_hat(psi) = weighted mean of m_i
#' @keywords internal
.mcglm_compute_mhat_bin <- function(psi, x, mu_fun, c1, c2, wt = NULL) {
  m <- .mcglm_compute_m_bin(psi, x, mu_fun, c1, c2)
  if (is.null(wt)) return(colMeans(m))
  colSums(wt * m) / sum(wt)
}

#' Compute I_hat(psi) from pre-built xi_hat
#' @keywords internal
.mcglm_compute_Ihat <- function(psi, xi_hat, mu_dot_fun, wt = NULL) {
  n <- nrow(xi_hat)
  eta_tilde <- as.numeric(xi_hat %*% psi)
  w <- mu_dot_fun(eta_tilde)
  if (is.null(wt)) {
    crossprod(xi_hat * w, xi_hat) / n
  } else {
    crossprod(xi_hat * (wt * w), xi_hat) / sum(wt)
  }
}

#' Analytical Jacobian M_hat for binary misclassification
#' @keywords internal
.mcglm_compute_Mhat_bin <- function(psi, x, mu_fun, mu_dot_fun, c1, c2,
                                    wt = NULL) {
  gamma <- psi[1]
  alpha <- psi[-1]
  r <- length(alpha)
  p <- 1L + r
  n <- nrow(x)
  N <- if (is.null(wt)) n else sum(wt)

  eta0 <- as.numeric(x %*% alpha)
  mu_dot1 <- mu_dot_fun(gamma + eta0)
  mu_dot0 <- mu_dot_fun(eta0)
  d_mu_dot <- mu_dot1 - mu_dot0

  M <- matrix(0, p, p)

  if (is.null(wt)) {
    M[1, 1]    <- -c1 * mean(mu_dot1)
    M[1, 2:p]  <- -c1 * colMeans(d_mu_dot * x)
    M[2:p, 1]    <- -c2 * colMeans(mu_dot1 * x)
    M[2:p, 2:p]  <- -c2 * crossprod(d_mu_dot * x, x) / n
  } else {
    M[1, 1]    <- -c1 * sum(wt * mu_dot1) / N
    M[1, 2:p]  <- -c1 * colSums(wt * d_mu_dot * x) / N
    M[2:p, 1]    <- -c2 * colSums(wt * mu_dot1 * x) / N
    M[2:p, 2:p]  <- -c2 * crossprod(x * (wt * d_mu_dot), x) / N
  }

  M
}

# ---- Multicategory misclassification helpers ----

#' Compute per-observation drift m_i(psi) for multicategory misclassification
#' @keywords internal
.mcglm_compute_m_multi <- function(psi, x, K, mu_fun, Pi, pi_z) {
  n  <- nrow(x)
  s  <- K - 1
  r  <- ncol(x)
  p  <- s + r

  gamma <- c(0, psi[seq_len(s)])
  alpha <- psi[(s + 1):p]
  eta_base <- as.numeric(x %*% alpha)

  mu_mat <- sapply(seq_len(K), function(ell) {
    mu_fun(eta_base + gamma[ell])
  })

  m <- matrix(0, n, p)

  for (k in seq_len(s)) {
    for (ell in seq_len(K)) {
      coeff <- pi_z[ell] * Pi[k + 1, ell]
      m[, k] <- m[, k] + coeff * (mu_mat[, ell] - mu_mat[, k + 1])
    }
  }

  for (ell in seq_len(K)) {
    for (j in seq_len(K)) {
      coeff <- pi_z[ell] * Pi[j, ell]
      diff_mu <- mu_mat[, ell] - mu_mat[, j]
      m[, (s + 1):p] <- m[, (s + 1):p] + coeff * diff_mu * x
    }
  }

  m
}

#' @keywords internal
.mcglm_compute_mhat_multi <- function(psi, x, K, mu_fun, Pi, pi_z, wt = NULL) {
  m <- .mcglm_compute_m_multi(psi, x, K, mu_fun, Pi, pi_z)
  if (is.null(wt)) return(colMeans(m))
  colSums(wt * m) / sum(wt)
}

#' Compute I_hat(psi) for multicategory case
#' @keywords internal
.mcglm_compute_Ihat_multi <- function(psi, xi_hat, z_hat, K, mu_dot_fun,
                                      wt = NULL) {
  n <- nrow(xi_hat)
  s <- K - 1
  r <- ncol(xi_hat) - s
  gamma <- c(0, psi[seq_len(s)])
  alpha <- psi[(s + 1):(s + r)]
  eta_base <- as.numeric(xi_hat[, (s + 1):(s + r), drop = FALSE] %*% alpha)

  eta_tilde <- eta_base + gamma[z_hat + 1]
  w <- mu_dot_fun(eta_tilde)
  if (is.null(wt)) {
    crossprod(xi_hat * w, xi_hat) / n
  } else {
    crossprod(xi_hat * (wt * w), xi_hat) / sum(wt)
  }
}

#' Numerical Jacobian M_hat for multicategory drift
#' @keywords internal
.mcglm_compute_Mhat_multi <- function(psi, x, K, mu_fun, Pi, pi_z, wt = NULL) {
  p <- length(psi)
  M <- matrix(0, p, p)
  h <- 1e-7
  m0 <- .mcglm_compute_mhat_multi(psi, x, K, mu_fun, Pi, pi_z, wt = wt)
  for (j in seq_len(p)) {
    psi_h <- psi
    psi_h[j] <- psi_h[j] + h
    m1 <- .mcglm_compute_mhat_multi(psi_h, x, K, mu_fun, Pi, pi_z, wt = wt)
    M[, j] <- (m1 - m0) / h
  }
  M
}
