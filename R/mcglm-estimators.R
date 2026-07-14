# ---------------------------------------------------------------------------
# Estimators: naive, BCA, BCM, CS -- binary and multicategory
# Used by mcglm()
# ---------------------------------------------------------------------------

# ========================== BINARY CASE ==================================

#' Fit the naive GLM estimator (binary misclassification)
#' @keywords internal
.mcglm_fit_naive_bin <- function(y, xi_hat, family, wt = NULL) {
  fam <- .normalize_family(family)
  dat <- data.frame(y = y, xi_hat)
  if (is.null(wt)) {
    fit <- stats::glm(y ~ . - 1, data = dat, family = fam)
  } else {
    fit <- stats::glm(y ~ . - 1, data = dat, family = fam, weights = wt)
  }
  list(
    coefficients = unname(stats::coef(fit)),
    fitted       = stats::fitted(fit),
    glm_fit      = fit
  )
}

#' Additive bias correction (BCA) for binary misclassification
#' @keywords internal
.mcglm_fit_bca_bin <- function(psi_naive, y, xi_hat, x, family, c1, c2,
                               iterate = FALSE, max_iter = 50, tol = 1e-8,
                               wt = NULL) {
  fam <- .normalize_family(family)
  I_hat_inv <- solve(.mcglm_compute_Ihat(psi_naive, xi_hat, fam$mu.eta,
                                          wt = wt))

  m_hat <- .mcglm_compute_mhat_bin(psi_naive, x, fam$linkinv, c1, c2,
                                    wt = wt)
  psi   <- psi_naive - I_hat_inv %*% m_hat
  if (!iterate) return(as.numeric(psi))

  for (iter in seq_len(max_iter)) {
    m_hat   <- .mcglm_compute_mhat_bin(psi, x, fam$linkinv, c1, c2,
                                        wt = wt)
    psi_new <- psi_naive - I_hat_inv %*% m_hat
    if (max(abs(psi_new - psi)) < tol) break
    psi <- as.numeric(psi_new)
  }
  as.numeric(psi)
}

#' Multiplicative bias correction (BCM) for binary misclassification
#' @keywords internal
.mcglm_fit_bcm_bin <- function(psi_naive, y, xi_hat, x, family, c1, c2,
                               iterate = FALSE, max_iter = 50, tol = 1e-8,
                               wt = NULL) {
  fam <- .normalize_family(family)
  n <- length(y)
  N <- if (is.null(wt)) n else sum(wt)

  psi <- psi_naive
  for (iter in seq_len(if (iterate) max_iter else 1L)) {
    eta_tilde <- as.numeric(xi_hat %*% psi)
    resid <- y - fam$linkinv(eta_tilde)
    if (is.null(wt)) {
      U_hat <- colMeans(xi_hat * resid)
    } else {
      U_hat <- colSums(wt * xi_hat * resid) / N
    }
    m_hat <- .mcglm_compute_mhat_bin(psi, x, fam$linkinv, c1, c2, wt = wt)
    Phi   <- U_hat - m_hat

    I_hat <- .mcglm_compute_Ihat(psi, xi_hat, fam$mu.eta, wt = wt)
    M_hat <- .mcglm_compute_Mhat_bin(psi, x, fam$linkinv, fam$mu.eta, c1, c2,
                                      wt = wt)
    step  <- solve(I_hat + M_hat, Phi)

    psi_new <- psi + step
    if (iterate && max(abs(psi_new - psi)) < tol) { psi <- psi_new; break }
    psi <- psi_new
  }
  psi
}

#' Corrected-score estimator for binary misclassification
#' @keywords internal
.mcglm_fit_cs_bin <- function(psi_init, y, xi_hat, x, family, c1, c2,
                              wt = NULL) {
  fam <- .normalize_family(family)
  n   <- length(y)
  N   <- if (is.null(wt)) n else sum(wt)

  phi_mean <- function(psi) {
    eta_tilde <- as.numeric(xi_hat %*% psi)
    resid     <- y - fam$linkinv(eta_tilde)
    if (is.null(wt)) {
      score <- colMeans(xi_hat * resid)
    } else {
      score <- colSums(wt * xi_hat * resid) / N
    }
    m_hat <- .mcglm_compute_mhat_bin(psi, x, fam$linkinv, c1, c2, wt = wt)
    score - m_hat
  }

  phi_jac <- function(psi) {
    I_hat <- .mcglm_compute_Ihat(psi, xi_hat, fam$mu.eta, wt = wt)
    M_hat <- .mcglm_compute_Mhat_bin(psi, x, fam$linkinv, fam$mu.eta, c1, c2,
                                      wt = wt)
    -(I_hat + M_hat)
  }

  sol <- nleqslv::nleqslv(psi_init, phi_mean, jac = phi_jac,
                           control = list(maxit = 500, ftol = 1e-12))
  if (sol$termcd > 2) {
    warning("CS solver did not converge (termcd = ", sol$termcd, ")")
  }
  sol$x
}


# ========================= MULTICATEGORY CASE ============================

#' Fit the naive GLM estimator (multicategory misclassification)
#' @keywords internal
.mcglm_fit_naive_multi <- function(y, xi_hat, family, wt = NULL) {
  fam <- .normalize_family(family)
  dat <- data.frame(y = y, xi_hat)
  if (is.null(wt)) {
    fit <- stats::glm(y ~ . - 1, data = dat, family = fam)
  } else {
    fit <- stats::glm(y ~ . - 1, data = dat, family = fam, weights = wt)
  }
  list(
    coefficients = unname(stats::coef(fit)),
    fitted       = stats::fitted(fit),
    glm_fit      = fit
  )
}

#' BCA for multicategory misclassification
#' @keywords internal
.mcglm_fit_bca_multi <- function(psi_naive, y, xi_hat, z_hat, x, K, family,
                                 Pi, pi_z, iterate = FALSE, max_iter = 50,
                                 tol = 1e-8, wt = NULL) {
  fam <- .normalize_family(family)
  I_hat_inv <- solve(.mcglm_compute_Ihat_multi(psi_naive, xi_hat, z_hat, K,
                                                fam$mu.eta, wt = wt))

  m_hat <- .mcglm_compute_mhat_multi(psi_naive, x, K, fam$linkinv, Pi, pi_z,
                                      wt = wt)
  psi   <- psi_naive - I_hat_inv %*% m_hat
  if (!iterate) return(as.numeric(psi))

  for (iter in seq_len(max_iter)) {
    m_hat   <- .mcglm_compute_mhat_multi(psi, x, K, fam$linkinv, Pi, pi_z, wt = wt)
    psi_new <- psi_naive - I_hat_inv %*% m_hat
    if (max(abs(psi_new - psi)) < tol) break
    psi <- as.numeric(psi_new)
  }
  as.numeric(psi)
}

#' BCM for multicategory misclassification
#' @keywords internal
.mcglm_fit_bcm_multi <- function(psi_naive, y, xi_hat, z_hat, x, K, family,
                                 Pi, pi_z, iterate = FALSE, max_iter = 50,
                                 tol = 1e-8, wt = NULL,
                                 jacobian = c("analytical", "numerical")) {
  jacobian <- match.arg(jacobian)
  fam <- .normalize_family(family)
  n   <- length(y)
  N   <- if (is.null(wt)) n else sum(wt)
  s   <- K - 1
  r   <- ncol(x)

  psi <- psi_naive
  for (iter in seq_len(if (iterate) max_iter else 1L)) {
    gamma <- c(0, psi[seq_len(s)])
    alpha <- psi[(s + 1):(s + r)]
    eta_base  <- as.numeric(x %*% alpha)
    eta_tilde <- eta_base + gamma[z_hat + 1]
    resid <- y - fam$linkinv(eta_tilde)
    if (is.null(wt)) {
      U_hat <- colMeans(xi_hat * resid)
    } else {
      U_hat <- colSums(wt * xi_hat * resid) / N
    }
    m_hat <- .mcglm_compute_mhat_multi(psi, x, K, fam$linkinv, Pi, pi_z, wt = wt)
    Phi   <- U_hat - m_hat

    I_hat <- .mcglm_compute_Ihat_multi(psi, xi_hat, z_hat, K, fam$mu.eta,
                                        wt = wt)
    M_hat <- .mcglm_compute_Mhat_multi(psi, x, K, fam$linkinv, Pi, pi_z, wt = wt,
                                        jacobian = jacobian,
                                        mu_dot_fun = fam$mu.eta)
    step  <- solve(I_hat + M_hat, Phi)

    psi_new <- psi + step
    if (iterate && max(abs(psi_new - psi)) < tol) { psi <- psi_new; break }
    psi <- psi_new
  }
  psi
}

#' Corrected-score estimator for multicategory misclassification
#' @keywords internal
.mcglm_fit_cs_multi <- function(psi_init, y, xi_hat, z_hat, x, K, family,
                                Pi, pi_z, wt = NULL,
                                jacobian = c("analytical", "numerical")) {
  jacobian <- match.arg(jacobian)
  fam <- .normalize_family(family)
  n   <- length(y)
  N   <- if (is.null(wt)) n else sum(wt)
  s   <- K - 1
  r   <- ncol(x)

  phi_mean <- function(psi) {
    gamma <- c(0, psi[seq_len(s)])
    alpha <- psi[(s + 1):(s + r)]
    eta_base <- as.numeric(x %*% alpha)
    eta_tilde <- eta_base + gamma[z_hat + 1]
    resid <- y - fam$linkinv(eta_tilde)
    if (is.null(wt)) {
      score <- colMeans(xi_hat * resid)
    } else {
      score <- colSums(wt * xi_hat * resid) / N
    }
    m_hat <- .mcglm_compute_mhat_multi(psi, x, K, fam$linkinv, Pi, pi_z, wt = wt)
    score - m_hat
  }

  phi_jac <- function(psi) {
    I_hat <- .mcglm_compute_Ihat_multi(psi, xi_hat, z_hat, K, fam$mu.eta,
                                        wt = wt)
    M_hat <- .mcglm_compute_Mhat_multi(psi, x, K, fam$linkinv, Pi, pi_z, wt = wt,
                                        jacobian = jacobian,
                                        mu_dot_fun = fam$mu.eta)
    -(I_hat + M_hat)
  }

  sol <- nleqslv::nleqslv(psi_init, phi_mean, jac = phi_jac,
                           control = list(maxit = 500, ftol = 1e-12))
  if (sol$termcd > 2) {
    warning("CS solver did not converge (termcd = ", sol$termcd, ")")
  }
  sol$x
}
