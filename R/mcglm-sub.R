# ---------------------------------------------------------------------------
# Subtraction-correction estimator ("sub")
#
# Reference: Yi, G. Y., Yan, Y., Liao, X. and Spiegelman, D. (2019).
# Parametric regression analysis with covariate misclassification in main
# study/validation study designs. Int. J. Biostat., 15(1), 20170002,
# Section 3.2, eq. (9).
#
# The naive score S~(psi) = xi_hat {Y - mu(psi' xi_hat)} is corrected by
# subtracting its conditional expectation given the observed (z_hat, x):
#   U_sub = S~ - E[S~ | z_hat, x] = xi_hat {Y - nu(z_hat, x; psi)},
#   nu(j, x; psi) = sum_l P(Z = l | z_hat = j) mu(gamma_l + alpha' x).
# For canonical-link GLMs this needs only the mean model and the predictive
# matrix W[j + 1, l + 1] = P(Z = l | z_hat = j), so it covers any K, the
# three supported families and no dispersion parameter.
#
# The package's "cs" estimator subtracts E[S~ | x] instead; conditioning on
# the finer sigma-field (z_hat, x) gives the same Jacobian and a smaller
# middle matrix, so "sub" is never less efficient than "cs" when the
# probabilities are known.
# ---------------------------------------------------------------------------

#' Expand a prevalence argument to the full K-vector
#' @keywords internal
.mcglm_pi_vector <- function(pi_z, K) {
  pi_z <- as.numeric(pi_z)
  if (K == 2L && length(pi_z) == 1L) c(1 - pi_z, pi_z) else pi_z
}

#' Predictive matrix W[j + 1, l + 1] = P(Z = l | Z_hat = j)
#'
#' @param Pi K x K column-stochastic misclassification matrix.
#' @param pi_z Prevalence (scalar for K = 2 or K-vector).
#' @return K x K matrix whose rows sum to one.
#' @keywords internal
.mcglm_predictive_matrix <- function(Pi, pi_z) {
  Pi <- as.matrix(Pi)
  K  <- nrow(Pi)
  joint <- sweep(Pi, 2, .mcglm_pi_vector(pi_z, K), "*")
  marg  <- rowSums(joint)
  if (any(!is.finite(marg)) || any(marg <= 0))
    stop("P(Z_hat = j) implied by Pi and pi_z is zero for proxy category ",
         paste(which(!(marg > 0)) - 1L, collapse = ", "),
         "; the predictive probabilities P(Z | Z_hat) are undefined.",
         call. = FALSE)
  joint / marg
}

#' Per-class mean functions evaluated at every row
#'
#' @return n x K matrix with column l + 1 equal to fun(gamma_l + alpha'x),
#'   gamma_0 = 0.
#' @keywords internal
.mcglm_class_eval <- function(psi, x, K, fun) {
  s <- K - 1L
  n <- nrow(x)
  gamma <- c(0, psi[seq_len(s)])
  eta <- as.numeric(x %*% psi[-seq_len(s)])
  out <- matrix(0, n, K)
  for (l in seq_len(K)) out[, l] <- fun(eta + gamma[l])
  out
}

#' Per-observation subtraction-corrected score and its derivative pieces
#'
#' @param Wi n x K matrix of \eqn{\Pr(Z = \ell \mid \hat Z_i, x_i)} (rows of
#'   the predictive matrix, or row-specific under a prevalence model).
#' @return List with \code{U} (n x p estimating-function rows),
#'   \code{D} (n x p matrix of d nu_i / d psi), \code{mu} (n x K class
#'   means) and \code{nu} (length-n corrected means).
#' @keywords internal
.mcglm_sub_parts <- function(psi, y, xi_hat, x, K, fam, Wi) {
  mu  <- .mcglm_class_eval(psi, x, K, fam$linkinv)
  mud <- .mcglm_class_eval(psi, x, K, fam$mu.eta)
  nu  <- rowSums(Wi * mu)
  D   <- cbind(Wi[, -1L, drop = FALSE] * mud[, -1L, drop = FALSE],
               rowSums(Wi * mud) * x)
  list(U = xi_hat * (y - nu), D = D, mu = mu, nu = nu)
}

#' Solve the subtraction-corrected estimating equation
#'
#' @param W Predictive matrix from \code{.mcglm_predictive_matrix}.
#' @return List with \code{coefficients}, \code{converged}, \code{termcd}
#'   and \code{iterations}.
#' @keywords internal
.mcglm_fit_sub <- function(psi_init, y, xi_hat, z_hat, x, K, family, W,
                           wt = NULL) {
  fam <- .normalize_family(family)
  w   <- if (is.null(wt)) rep(1, length(y)) else wt
  N   <- sum(w)
  Wi  <- W[z_hat + 1L, , drop = FALSE]

  score_mean <- function(psi)
    colSums(w * .mcglm_sub_parts(psi, y, xi_hat, x, K, fam, Wi)$U) / N
  score_jac <- function(psi) {
    D <- .mcglm_sub_parts(psi, y, xi_hat, x, K, fam, Wi)$D
    -crossprod(xi_hat * w, D) / N
  }

  sol <- nleqslv::nleqslv(psi_init, score_mean, jac = score_jac,
                          control = list(maxit = 500, ftol = 1e-12))
  if (sol$termcd > 2)
    warning("SUB solver did not converge (termcd = ", sol$termcd, ")")
  list(coefficients = sol$x, converged = sol$termcd <= 2,
       termcd = sol$termcd, iterations = sol$iter)
}

#' Subtraction correction with estimated misclassification probabilities
#'
#' Solves the SUB equation with the predictive matrix from an
#' \code{\link{estimate_mc}} object (true-category score on internal
#' validation rows, per \code{control$beta_equation}) and returns the
#' stacked \eqn{(\psi, \eta)} sandwich via \code{.mcglm_fit_validated}.
#' @keywords internal
.mcglm_fit_sub_validated <- function(psi_init, y, xi_hat, z_hat, x, K,
                                     family, est, wt = NULL,
                                     control = control_mc()) {
  fam <- .normalize_family(family)
  U_fun <- function(psi, par)
    .mcglm_sub_parts(psi, y, xi_hat, x, K, fam, par$Wi)$U
  J_fun <- function(psi, par, w)
    -crossprod(xi_hat * w,
               .mcglm_sub_parts(psi, y, xi_hat, x, K, fam, par$Wi)$D)
  .mcglm_fit_validated(psi_init, est, y, x, z_hat, K, fam, U_fun, J_fun,
                       wt = wt, control = control, label = "SUB")
}

#' Sandwich variance of the subtraction-corrected estimator
#'
#' Treats the predictive matrix \code{W} as known:
#' \eqn{J^{-1} S J^{-\top} / N} with \eqn{J = -N^{-1}\sum_i \hat\xi_i
#' \partial\nu_i/\partial\psi^\top} and \eqn{S = N^{-1}\sum_i U_i U_i^\top}.
#' @keywords internal
.mcglm_vcov_sub <- function(psi, y, xi_hat, z_hat, x, K, family, W,
                            wt = NULL) {
  fam <- .normalize_family(family)
  w   <- if (is.null(wt)) rep(1, length(y)) else wt
  N   <- sum(w)
  parts <- .mcglm_sub_parts(psi, y, xi_hat, x, K, fam,
                            W[z_hat + 1L, , drop = FALSE])
  S <- crossprod(parts$U * w, parts$U) / N
  J <- -crossprod(xi_hat * w, parts$D) / N
  J_inv <- solve(J)
  J_inv %*% S %*% t(J_inv) / N
}
