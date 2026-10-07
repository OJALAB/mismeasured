# ---------------------------------------------------------------------------
# Latent prevalence from proxy frequencies
#
# Given the misclassification matrix Pi and the proxy distribution p, the
# latent prevalence solves p = Pi pi. The maximum-likelihood estimate over
# the simplex is found by the EM (maximum-likelihood label shift) iteration
#   pi_l <- sum_j p_j P(Z = l | Z_hat = j; Pi, pi)
# (Saerens, Latinne and Decaestecker, 2002). It equals Pi^{-1} p whenever
# that lies inside the simplex and otherwise returns the boundary MLE
# instead of clamping.
# ---------------------------------------------------------------------------

#' EM estimate of the latent prevalence
#'
#' @param p Proxy proportions (length K, summing to one).
#' @param Pi K x K column-stochastic misclassification matrix.
#' @param maxit,tol Iteration limit and convergence tolerance (maximum
#'   absolute change in the prevalence).
#' @param boundary_tol Prevalences below this are reported as boundary.
#' @return List with \code{pi}, \code{converged}, \code{iterations} and
#'   \code{boundary}.
#' @keywords internal
.mc_prevalence_em <- function(p, Pi, maxit = 1000L, tol = 1e-10,
                              boundary_tol = 1e-6) {
  K <- nrow(Pi)
  inv <- tryCatch(as.numeric(solve(Pi, p)), error = function(e) NULL)
  if (!is.null(inv) && all(is.finite(inv)) && all(inv > 0)) {
    pi <- inv / sum(inv)
  } else {
    start <- if (is.null(inv) || any(!is.finite(inv))) rep(1 / K, K) else
      pmax(inv, 0)
    if (sum(start) <= 0) start <- rep(1 / K, K)
    pi <- 0.5 * start / sum(start) + 0.5 / K
  }
  converged <- FALSE
  for (it in seq_len(maxit)) {
    joint <- sweep(Pi, 2, pi, "*")
    marg <- rowSums(joint)
    post <- joint / ifelse(marg > 0, marg, 1)
    pi_new <- colSums(p * post)
    pi_new <- pi_new / sum(pi_new)
    if (max(abs(pi_new - pi)) < tol) {
      pi <- pi_new
      converged <- TRUE
      break
    }
    pi <- pi_new
  }
  list(pi = pi, converged = converged, iterations = it,
       boundary = any(pi < boundary_tol))
}
