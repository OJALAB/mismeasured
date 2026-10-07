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

# ---------------------------------------------------------------------------
# Covariate-dependent prevalence P(Z = l | x) (Yi et al., 2019, Section 4.3)
#
# A multinomial logit (binary logit for K = 2) with coefficients alpha, a
# (K - 1) x q matrix (row l: class l against the baseline class 0), on a
# design matrix built from a one-sided formula. Fits use nnet::multinom.
# ---------------------------------------------------------------------------

#' Design matrix of a prevalence model
#'
#' @param formula One-sided formula, e.g. \code{~ x1 + region}.
#' @param data Data frame.
#' @param terms,xlevels Terms and factor levels of an earlier design; when
#'   given, \code{data} is coded consistently with it.
#' @return List with \code{x} (model matrix), \code{terms}, \code{xlevels}.
#' @keywords internal
.mc_prev_design <- function(formula, data, terms = NULL, xlevels = NULL) {
  if (is.null(data))
    stop("prevalence_model needs the data containing its variables (use ",
         "the formula interface of mcglm(), the data argument of ",
         "estimate_mc(), and validation_sample(data = ) for an external ",
         "audit).", call. = FALSE)
  if (is.null(terms)) {
    if (length(formula) != 2L)
      stop("prevalence_model must be a one-sided formula, e.g. ~ x1 + region.",
           call. = FALSE)
    mf <- stats::model.frame(formula, data, na.action = stats::na.fail)
    terms <- stats::terms(mf)
    xlevels <- stats::.getXlevels(terms, mf)
  } else {
    mf <- stats::model.frame(terms, data, xlev = xlevels,
                             na.action = stats::na.fail)
  }
  list(x = stats::model.matrix(terms, mf), terms = terms, xlevels = xlevels)
}

#' Class probabilities of a prevalence model
#'
#' @param alpha (K - 1) x q coefficient matrix.
#' @param w n x q design matrix.
#' @return n x K matrix of \eqn{\Pr(Z = \ell \mid x_i)}.
#' @keywords internal
.mc_prev_probs <- function(alpha, w) {
  eta <- cbind(0, w %*% t(alpha))
  eta <- eta - eta[cbind(seq_len(nrow(eta)), max.col(eta, ties.method = "first"))]
  ex <- exp(eta)
  ex / rowSums(ex)
}

#' Weighted multinomial logit fit with nnet::multinom
#'
#' @param codes 0-based class codes.
#' @param w Design matrix (including any intercept column).
#' @param weights Non-negative weights.
#' @param K Number of classes.
#' @return (K - 1) x q coefficient matrix.
#' @keywords internal
.mc_multinom <- function(codes, w, weights, K) {
  keep <- weights > 0
  response <- factor(codes[keep], levels = seq_len(K) - 1L)
  design <- w[keep, , drop = FALSE]
  wts <- weights[keep]
  fit <- nnet::multinom(response ~ design - 1, weights = wts, trace = FALSE,
                        maxit = 1000L, reltol = 1e-12,
                        MaxNWts = max(1000L, 2L * K * ncol(w)))
  matrix(stats::coef(fit), K - 1L, ncol(w))
}

#' Likelihood-ratio test of independence between Z and covariates
#'
#' Compares multinomial logits of the true category on the non-constant
#' columns of \code{x} and on an intercept only, within the validation
#' sample. Used to warn when the constant-prevalence methods are at risk.
#' @return \code{NULL} when \code{x} has no non-constant column, otherwise
#'   a list with \code{statistic}, \code{df} and \code{p.value}.
#' @keywords internal
.mc_independence_test <- function(z, x, K) {
  keep <- vapply(seq_len(ncol(x)), function(k) length(unique(x[, k])) > 1L,
                 logical(1))
  if (!any(keep)) return(NULL)
  X <- x[, keep, drop = FALSE]
  response <- factor(z, levels = seq_len(K) - 1L)
  f1 <- nnet::multinom(response ~ X, trace = FALSE, maxit = 1000L)
  f0 <- nnet::multinom(response ~ 1, trace = FALSE, maxit = 1000L)
  stat <- max(f0$deviance - f1$deviance, 0)
  df <- (K - 1L) * ncol(X)
  list(statistic = stat, df = df,
       p.value = stats::pchisq(stat, df, lower.tail = FALSE))
}
