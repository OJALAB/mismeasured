# ---------------------------------------------------------------------------
# Propagating estimated misclassification probabilities into the variance of
# the plug-in estimators: "cs", "bca", "bcm" and "onestep" (fix_omega).
#
# These estimators take the probabilities as plug-ins; their point
# estimates do not change when the probabilities come from a validation
# sample (no true-category rows). Their covariance adds the estimation
# uncertainty of the probabilities through the stacked estimating
# functions of estimate_mc() -- column-conditional Pi and the prevalence --
# so any validation design (weights, strata, Hajek/HT, prevalence = "em")
# is covered:
#   cs, onestep: Z-estimators, stacked with eta at the fitted psi;
#   bca, bcm:    written as estimating equations in (psi_naive, psi, eta)
#                and stacked with the naive score and eta, so that the
#                averages they use (information, drift, score) contribute.
# ---------------------------------------------------------------------------

#' Check supplied probabilities against a validation-sample estimate
#'
#' Plug-in methods accept probabilities together with a validation sample
#' only when they equal the sample's estimates (the variance assumes so).
#' @keywords internal
.mc_check_supplied <- function(est, K, Pi = NULL, pi_z = NULL, p01 = NULL,
                               p10 = NULL, c1 = NULL, c2 = NULL) {
  P <- unname(est$Pi)
  pr <- unname(est$pi)
  ok <- function(a, b) is.null(a) ||
    (length(a) == length(b) && max(abs(as.numeric(a) - b)) <= 1e-8)
  good <- ok(Pi, P) &&
    ok(pi_z, if (K == 2L && length(pi_z) == 1L) pr[2L] else pr) &&
    ok(p01, if (K == 2L) P[2L, 1L] else NA) &&
    ok(p10, if (K == 2L) P[1L, 2L] else NA) &&
    ok(c1, if (K == 2L) P[2L, 1L] * pr[1L] else NA) &&
    ok(c2, if (K == 2L) P[2L, 1L] * pr[1L] - P[1L, 2L] * pr[2L] else NA)
  if (!good)
    stop("Supplied misclassification probabilities differ from the ",
         "empirical estimates of the validation sample; omit them to use ",
         "the estimates.", call. = FALSE)
  invisible(TRUE)
}

#' Stacked covariance of "cs" with estimated probabilities
#'
#' The drift-corrected score
#' \eqn{\hat\xi_i\{Y_i - \mu(\psi^\top\hat\xi_i)\} - m_i(\psi; \Pi, \pi)}
#' (the multicategory drift, which reduces to the \eqn{c_1, c_2} form for
#' \eqn{K = 2}) is stacked with the estimating functions of the estimate
#' object at the fitted \eqn{\hat\psi}.
#' @keywords internal
.mcglm_cs_validated_vcov <- function(psi, y, xi_hat, z_hat, x, K, family,
                                     est, wt = NULL,
                                     control = control_mc()) {
  fam <- .normalize_family(family)
  U_fun <- function(b, par)
    xi_hat * (y - fam$linkinv(drop(xi_hat %*% b))) -
      .mcglm_compute_m_multi(b, x, K, fam$linkinv, par$Pi, par$pi)
  J_fun <- function(b, par, w)
    .mc_num_jacobian(function(bb) colSums(w * U_fun(bb, par)), b)
  .mcglm_fit_validated(psi, est, y, x, z_hat, K, fam, U_fun, J_fun,
                       wt = wt, control = control, label = "CS",
                       psi_fixed = psi, beta_equation = "none")$vcov
}

#' Stacked covariance of "bca" / "bcm" with estimated probabilities
#'
#' Both are written as estimating equations in \eqn{(\psi_n, \psi, \eta)},
#' with \eqn{\psi_n} the naive estimate (naive score rows
#' \eqn{S_i(\psi_n)}), so that the sampling variability of every average
#' they use (information, drift, score) enters the sandwich:
#' \itemize{
#'   \item BCA: \eqn{\dot\mu_i \hat\xi_i\hat\xi_i^\top(\psi - \psi_n) +
#'     m_i(\psi_n)} (iterated: \eqn{m_i(\psi)});
#'   \item BCM: \eqn{(\dot\mu_i \hat\xi_i\hat\xi_i^\top + M_i)(\psi - \psi_n)
#'     - \{S_i(\psi_n) - m_i(\psi_n)\}}, with \eqn{M_i = \partial m_i /
#'     \partial\psi^\top} at \eqn{\psi_n} (iterated: the corrected score
#'     \eqn{S_i(\psi) - m_i(\psi)});
#' }
#' each summing to zero at the fitted values. They are stacked with the
#' estimating functions of the estimate object.
#'
#' @param type \code{"bca"} or \code{"bcm"}.
#' @keywords internal
.mcglm_bc_validated_vcov <- function(type, psi_naive, psi, y, xi_hat, z_hat,
                                     x, K, family, est, wt = NULL,
                                     iterate = FALSE) {
  fam <- .normalize_family(family)
  vb <- est$validation
  n <- length(y)
  p <- length(psi)
  q <- length(est$eta)
  w <- if (is.null(wt)) rep(1, n) else wt

  drift <- function(b, par)
    .mcglm_compute_m_multi(b, x, K, fam$linkinv, par$Pi, par$pi)
  score <- function(b) xi_hat * (y - fam$linkinv(drop(xi_hat %*% b)))
  info_times <- function(b, d)
    xi_hat * (fam$mu.eta(drop(xi_hat %*% b)) * drop(xi_hat %*% d))
  rows <- function(theta) {
    bn <- theta[seq_len(p)]
    b  <- theta[p + seq_len(p)]
    par <- est$map(theta[2L * p + seq_len(q)])
    d <- b - bn
    h <- if (type == "bca") {
      info_times(bn, d) + drift(if (iterate) b else bn, par)
    } else if (iterate) {
      score(b) - drift(b, par)
    } else {
      # M_i d by a central difference of the per-row drift along d
      e <- 1e-6
      Md <- (drift(bn + e * d, par) - drift(bn - e * d, par)) / (2 * e)
      info_times(bn, d) + Md - (score(bn) - drift(bn, par))
    }
    cbind(score(bn), h)
  }
  r <- est$rows(est$eta)
  total <- function(theta) {
    out <- c(colSums(w * rows(theta)), numeric(q))
    rr <- est$rows(theta[2L * p + seq_len(q)])
    out[2L * p + seq_len(q)] <- colSums(vb$w_v * rr$v) +
      (if (is.null(rr$m)) 0 else colSums(w * rr$m))
    out
  }
  theta <- c(psi_naive, psi, unname(est$eta))
  A <- .mc_num_jacobian(total, theta)
  G_m <- cbind(rows(theta), if (is.null(r$m)) matrix(0, n, q) else r$m)
  G_v <- cbind(matrix(0, vb$n, 2L * p), r$v)
  V <- .mc_sandwich_or_na(A, .mc_meat(G_m, G_v, vb, w, n), toupper(type))
  V[p + seq_len(p), p + seq_len(p), drop = FALSE]
}

#' Stacked covariance of "onestep" (fix_omega = TRUE) with estimated
#' probabilities
#'
#' The one-step estimate maximises the mixture likelihood with weights
#' \eqn{\Pi_{\hat z \ell}\pi_\ell} plugged in; its score is the
#' posterior-weighted score of \code{.mcglm_ec_parts}, stacked with the
#' estimate object at the fitted \eqn{\hat\psi} (and \eqn{\hat\sigma} for
#' gaussian, re-solved from its score equation).
#' @keywords internal
.mcglm_onestep_validated_vcov <- function(psi, y, xi_hat, z_hat, x, K,
                                          family, est, wt = NULL,
                                          control = control_mc()) {
  model <- .mcglm_ec_model(y, x, K, family)
  w <- if (is.null(wt)) rep(1, length(y)) else wt
  par0 <- .mc_with_rows(est$map(est$eta), z_hat, est$w_main)
  theta <- psi
  if (model$has_sigma) {
    tau_score <- function(tau)
      sum(w * .mcglm_ec_parts(c(psi, tau), model, log(par0$Wi))$U[, model$p + 1L])
    res <- drop(y - xi_hat %*% psi)
    t0 <- 0.5 * log(sum(w * res^2) / sum(w))
    theta <- c(psi, stats::uniroot(tau_score, c(t0 - 3, t0 + 3),
                                   extendInt = "yes", tol = 1e-12)$root)
  }
  U_fun <- function(th, par) .mcglm_ec_parts(th, model, log(par$Wi))$U
  J_fun <- function(th, par, ww)
    .mcglm_ec_parts(th, model, log(par$Wi), w = ww, jac = TRUE)$J
  V <- .mcglm_fit_validated(theta, est, y, x, z_hat, K, model$family, U_fun,
                            J_fun, wt = wt, control = control,
                            label = "onestep", psi_fixed = theta,
                            beta_equation = "none")$vcov
  V[seq_len(model$p), seq_len(model$p), drop = FALSE]
}
