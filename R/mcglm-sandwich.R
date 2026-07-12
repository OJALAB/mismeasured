# ---------------------------------------------------------------------------
# sandwich-package integration: estfun() and bread() methods for mcglm
# Registered lazily on the sandwich generics (sandwich in Suggests).
# ---------------------------------------------------------------------------

#' Extract empirical estimating functions from an \code{mcglm} fit
#'
#' Returns the \eqn{n \times p} matrix of per-observation estimating
#' functions evaluated at the coefficients of the chosen \code{method}.
#' For \code{method = "cs"} these are the corrected-score contributions
#' \eqn{\phi_i(\hat\psi) = \hat\xi_i\{y_i - \mu(\tilde\eta_i)\} - m_i(\hat\psi)};
#' for \code{"naive"}, \code{"bca"} and \code{"bcm"} they are the proxy-score
#' contributions \eqn{\hat\xi_i\{y_i - \mu(\tilde\eta_i)\}} (whose sandwich is
#' the asymptotic variance of all three estimators in the drifting regime).
#'
#' Together with \code{\link{bread.mcglm}} this makes \code{mcglm} fits work
#' with the \pkg{sandwich} ecosystem (\code{meat()},
#' \code{lmtest::coeftest(vcov. = )}, bootstrap tools, ...).
#'
#' @section Assembling the sandwich:
#' Prefer \code{vcov(fit, method = )}: it already \emph{is} the correct
#' sandwich \eqn{J^{-1} S J^{-\top}/n}. Do not pass \code{method} to
#' \code{sandwich::sandwich()} -- it forwards extra arguments to the meat
#' but \emph{not} to the bread, silently mixing estimators -- and note
#' that \code{sandwich::sandwich()} assembles \code{bread \%*\% meat \%*\%
#' bread} without a transpose, which is only valid for the symmetric
#' bread of \code{"naive"}/\code{"bca"}/\code{"bcm"}; the corrected-score
#' bread \eqn{(\hat I + \hat M)^{-1}} is not symmetric. To assemble by
#' hand: \code{B \%*\% meat(fit, method = m) \%*\% t(B) / n} with
#' \code{B = bread(fit, method = m)}, which reproduces
#' \code{vcov(fit, method = m)} exactly for unweighted fits.
#'
#' @param x An \code{mcglm} object.
#' @param method Estimation method; defaults to the last fitted one.
#'   One of \code{"naive"}, \code{"bca"}, \code{"bcm"}, \code{"cs"}.
#' @param ... Unused.
#' @return An \eqn{n \times p} numeric matrix. With frequency weights the
#'   rows are multiplied by the weights (the \pkg{sandwich} convention for
#'   \code{glm}), in which case \code{sandwich::sandwich()} treats them as
#'   probability weights and differs from \code{vcov()}, which treats them
#'   as frequencies.
#' @seealso \code{\link{bread.mcglm}}, \code{\link{vcov.mcglm}}
#' @exportS3Method sandwich::estfun
estfun.mcglm <- function(x, method = NULL, ...) {
  if (isTRUE(x$is_multinomial))
    stop("estfun() is not available for multinomial mcglm fits.")
  if (is.null(method)) method <- .mcglm_default_method(x)
  method <- match.arg(method, c("naive", "bca", "bcm", "cs"))
  psi <- unname(stats::coef(x, method = method))
  fam <- .mcglm_get_link_funs(x$family)

  eta   <- .mcglm_linear_predictor(x, psi)
  resid <- x$y - fam$mu(eta)
  ef    <- x$xi_hat * resid

  if (method == "cs") {
    mc <- x$misclass
    if (x$K == 2L) {
      ef <- ef - .mcglm_compute_m_bin(psi, x$x, fam$mu, mc$c1, mc$c2)
    } else {
      ef <- ef - .mcglm_compute_m_multi(psi, x$x, x$K, fam$mu, mc$Pi, mc$pi_z)
    }
  }

  if (!is.null(x$weights)) ef <- x$weights * ef
  colnames(ef) <- names(stats::coef(x, method = method))
  ef
}

#' Bread matrix for an \code{mcglm} fit
#'
#' Returns the inverse of the (negated, averaged) Jacobian of the estimating
#' functions in \code{\link{estfun.mcglm}}: \eqn{(\hat I + \hat M)^{-1}} for
#' \code{method = "cs"} and \eqn{\hat I^{-1}} for \code{"naive"},
#' \code{"bca"}, \code{"bcm"}, each evaluated at that method's coefficients.
#'
#' @param x An \code{mcglm} object.
#' @param method Estimation method; defaults to the last fitted one.
#' @param ... Unused.
#' @return A \eqn{p \times p} matrix, such that
#'   \code{sandwich::sandwich(fit)} assembles
#'   \eqn{n^{-1} B \,\hat S\, B^\top} with \eqn{\hat S} from
#'   \code{estfun.mcglm}.
#' @seealso \code{\link{estfun.mcglm}}, \code{\link{vcov.mcglm}}
#' @exportS3Method sandwich::bread
bread.mcglm <- function(x, method = NULL, ...) {
  if (isTRUE(x$is_multinomial))
    stop("bread() is not available for multinomial mcglm fits.")
  if (is.null(method)) method <- .mcglm_default_method(x)
  method <- match.arg(method, c("naive", "bca", "bcm", "cs"))
  psi <- unname(stats::coef(x, method = method))
  fam <- .mcglm_get_link_funs(x$family)
  wt  <- x$weights

  if (x$K == 2L) {
    J <- .mcglm_compute_Ihat(psi, x$xi_hat, fam$mu_dot, wt = wt)
    if (method == "cs") {
      mc <- x$misclass
      J <- J + .mcglm_compute_Mhat_bin(psi, x$x, fam$mu, fam$mu_dot,
                                       mc$c1, mc$c2, wt = wt)
    }
  } else {
    J <- .mcglm_compute_Ihat_multi(psi, x$xi_hat, x$z_hat, x$K, fam$mu_dot,
                                   wt = wt)
    if (method == "cs") {
      mc <- x$misclass
      J <- J + .mcglm_compute_Mhat_multi(psi, x$x, x$K, fam$mu, mc$Pi,
                                         mc$pi_z, wt = wt,
                                         mu_dot_fun = fam$mu_dot)
    }
  }
  # the internal helpers average with weight-sum N; rescale to the
  # per-observation (1/n) convention that sandwich::sandwich() assumes
  if (!is.null(wt)) J <- J * (sum(wt) / x$n)

  B <- solve(J)
  nms <- names(stats::coef(x, method = method))
  dimnames(B) <- list(nms, nms)
  B
}

#' Linear predictor at arbitrary coefficients (binary or multicategory:
#' the dummy encoding in xi_hat picks out the z-specific gamma, so a
#' single matrix product covers both cases)
#' @keywords internal
.mcglm_linear_predictor <- function(object, psi) {
  as.numeric(object$xi_hat %*% psi)
}

#' Parse an \code{mc()} model formula
#'
#' Public interface to the formula parsing used by \code{\link{mcglm}}:
#' splits \code{y ~ mc(z, Pi) + x1 + ...} into the response, the 0-based
#' proxy covariate, the model matrix of the remaining terms, and the
#' misclassification matrix. Intended for packages building on
#' \pkg{mismeasured} (e.g. survey-integration wrappers) that must
#' intercept the \code{mc()} term before standard formula machinery
#' (such as \code{\link[stats]{model.matrix}}) sees it.
#'
#' @param formula A formula with exactly one \code{\link{mc}} term on the
#'   right-hand side, e.g. \code{y ~ mc(z, Pi) + x1}.
#' @param data A data frame containing the variables.
#' @param env Environment in which to evaluate the \code{mc()} matrix
#'   argument (default: the caller's environment).
#' @return A list with components \code{y} (numeric response),
#'   \code{z_hat} (integer proxy, coded \code{0, ..., K-1}),
#'   \code{x} (model matrix of the remaining terms, including an
#'   intercept), \code{Pi} (the \eqn{K \times K} misclassification
#'   matrix, or \code{NULL} if not supplied inside \code{mc()}), and
#'   \code{K} (number of categories).
#'
#' @examples
#' Pi <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
#' df <- data.frame(y = rpois(50, 1), z = rbinom(50, 1, 0.4),
#'                  x1 = rnorm(50))
#' parts <- mc_parse_formula(y ~ mc(z, Pi) + x1, df)
#' str(parts)
#'
#' @seealso \code{\link{mc}}, \code{\link{mcglm}}
#' @export
mc_parse_formula <- function(formula, data, env = parent.frame()) {
  .mcglm_parse_formula(formula, data, env)
}
