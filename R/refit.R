# =========================================================================
# refit generic and method
# =========================================================================

#' Refit a SIMEX model with a different extrapolation function
#'
#' Re-extrapolates the stored simulation results using a different fitting
#' method, without re-running the simulation step.
#'
#' @param object a \code{simex} object.
#' @param ... additional arguments, including \code{extrapolation}
#'   (\code{"quadratic"}, \code{"linear"}, or \code{"loglinear"}) and
#'   \code{jackknife} (logical or character).
#'
#' @return An updated \code{simex} object with new extrapolation results.
#'
#' @export
refit <- function(object, ...) UseMethod("refit")

#' @export
refit.simex <- function(object, extrapolation = "quadratic",
                        jackknife = TRUE, ...) {
  extrapolation <- match.arg(extrapolation,
                             c("quadratic", "linear", "loglinear"))

  # Improved MC-SIMEX uses exact correction — nothing to re-extrapolate

  if (!is.null(object$method) && object$method == "improved")
    stop("refit() is not applicable for improved MC-SIMEX ",
         "(exact correction, no extrapolation to change).", call. = FALSE)

  lambda <- object$lambda
  estimates <- object$SIMEX.estimates

  # Re-extrapolate
  extrap_result <- extrapolate_simex(lambda, estimates, extrapolation)
  object$coefficients <- extrap_result$coefficients
  object$extrapolation <- extrap_result$model
  object$extrapolation.method <- extrapolation

  # Re-estimate variance
  if (!isFALSE(jackknife) && !is.null(object$theta)) {
    jk_method <- if (is.character(jackknife)) jackknife else extrapolation

    if (object$error.type == "me") {
      vcov_jk <- jackknife_variance_simex(
        object$theta, object$naive.model, lambda[-1], jk_method,
        object$vcov.model
      )
    } else {
      # Standard MC-SIMEX — build xi_hat for variance estimation
      fam <- get_link_funs(object$family)
      y <- object$naive.model$y
      wt <- if (!is.null(object$naive.model$prior.weights))
        object$naive.model$prior.weights else rep(1.0, object$n)
      xi_hat <- .refit_build_xi_hat(object)

      vcov_jk <- jackknife_variance_mcsimex(
        object$theta, object$naive.coefficients, object$naive.model,
        lambda[-1], jk_method, fam, xi_hat, y, wt
      )
    }
    p_names <- names(object$coefficients)
    rownames(vcov_jk) <- colnames(vcov_jk) <- p_names
    object$vcov <- vcov_jk
  }

  # Update fitted values and residuals
  if (object$error.type == "me") {
    X <- object$naive.model$x
    eta <- as.numeric(X %*% object$coefficients)
  } else {
    xi_hat <- .refit_build_xi_hat(object)
    eta <- as.numeric(xi_hat %*% object$coefficients)
  }
  object$fitted.values <- object$family$linkinv(eta)
  object$residuals <- object$naive.model$y - object$fitted.values

  object$call$extrapolation <- extrapolation
  object
}


#' Build xi_hat for MC-SIMEX refit, handling single-mc, multi-mc, and response-mc
#' @keywords internal
.refit_build_xi_hat <- function(object) {
  # Use stored xi.hat if available (always set since v0.3.0)
  if (!is.null(object$xi.hat))
    return(object$xi.hat)

  # Fallback for response-only mc (xi_hat = naive model's X)
  object$naive.model$x
}
