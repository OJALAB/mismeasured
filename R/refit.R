# =========================================================================
# refit generic and methods
# =========================================================================

#' Refit a SIMEX or MC-SIMEX model with a different extrapolation function
#'
#' Re-extrapolates the stored simulation results using a different fitting
#' method, without re-running the simulation step. Pass \code{extrapolation}
#' and \code{jackknife} arguments to the method.
#'
#' @param object a \code{simex} or \code{mcsimex} object.
#' @param ... additional arguments passed to methods, including
#'   \code{extrapolation} (\code{"quadratic"}, \code{"linear"}, or
#'   \code{"loglinear"}) and \code{jackknife} (logical or character).
#'
#' @return An updated object of the same class with new extrapolation results.
#'
#' @examples
#' \dontrun{
#' fit_quad <- simex(y ~ x, data = df, SIMEXvariable = "x",
#'                   measurement.error = 0.3, extrapolation = "quadratic")
#' fit_lin <- refit(fit_quad, extrapolation = "linear")
#' }
#'
#' @export
refit <- function(object, ...) UseMethod("refit")


#' @export
refit.simex <- function(object, extrapolation = "quadratic",
                        jackknife = TRUE, ...) {
  extrapolation <- match.arg(extrapolation,
                             c("quadratic", "linear", "loglinear"))

  lambda <- object$lambda
  estimates <- object$SIMEX.estimates

  # Re-extrapolate
  extrap_result <- extrapolate_simex(lambda, estimates, extrapolation)
  object$coefficients <- extrap_result$coefficients
  object$extrapolation <- extrap_result$model
  object$extrapolation.method <- extrapolation

  # Re-estimate variance if requested
  if (!isFALSE(jackknife) && !is.null(object$theta)) {
    jk_method <- if (is.character(jackknife)) jackknife else extrapolation
    vcov_jk <- jackknife_variance_simex(
      object$theta, object$naive.model, lambda[-1], jk_method
    )
    p_names <- names(object$coefficients)
    rownames(vcov_jk) <- colnames(vcov_jk) <- p_names
    object$vcov <- vcov_jk
  }

  # Update fitted values and residuals
  X <- object$naive.model$x
  family <- object$family
  eta <- as.numeric(X %*% object$coefficients)
  object$fitted.values <- family$linkinv(eta)
  object$residuals <- object$naive.model$y - object$fitted.values

  object$call$extrapolation <- extrapolation
  object
}


#' @export
refit.mcsimex <- function(object, extrapolation = "quadratic",
                          jackknife = TRUE, ...) {
  extrapolation <- match.arg(extrapolation,
                             c("quadratic", "linear", "loglinear"))

  lambda <- object$lambda
  estimates <- object$SIMEX.estimates

  # Re-extrapolate
  extrap_result <- extrapolate_simex(lambda, estimates, extrapolation)
  object$coefficients <- extrap_result$coefficients
  object$extrapolation <- extrap_result$model
  object$extrapolation.method <- extrapolation

  # Re-estimate variance if requested
  if (!isFALSE(jackknife) && !is.null(object$theta)) {
    jk_method <- if (is.character(jackknife)) jackknife else extrapolation
    K <- nrow(object$mc.matrix[[object$SIMEXvariable]])
    fam <- get_link_funs(object$family)

    # Rebuild design matrix components
    data <- object$data
    z_hat <- as.integer(data[[object$SIMEXvariable]]) - 1L
    other_terms <- setdiff(all.vars(object$formula)[-1], object$SIMEXvariable)
    if (length(other_terms) > 0) {
      x_mat <- model.matrix(reformulate(other_terms, intercept = TRUE),
                            data = data)
    } else {
      x_mat <- matrix(1, nrow = object$n, ncol = 1)
    }
    xi_hat <- .build_xi_hat(z_hat, x_mat, K)
    y <- object$naive.model$y
    wt <- if (!is.null(object$naive.model$prior.weights))
      object$naive.model$prior.weights else rep(1.0, object$n)

    vcov_jk <- jackknife_variance_mcsimex(
      object$theta, object$naive.coefficients, object$naive.model,
      lambda[-1], jk_method, fam, xi_hat, y, wt
    )
    p_names <- names(object$coefficients)
    rownames(vcov_jk) <- colnames(vcov_jk) <- p_names
    object$vcov <- vcov_jk
  }

  # Update fitted values and residuals
  data <- object$data
  z_hat <- as.integer(data[[object$SIMEXvariable]]) - 1L
  K <- nrow(object$mc.matrix[[object$SIMEXvariable]])
  other_terms <- setdiff(all.vars(object$formula)[-1], object$SIMEXvariable)
  if (length(other_terms) > 0) {
    x_mat <- model.matrix(reformulate(other_terms, intercept = TRUE),
                          data = data)
  } else {
    x_mat <- matrix(1, nrow = object$n, ncol = 1)
  }
  xi_hat <- .build_xi_hat(z_hat, x_mat, K)
  eta <- as.numeric(xi_hat %*% object$coefficients)
  object$fitted.values <- object$family$linkinv(eta)
  object$residuals <- object$naive.model$y - object$fitted.values

  object$call$extrapolation <- extrapolation
  object
}
