# =========================================================================
# SIMEX: Simulation-Extrapolation for continuous measurement error
# =========================================================================

#' SIMEX estimator for GLMs with measurement error
#'
#' Implements the Simulation-Extrapolation algorithm of Cook and Stefanski (1994)
#' for correcting attenuation bias due to measurement error in covariates.
#' The simulation step runs in C++ for high performance.
#'
#' @param formula a formula describing the model (as in \code{\link{glm}}).
#' @param family a family object (e.g. \code{gaussian()}, \code{binomial()}).
#'   Default: \code{gaussian()}.
#' @param data a data frame containing the variables in the formula.
#' @param SIMEXvariable character vector naming the variable(s) measured with
#'   error. These must appear as untransformed main effects in the formula.
#' @param measurement.error known measurement error standard deviations.
#'   Can be: a scalar (homoscedastic, one variable), a numeric vector of length
#'   \code{length(SIMEXvariable)} (homoscedastic, multiple variables), or a
#'   matrix with \code{nrow(data)} rows and \code{length(SIMEXvariable)} columns
#'   (heteroscedastic).
#' @param lambda numeric vector of SIMEX exponents (default:
#'   \code{c(0.5, 1, 1.5, 2)}).
#' @param B number of simulation replicates per lambda level (default: 200).
#' @param extrapolation extrapolation method: \code{"quadratic"} (default),
#'   \code{"linear"}, or \code{"loglinear"}.
#' @param jackknife logical or character. If \code{TRUE}, compute jackknife
#'   variance using the same extrapolation method. If a character string, use
#'   that method. \code{FALSE} to skip.
#' @param weights optional prior weights (as in \code{\link{glm}}).
#' @param seed random seed for reproducibility (default: 42).
#'
#' @return An object of class \code{"simex"} with components:
#'   \describe{
#'     \item{coefficients}{Named vector of SIMEX-corrected coefficients.}
#'     \item{naive.coefficients}{Coefficients from the naive (uncorrected) model.}
#'     \item{residuals}{Working residuals.}
#'     \item{fitted.values}{Fitted values from the corrected model.}
#'     \item{family}{The family object.}
#'     \item{formula}{The model formula.}
#'     \item{call}{The matched call.}
#'     \item{data}{The data used.}
#'     \item{vcov}{Jackknife variance-covariance matrix (if requested).}
#'     \item{SIMEX.estimates}{Matrix of averaged estimates at each lambda.}
#'     \item{theta}{List of B x p matrices of simulation estimates.}
#'     \item{extrapolation}{Fitted extrapolation model(s).}
#'     \item{lambda}{Lambda vector used (including 0).}
#'     \item{B}{Number of replicates.}
#'     \item{SIMEXvariable}{Name(s) of the error-prone variable(s).}
#'     \item{measurement.error}{The measurement error specification.}
#'     \item{naive.model}{The naive \code{glm} fit.}
#'     \item{extrapolation.method}{Extrapolation method used.}
#'   }
#'
#' @details
#' The SIMEX algorithm has two steps:
#'
#' \strong{Simulation step}: For each \eqn{\lambda} in the lambda grid and
#' \eqn{b = 1, \ldots, B}, extra measurement error is added:
#' \eqn{X_b(\lambda) = X_{measured} + \sqrt{\lambda} \cdot \varepsilon_b \cdot \sigma_{me}},
#' where \eqn{\varepsilon_b \sim N(0,1)}. The model is refitted to each
#' contaminated dataset. The entire loop runs in C++.
#'
#' \strong{Extrapolation step}: The averaged coefficient estimates are modelled
#' as a function of lambda and extrapolated to \eqn{\lambda = -1}.
#'
#' @references
#' Cook, J.R. and Stefanski, L.A. (1994). Simulation-extrapolation estimation
#' in parametric measurement error models. \emph{Journal of the American
#' Statistical Association}, 89, 1314--1328.
#'
#' Carroll, R.J., Ruppert, D., Stefanski, L.A. and Crainiceanu, C. (2006).
#' \emph{Measurement error in nonlinear models: A modern perspective.}
#' Second Edition. London: Chapman and Hall.
#'
#' @examples
#' # Simulate data with measurement error
#' set.seed(42)
#' n <- 500
#' x_true <- rnorm(n)
#' y <- 1 + 2 * x_true + rnorm(n, sd = 0.5)
#' x_measured <- x_true + rnorm(n, sd = 0.5)
#'
#' df <- data.frame(y = y, x = x_measured)
#' fit <- simex(y ~ x, family = gaussian(), data = df,
#'              SIMEXvariable = "x", measurement.error = 0.5, B = 100)
#' summary(fit)
#'
#' @export
simex <- function(formula, family = gaussian(), data,
                  SIMEXvariable, measurement.error,
                  lambda = c(0.5, 1, 1.5, 2),
                  B = 200L,
                  extrapolation = c("quadratic", "linear", "loglinear"),
                  jackknife = TRUE,
                  weights = NULL,
                  seed = 42L) {

  cl <- match.call()
  extrapolation <- match.arg(extrapolation)

  # --- resolve family ---
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())()
  if (is.function(family)) family <- family()
  if (is.null(family$family))
    stop("'family' not recognized")

  dist_code <- family_to_dist_code(family)

  # --- validate inputs ---
  if (!is.character(SIMEXvariable))
    stop("SIMEXvariable must be a character vector")

  for (sv in SIMEXvariable) {
    if (!sv %in% names(data))
      stop("SIMEXvariable '", sv, "' not found in data")
    if (is.factor(data[[sv]]))
      stop("SIMEXvariable '", sv, "' is a factor. Use mcsimex() for discrete misclassification.")
  }

  lambda <- sort(as.numeric(lambda))
  if (any(lambda <= 0))
    stop("lambda values must be positive")
  B <- as.integer(B)
  n <- nrow(data)
  n_simex <- length(SIMEXvariable)

  # --- process measurement.error ---
  measurement.error <- as.matrix(measurement.error)
  if (nrow(measurement.error) == 1 && n > 1) {
    measurement.error <- matrix(measurement.error, nrow = n,
                                ncol = ncol(measurement.error), byrow = TRUE)
  }
  if (ncol(measurement.error) == 1 && n_simex > 1) {
    measurement.error <- matrix(measurement.error, nrow = n, ncol = n_simex)
  }
  if (nrow(measurement.error) != n || ncol(measurement.error) != n_simex)
    stop("measurement.error dimensions must match data rows x SIMEXvariable count")
  if (any(measurement.error < 0))
    stop("measurement.error must be non-negative")
  if (any(colSums(measurement.error) == 0))
    stop("measurement.error must not be zero for all observations")

  # --- fit naive model ---
  if (!is.null(weights)) {
    naive_fit <- glm(formula, family = family, data = data, weights = weights,
                     x = TRUE, y = TRUE)
  } else {
    naive_fit <- glm(formula, family = family, data = data,
                     x = TRUE, y = TRUE)
  }

  y <- naive_fit$y
  X <- naive_fit$x  # model matrix
  p <- ncol(X)
  p_names <- colnames(X)
  psi_naive <- unname(coef(naive_fit))
  names(psi_naive) <- p_names

  # Weights
  wt <- if (is.null(weights)) rep(1.0, n) else as.numeric(weights)

  # --- identify SIMEX columns in model matrix ---
  simex_cols <- integer(n_simex)
  for (j in seq_len(n_simex)) {
    idx <- which(p_names == SIMEXvariable[j])
    if (length(idx) != 1)
      stop("SIMEXvariable '", SIMEXvariable[j],
           "' must appear as exactly one untransformed main effect in the model matrix. ",
           "Interactions and transformations are not supported in the C++ fast path.")
    simex_cols[j] <- idx - 1L  # 0-based for C++
  }

  # --- simulation step (C++) ---
  sim_result <- simex_sim_cpp(y, X, simex_cols, measurement.error,
                              dist_code, lambda, B, wt, as.integer(seed))

  n_lambda <- length(lambda)

  # Reshape: list of B x p matrices
  theta_list <- vector("list", n_lambda)
  avg_estimates <- matrix(0, n_lambda, p)
  for (l in seq_len(n_lambda)) {
    rows <- ((l - 1) * B + 1):(l * B)
    theta_mat <- sim_result[rows, , drop = FALSE]
    colnames(theta_mat) <- p_names
    theta_list[[l]] <- theta_mat
    avg_estimates[l, ] <- colMeans(theta_mat)
  }
  names(theta_list) <- paste0("lambda_", lambda)

  # Stack estimates
  all_lambda <- c(0, lambda)
  all_estimates <- rbind(psi_naive, avg_estimates)
  rownames(all_estimates) <- paste0("lambda_", all_lambda)
  colnames(all_estimates) <- p_names

  # --- extrapolation step ---
  extrap_result <- extrapolate_simex(all_lambda, all_estimates, extrapolation)

  # --- jackknife variance ---
  vcov_jk <- NULL
  if (!isFALSE(jackknife)) {
    jk_method <- if (is.character(jackknife)) jackknife else extrapolation
    vcov_jk <- jackknife_variance_simex(theta_list, naive_fit, lambda, jk_method)
    rownames(vcov_jk) <- colnames(vcov_jk) <- p_names
  }

  # --- fitted values and residuals ---
  corrected_coefs <- extrap_result$coefficients
  eta_corrected <- as.numeric(X %*% corrected_coefs)
  fitted_vals <- family$linkinv(eta_corrected)
  resids <- y - fitted_vals

  # --- build output ---
  out <- structure(
    list(
      coefficients       = corrected_coefs,
      naive.coefficients = psi_naive,
      residuals          = resids,
      fitted.values      = fitted_vals,
      family             = family,
      formula            = formula,
      call               = cl,
      data               = data,
      vcov               = vcov_jk,
      SIMEX.estimates    = all_estimates,
      theta              = theta_list,
      extrapolation      = extrap_result$model,
      lambda             = all_lambda,
      B                  = B,
      n                  = n,
      p                  = p,
      SIMEXvariable      = SIMEXvariable,
      measurement.error  = measurement.error,
      naive.model        = naive_fit,
      extrapolation.method = extrapolation
    ),
    class = "simex"
  )
  out
}


# =========================================================================
# S3 methods for class "simex"
# =========================================================================

#' @export
print.simex <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\nSIMEX corrected coefficients:\n")
  print.default(format(coef(x), digits = digits), print.gap = 2,
                quote = FALSE)
  cat("\n")
  invisible(x)
}


#' @export
summary.simex <- function(object, ...) {
  est <- coef(object)
  p.names <- names(est)
  n <- object$n
  p <- length(est)
  rdf <- n - p

  ans <- list()
  ans$call <- object$call
  ans$family <- object$family
  ans$formula <- object$formula
  ans$SIMEXvariable <- object$SIMEXvariable
  ans$measurement.error <- object$measurement.error
  ans$lambda <- object$lambda
  ans$B <- object$B
  ans$extrapolation.method <- object$extrapolation.method
  ans$SIMEX.estimates <- object$SIMEX.estimates
  ans$naive.coefficients <- object$naive.coefficients
  ans$residuals <- object$residuals
  ans$n <- n
  ans$df.residual <- rdf

  if (!is.null(object$vcov)) {
    se <- sqrt(pmax(diag(object$vcov), 0))
    tval <- est / se
    pval <- 2 * pt(abs(tval), rdf, lower.tail = FALSE)
    ans$coefficients <- cbind(
      Estimate = est, `Std. Error` = se,
      `t value` = tval, `Pr(>|t|)` = pval
    )
    rownames(ans$coefficients) <- p.names
  } else {
    ans$coefficients <- cbind(Estimate = est)
    rownames(ans$coefficients) <- p.names
  }

  class(ans) <- "summary.simex"
  ans
}


#' @export
print.summary.simex <- function(x,
                                digits = max(3, getOption("digits") - 3),
                                ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\nFamily:", x$family$family, "\n")
  cat("SIMEX variable(s):", paste(x$SIMEXvariable, collapse = ", "), "\n")

  # Summarize measurement error
  me_summary <- apply(x$measurement.error, 2, function(col) {
    if (length(unique(col)) == 1) unique(col) else "heteroscedastic"
  })
  cat("Measurement error:", paste(me_summary, collapse = ", "), "\n")
  cat("Extrapolation:", x$extrapolation.method, "\n")
  cat("Lambda grid:", paste(x$lambda, collapse = ", "), "\n")
  cat("B =", x$B, ", n =", x$n, "\n")

  if (!is.null(x$residuals)) {
    cat("\nResiduals:\n")
    rq <- quantile(x$residuals, na.rm = TRUE)
    names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
    print(rq, digits = digits)
  }

  cat("\nNaive coefficients:\n")
  print(round(x$naive.coefficients, digits))

  cat("\nSIMEX corrected coefficients:\n")
  if (ncol(x$coefficients) > 1) {
    printCoefmat(x$coefficients, digits = digits)
  } else {
    print(round(x$coefficients, digits))
  }
  cat("\n")
  invisible(x)
}


#' @export
plot.simex <- function(x,
                       xlab = expression(1 + lambda),
                       ylab = NULL,
                       ask = FALSE,
                       show = NULL, ...) {
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  if (ask) par(ask = TRUE)

  b <- x$SIMEX.estimates
  p <- ncol(b)
  p.names <- colnames(b)
  lam <- x$lambda

  if (is.null(ylab)) ylab <- p.names
  if (is.null(show)) show <- rep(TRUE, p)

  a <- seq(-1, max(lam), by = 0.01)

  for (j in seq_len(p)) {
    if (!show[j]) next

    fit_j <- x$extrapolation[[j]]
    nd <- data.frame(lambda = a)
    if (x$extrapolation.method == "loglinear") {
      offset <- attr(fit_j, "loglinear_offset")
      if (is.null(offset)) offset <- 0
      d <- exp(predict(fit_j, newdata = nd)) - offset
    } else {
      d <- predict(fit_j, newdata = nd)
    }

    plot(lam + 1, b[, j], xlab = xlab, ylab = ylab[j], type = "n",
         main = p.names[j], ...)
    points(lam[-1] + 1, b[-1, j], pch = 19)
    points(lam[1] + 1, b[1, j], pch = 1)
    lines(a[a >= 0] + 1, d[a >= 0])
    lines(a[a < 0] + 1, d[a < 0], lty = 2)
    points(0, x$coefficients[j], pch = 4, cex = 1.5, col = "red")
  }
}


#' @export
predict.simex <- function(object, newdata = NULL, type = "response", ...) {
  if (is.null(newdata)) {
    if (type == "link") {
      return(object$family$linkfun(object$fitted.values))
    }
    return(object$fitted.values)
  }

  # Use the naive model's formula to build the model matrix for new data
  tt <- terms(object$formula)
  tt <- delete.response(tt)
  mf <- model.frame(tt, data = newdata)
  X_new <- model.matrix(tt, mf)

  eta <- as.numeric(X_new %*% object$coefficients)
  if (type == "link") return(eta)
  object$family$linkinv(eta)
}


#' @export
coef.simex <- function(object, ...) object$coefficients

#' @export
vcov.simex <- function(object, ...) object$vcov

#' @export
residuals.simex <- function(object, ...) object$residuals

#' @export
fitted.simex <- function(object, ...) object$fitted.values

#' @export
nobs.simex <- function(object, ...) object$n

#' @export
family.simex <- function(object, ...) object$family

#' @export
formula.simex <- function(x, ...) x$formula

#' @export
confint.simex <- function(object, parm, level = 0.95, ...) {
  cf <- coef(object)
  if (is.null(object$vcov))
    stop("No variance estimates available. Rerun with jackknife = TRUE.")

  ses <- sqrt(pmax(diag(object$vcov), 0))
  a <- (1 - level) / 2
  fac <- qnorm(c(a, 1 - a))

  pnames <- names(cf)
  if (missing(parm)) {
    parm <- pnames
  } else if (is.numeric(parm)) {
    parm <- pnames[parm]
  }

  ci <- cbind(cf[parm] + fac[1] * ses[parm],
              cf[parm] + fac[2] * ses[parm])
  colnames(ci) <- paste0(format(100 * c(a, 1 - a), trim = TRUE), " %")
  ci
}
