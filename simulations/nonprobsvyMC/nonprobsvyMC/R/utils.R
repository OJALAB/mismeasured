# ---------------------------------------------------------------------------
# Internal helpers for nonprob_mc()
# ---------------------------------------------------------------------------

#' Name of the variable inside the mc() term of a formula
#' @keywords internal
.mc_variable_name <- function(formula) {
  found <- NULL
  walk <- function(node) {
    if (is.call(node)) {
      op <- deparse(node[[1]])
      if (op == "mc") {
        found <<- deparse(node[[2]])
        return(invisible(NULL))
      }
      for (i in seq_along(node)[-1]) walk(node[[i]])
    }
    invisible(NULL)
  }
  walk(formula[[3]])
  if (is.null(found))
    stop("'outcome' must contain exactly one mc() term, ",
         "e.g. y ~ mc(z, Pi) + x1.")
  found
}

#' Design matrix in mcglm coefficient order: (z dummies, x)
#'
#' For K = 2 the proxy enters as a single column; for K > 2 as K - 1
#' dummies with baseline 0, matching the (gamma, alpha) coefficient
#' layout of \code{mcglm}.
#' @keywords internal
.build_xi <- function(z, x, K) {
  if (K == 2L) return(cbind(z, x))
  n <- length(z)
  d <- matrix(0, n, K - 1L)
  for (k in seq_len(K - 1L)) d[, k] <- as.numeric(z == k)
  cbind(d, x)
}

#' Outcome-model component of the mass-imputation variance
#'
#' Corrected-score linearization
#' \eqn{G_P^\top B \hat S B^\top G_P / n_{NP}} with \eqn{B} and
#' \eqn{\hat S} the sandwich bread and meat of the \code{mcglm} fit and
#' \eqn{G_P} the design-weighted mean of \eqn{\dot\mu\,\xi} over the
#' probability sample.
#' @keywords internal
.var_nonprob_mc <- function(model_fitted, mc_method, parts_rand,
                            svydesign, family_outcome) {
  B <- bread(model_fitted, method = mc_method)
  S <- meat(model_fitted, method = mc_method)

  xi_rand <- .build_xi(parts_rand$z_hat, parts_rand$x, model_fitted$K)
  psi     <- coef(model_fitted, method = mc_method)
  eta     <- as.numeric(xi_rand %*% psi)

  fam   <- switch(family_outcome,
                  gaussian = gaussian(),
                  binomial = binomial(),
                  poisson  = poisson())
  w     <- weights(svydesign)
  G_P   <- colSums(w * fam$mu.eta(eta) * xi_rand) / sum(w)

  as.numeric(t(G_P) %*% B %*% S %*% t(B) %*% G_P) / stats::nobs(model_fitted)
}
