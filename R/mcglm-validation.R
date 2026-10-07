# ---------------------------------------------------------------------------
# Main study / validation study designs for the regression estimators
# (Yi, Yan, Liao and Spiegelman, 2019, Section 4).
#
# The misclassification nuisance eta comes from estimate_mc() (an
# "mc_estimate" object with design-weighted estimating functions). The
# regression parameter psi solves, over the main-study units,
#
#   external validation:          sum_M w_i U_i(psi, eta) = 0,
#   internal, beta_equation "yi":  sum_{M \ V} w_i U_i + sum_V w_i S_i = 0,
#   internal, "weighted":          sum_M w_i U_i + sum_V w_i d_i (S_i - U_i) = 0,
#
# where U is the method's corrected estimating function, S the score at the
# true category, w main-study frequency weights and d design weights.
# "yi" (eq. 23) needs selection into V to depend on (Z_hat, x) only;
# "weighted" is unbiased for any known inclusion probabilities.
#
# Inference stacks theta = (psi, eta):
#   Var(theta_hat) = A^{-1} B A^{-T},  A = sum_u w_u d g_u / d theta',
# with B from .mc_meat() (validation rows merged into their main unit for
# internal designs; stratum-centred for stratified external designs).
# ---------------------------------------------------------------------------

# Methods whose misclassification nuisance is estimated from a validation
# sample (rather than supplied) when mcglm(validation = ) is used.
.mcglm_validated_methods <- c("sub", "ec", "il", "cs_akn")

#' Parse and check a validation-sample description
#'
#' @param validation \code{list(z, z_hat)} (external),
#'   \code{list(z, index)} (internal) or a \code{\link{validation_sample}}.
#' @param z_hat Main-study proxy codes (0-based).
#' @param K Number of categories.
#' @return The bound sample from \code{.mc_bind_validation}.
#' @keywords internal
.mcglm_parse_validation <- function(validation, z_hat, K) {
  .mc_bind_validation(as_validation_sample(validation), z_hat, K)
}

#' Score at the true category for internal validation rows
#'
#' @return List with \code{U} (rows \eqn{\xi_i\{Y_i - \mu(\psi^\top\xi_i)\}})
#'   and \code{xi}, \code{mud} for the Jacobian.
#' @keywords internal
.mcglm_true_parts <- function(psi, y, z, x, K, fam) {
  xi  <- .mcglm_build_xi_hat(z, x, K)
  eta <- as.numeric(xi %*% psi)
  list(U = xi * (y - fam$linkinv(eta)), xi = xi, mud = fam$mu.eta(eta))
}

#' Sum-scale Jacobian of the true-score rows
#' @keywords internal
.mcglm_true_jacobian <- function(parts, w) {
  -crossprod(parts$xi * (w * parts$mud), parts$xi)
}

#' Row-specific predictive probabilities for the main study
#'
#' Adds \code{Wi} (n x K, \eqn{\Pr(Z = \ell \mid \hat Z_i, x_i)}) to the
#' parameters returned by an estimate's \code{map()}: the proxy's row of the
#' predictive matrix for a constant prevalence, or
#' \eqn{\Pi_{\hat z_i \ell}\pi_\ell(x_i) / \sum_m \Pi_{\hat z_i m}\pi_m(x_i)}
#' under a prevalence model (with \code{P}, the n x K prevalences).
#' @keywords internal
.mc_with_rows <- function(par, z_hat, w_main = NULL) {
  if (is.null(par$alpha)) {
    par$Wi <- par$W[z_hat + 1L, , drop = FALSE]
    return(par)
  }
  par$P <- .mc_prev_probs(par$alpha, w_main)
  joint <- par$Pi[z_hat + 1L, , drop = FALSE] * par$P
  par$Wi <- joint / rowSums(joint)
  par
}

#' Check that an estimate object belongs to the main study being fitted
#' @keywords internal
.mcglm_check_estimate <- function(est, z_hat, K, wt) {
  if (est$K != K)
    stop("The estimate_mc() object has K = ", est$K, " but the model has K = ",
         K, ".", call. = FALSE)
  if (!is.null(est$prevalence$model) &&
      (is.null(est$w_main) || nrow(est$w_main) != length(z_hat)))
    stop("The estimate_mc() object's prevalence model has no covariates for ",
         "the data being fitted.", call. = FALSE)
  if (isTRUE(est$main_dependent)) {
    if (!identical(as.integer(est$z_hat), as.integer(z_hat)))
      stop("The estimate_mc() object was computed with different main-study ",
           "proxies (internal validation or prevalence = 'em'/'inverse' ",
           "use them); pass the fitted data's z_hat to estimate_mc().",
           call. = FALSE)
    same_w <- if (is.null(est$main_weights)) is.null(wt) || all(wt == 1) else
      !is.null(wt) && isTRUE(all.equal(est$main_weights, wt))
    if (!same_w)
      stop("The estimate_mc() object was computed with different main-study ",
           "weights.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Fit a corrected estimating equation with estimated misclassification
#'
#' Generic solver and stacked sandwich for a method whose corrected
#' estimating function depends on the misclassification probabilities.
#'
#' @param psi_init Starting values (typically the naive estimate).
#' @param est An \code{"mc_estimate"} object.
#' @param U_fun \code{function(psi, par)} returning the n x p matrix of
#'   corrected estimating-function rows for every main-study row, where
#'   \code{par = est$map(eta)} has \code{Pi}, \code{pi}, \code{W}.
#' @param J_fun \code{function(psi, par, w)}: sum-scale Jacobian of
#'   \code{colSums(w * U_fun(psi, par))} with respect to \code{psi}.
#' @param control A \code{\link{control_mc}} object (\code{beta_equation},
#'   \code{variance}).
#' @param label Method label for messages.
#' @param true_fun Optional list of \code{U(psi)} (rows for the internal
#'   validation units at their true category) and \code{J(psi, w)} (their
#'   sum-scale Jacobian). Defaults to the canonical GLM score
#'   \eqn{\xi_i\{Y_i - \mu(\psi^\top\xi_i)\}}.
#' @return List with \code{coefficients}, convergence fields, \code{vcov}
#'   (psi block), \code{nuisance} and \code{beta_equation}.
#' @keywords internal
.mcglm_fit_validated <- function(psi_init, est, y, x, z_hat, K, fam,
                                 U_fun, J_fun, wt = NULL,
                                 control = control_mc(), label = "method",
                                 true_fun = NULL) {
  .mcglm_check_estimate(est, z_hat, K, wt)
  vb  <- est$validation
  n   <- length(y)
  p   <- length(psi_init)
  q   <- length(est$eta)
  wt_m <- if (is.null(wt)) rep(1, n) else wt
  internal <- vb$type == "internal"

  beq <- control$beta_equation
  if (beq == "auto") beq <- if (vb$user_weights) "weighted" else "yi"
  if (internal) {
    idx <- vb$index
    y_v <- y[idx]
    x_v <- x[idx, , drop = FALSE]
    d_main <- numeric(n)
    d_main[idx] <- vb$d
    w_U <- wt_m
    if (beq == "yi") w_U[idx] <- 0
    if (is.null(true_fun))
      true_fun <- list(
        U = function(psi) .mcglm_true_parts(psi, y_v, vb$z, x_v, K, fam)$U,
        J = function(psi, w)
          .mcglm_true_jacobian(.mcglm_true_parts(psi, y_v, vb$z, x_v, K,
                                                 fam), w))
  }

  psi_rows <- function(psi, par) {
    G <- U_fun(psi, par)
    if (internal) {
      S <- true_fun$U(psi)
      G[idx, ] <- if (beq == "yi") S else
        G[idx, , drop = FALSE] + vb$d * (S - G[idx, , drop = FALSE])
    }
    G
  }
  psi_jac <- function(psi, par) {
    if (!internal) return(J_fun(psi, par, wt_m))
    if (beq == "yi")
      J_fun(psi, par, w_U) + true_fun$J(psi, wt_m[idx])
    else
      J_fun(psi, par, wt_m) - J_fun(psi, par, d_main * wt_m) +
        true_fun$J(psi, vb$d * wt_m[idx])
  }

  par0 <- .mc_with_rows(est$map(est$eta), z_hat, est$w_main)
  N <- sum(wt_m)
  sol <- nleqslv::nleqslv(
    psi_init,
    function(psi) colSums(wt_m * psi_rows(psi, par0)) / N,
    jac = function(psi) psi_jac(psi, par0) / N,
    control = list(maxit = 500, ftol = 1e-12))
  if (sol$termcd > 2)
    warning(label, " solver did not converge (termcd = ", sol$termcd, ")")
  psi <- sol$x

  # Stacked sandwich over (psi, eta).
  A_pp <- psi_jac(psi, par0)
  A_pe <- .mc_num_jacobian(
    function(e) colSums(wt_m * psi_rows(psi, .mc_with_rows(est$map(e), z_hat,
                                                            est$w_main))),
    unname(est$eta))
  A <- rbind(cbind(A_pp, A_pe), cbind(matrix(0, q, p), est$A))
  r <- est$rows(est$eta)
  G_m <- cbind(psi_rows(psi, par0),
               if (is.null(r$m)) matrix(0, n, q) else r$m)
  G_v <- cbind(matrix(0, vb$n, p), r$v)
  B <- .mc_meat(G_m, G_v, vb, wt_m, n)

  V_psi <- .mc_sandwich_or_na(
    if (control$variance == "conditional") A_pp else A,
    if (control$variance == "conditional")
      B[seq_len(p), seq_len(p), drop = FALSE] else B,
    label)[seq_len(p), seq_len(p), drop = FALSE]

  list(coefficients = psi, converged = sol$termcd <= 2,
       termcd = sol$termcd, iterations = sol$iter, vcov = V_psi,
       beta_equation = if (internal) beq else NA_character_,
       nuisance = list(eta = est$eta, vcov = est$vcov, Pi = est$Pi,
                       pi_z = est$pi, W = est$W,
                       prevalence = est$prevalence$method,
                       variance = control$variance,
                       beta_equation = if (internal) beq else NA_character_))
}
