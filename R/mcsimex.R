# =========================================================================
# MC-SIMEX: Misclassification Simulation-Extrapolation
# =========================================================================

#' MC-SIMEX estimator for GLMs with misclassified covariates
#'
#' Implements the Misclassification SIMEX method of Kuchenhoff, Mwalili and
#' Lesaffre (2006), as well as the improved MC-SIMEX of Sevilimedu and Yu
#' (2026). The simulation step runs in C++ for high performance.
#' Supports Gaussian, Poisson, and Binomial GLM families.
#'
#' @param formula a formula describing the model (as in \code{\link{glm}}).
#' @param family a family object (e.g. \code{poisson()}, \code{binomial()}).
#'   Default: \code{gaussian()}.
#' @param data a data frame containing the variables in the formula.
#' @param SIMEXvariable character vector naming the misclassified factor
#'   variable(s) in \code{data}.
#' @param mc.matrix misclassification matrix (K x K) or a named list of
#'   matrices when multiple variables are misclassified. Columns must sum to 1:
#'   \code{mc.matrix[j,l] = P(Z* = j | Z = l)}.
#' @param method character: \code{"standard"} (default) for the original
#'   MC-SIMEX with approximate extrapolation, or \code{"improved"} for the
#'   Sevilimedu & Yu (2026) method using the exact correction factor.
#'   The improved method requires K = 2 (binary misclassified covariate) and
#'   defaults to \code{lambda = 1, B = 1} for maximum efficiency.
#' @param lambda numeric vector of SIMEX exponents. For \code{method =
#'   "standard"}: default \code{c(0.5, 1, 1.5, 2)}. For \code{method =
#'   "improved"}: default \code{1} (single value sufficient). Can also be
#'   \code{"optimal"} to automatically select the lambda that minimizes
#'   the correction factor variance.
#' @param B number of simulation replicates per lambda level. For
#'   \code{method = "standard"}: default 200. For \code{method = "improved"}:
#'   default 1 (single replicate sufficient).
#' @param extrapolation extrapolation method for \code{method = "standard"}:
#'   \code{"quadratic"} (default), \code{"linear"}, or \code{"loglinear"}.
#'   Ignored when \code{method = "improved"}.
#' @param jackknife logical or character. If \code{TRUE}, compute variance
#'   estimates. For \code{method = "improved"}, uses the closed-form
#'   variance from Sevilimedu & Yu (2026). Default \code{TRUE}.
#' @param weights optional prior weights (as in \code{\link{glm}}).
#' @param seed random seed for reproducibility (default: 42).
#'
#' @return An object of class \code{"mcsimex"} with components mimicking
#'   \code{\link{glm}} output:
#'   \describe{
#'     \item{coefficients}{Named vector of SIMEX-corrected coefficients.}
#'     \item{naive.coefficients}{Coefficients from the naive (uncorrected) model.}
#'     \item{residuals}{Working residuals.}
#'     \item{fitted.values}{Fitted values from the corrected model.}
#'     \item{family}{The family object.}
#'     \item{formula}{The model formula.}
#'     \item{call}{The matched call.}
#'     \item{data}{The data used.}
#'     \item{vcov}{Variance-covariance matrix (if requested).}
#'     \item{SIMEX.estimates}{Matrix of averaged estimates at each lambda level.}
#'     \item{theta}{List of B x p matrices of raw simulation estimates.}
#'     \item{extrapolation}{Fitted extrapolation model(s) (standard only).}
#'     \item{lambda}{Lambda vector used (including 0).}
#'     \item{B}{Number of simulation replicates.}
#'     \item{SIMEXvariable}{Name(s) of the misclassified variable(s).}
#'     \item{mc.matrix}{The misclassification matrix/matrices.}
#'     \item{naive.model}{The naive \code{glm} fit.}
#'     \item{extrapolation.method}{Extrapolation method used.}
#'     \item{method}{Either \code{"standard"} or \code{"improved"}.}
#'     \item{c.lambda}{Correction factors at each lambda (improved only).}
#'     \item{pi.x}{Estimated marginal probability P(X=1) (improved only).}
#'   }
#'
#' @details
#' \strong{Standard MC-SIMEX} (Kuchenhoff et al., 2006):
#'
#' The algorithm has two steps: (1) simulate pseudo-datasets with increased
#' misclassification at each lambda level, refit the model B times, and average;
#' (2) extrapolate the averaged coefficients to lambda = -1 using an approximate
#' function (quadratic, linear, or loglinear).
#'
#' \strong{Improved MC-SIMEX} (Sevilimedu & Yu, 2026):
#'
#' Instead of using an approximate extrapolation function, the improved method
#' derives the exact relationship between the naive coefficient and the true
#' coefficient via the misclassification matrix. For a binary covariate X
#' observed as W with misclassification matrix Pi:
#'
#' \deqn{\beta_{na}(\lambda) = \frac{(\pi_{00}(\lambda) - \pi_x \delta_\lambda)
#' (1 - \pi_{00}(\lambda) + \delta_\lambda \pi_x)}{\pi_x(1-\pi_x)\delta_\lambda}
#' \beta_x}
#'
#' where \eqn{\delta_\lambda = \pi_{00}(\lambda) + \pi_{11}(\lambda) - 1} and
#' \eqn{\pi_{00}(\lambda), \pi_{11}(\lambda)} are diagonal elements of
#' \eqn{\Pi^{1+\lambda}}. The correction factor \eqn{c_\lambda} inverts this
#' relationship, and the corrected estimate is simply
#' \eqn{\hat\beta_x = c_\lambda \hat\beta_{na}(\lambda)}.
#'
#' Because the extrapolation is exact, only a single lambda and a single
#' replicate (K=1, B=1) suffice, reducing computation by a factor of ~100.
#'
#' @references
#' Kuchenhoff, H., Mwalili, S. M. and Lesaffre, E. (2006).
#' A general method for dealing with misclassification in regression: The
#' misclassification SIMEX. \emph{Biometrics}, 62(1), 85--96.
#'
#' Sevilimedu, V. and Yu, L. (2026). An improved misclassification simulation
#' extrapolation (MC-SIMEX) algorithm. \emph{Statistics in Medicine}, 45,
#' e70418. \doi{10.1002/sim.70418}
#'
#' @examples
#' # --- Standard MC-SIMEX ---
#' set.seed(42)
#' n <- 2000
#' x <- rnorm(n)
#' z <- rbinom(n, 1, 0.4)
#' y <- rpois(n, exp(0.5 + 0.8 * z + 0.3 * x))
#'
#' p01 <- 0.10; p10 <- 0.15
#' z_star <- z
#' z_star[z == 0] <- rbinom(sum(z == 0), 1, p01)
#' z_star[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)
#'
#' df <- data.frame(y = y, z_star = factor(z_star), x = x)
#' Pi <- matrix(c(1 - p01, p01, p10, 1 - p10), 2, 2)
#'
#' fit_std <- mcsimex(y ~ z_star + x, family = poisson(), data = df,
#'                    SIMEXvariable = "z_star", mc.matrix = Pi, B = 100)
#'
#' # --- Improved MC-SIMEX (much faster) ---
#' fit_imp <- mcsimex(y ~ z_star + x, family = poisson(), data = df,
#'                    SIMEXvariable = "z_star", mc.matrix = Pi,
#'                    method = "improved")
#' summary(fit_imp)
#'
#' @export
mcsimex <- function(formula, family = gaussian(), data,
                    SIMEXvariable, mc.matrix,
                    method = c("standard", "improved"),
                    lambda = NULL,
                    B = NULL,
                    extrapolation = c("quadratic", "linear", "loglinear"),
                    jackknife = TRUE,
                    weights = NULL,
                    seed = 42L) {

  cl <- match.call()
  method <- match.arg(method)
  extrapolation <- match.arg(extrapolation)

  # --- defaults depend on method ---
  if (is.null(lambda)) {
    lambda <- if (method == "improved") 1 else c(0.5, 1, 1.5, 2)
  }
  if (is.null(B)) {
    B <- if (method == "improved") 1L else 200L
  }

  # --- resolve family ---
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())()
  if (is.function(family)) family <- family()
  if (is.null(family$family))
    stop("'family' not recognized")

  dist_code <- family_to_dist_code(family)

  # --- validate SIMEXvariable ---
  if (!is.character(SIMEXvariable))
    stop("SIMEXvariable must be a character string")
  if (length(SIMEXvariable) != 1L)
    stop("Currently only one SIMEXvariable is supported for MC-SIMEX")

  # --- ensure factor ---
  if (!is.factor(data[[SIMEXvariable]]))
    data[[SIMEXvariable]] <- factor(data[[SIMEXvariable]])

  # --- validate mc.matrix ---
  if (is.matrix(mc.matrix)) {
    mc_list <- list(mc.matrix)
    names(mc_list) <- SIMEXvariable
  } else if (is.list(mc.matrix)) {
    mc_list <- mc.matrix
  } else {
    stop("mc.matrix must be a matrix or named list of matrices")
  }

  Pi <- as.matrix(mc_list[[SIMEXvariable]])
  K <- nrow(Pi)
  stopifnot(ncol(Pi) == K)

  col_sums <- colSums(Pi)
  if (any(abs(col_sums - 1) > 1e-6))
    stop("Columns of mc.matrix must sum to 1 (got: ",
         paste(round(col_sums, 4), collapse = ", "), ")")

  # --- improved method requires K = 2 ---
  if (method == "improved" && K != 2L)
    stop("method = 'improved' requires a binary misclassified covariate (K = 2)")

  # --- handle lambda = "optimal" for improved method ---
  optimal_lambda <- FALSE
  if (method == "improved" && is.character(lambda) && lambda == "optimal") {
    optimal_lambda <- TRUE
    lambda <- 1  # placeholder, will be replaced below
  }

  # --- fit naive model ---
  if (!is.null(weights)) {
    naive_fit <- glm(formula, family = family, data = data, weights = weights,
                     x = TRUE, y = TRUE)
  } else {
    naive_fit <- glm(formula, family = family, data = data,
                     x = TRUE, y = TRUE)
  }

  # --- extract data components for C++ ---
  y <- naive_fit$y
  n <- length(y)
  z_factor <- data[[SIMEXvariable]]
  z_levels <- levels(z_factor)

  if (length(z_levels) != K)
    stop("Number of levels in SIMEXvariable (", length(z_levels),
         ") does not match mc.matrix dimensions (", K, ")")

  z_hat <- as.integer(z_factor) - 1L

  other_terms <- setdiff(all.vars(formula)[-1], SIMEXvariable)
  if (length(other_terms) > 0) {
    other_formula <- reformulate(other_terms, intercept = TRUE)
    x_mat <- model.matrix(other_formula, data = data)
  } else {
    x_mat <- matrix(1, nrow = n, ncol = 1)
    colnames(x_mat) <- "(Intercept)"
  }

  wt <- if (is.null(weights)) rep(1.0, n) else as.numeric(weights)

  B <- as.integer(B)

  # --- naive coefficients ---
  xi_hat <- .build_xi_hat(z_hat, x_mat, K)
  p <- ncol(xi_hat)

  fam <- get_link_funs(family)
  dat_naive <- data.frame(y = y, xi_hat)
  naive_refit <- glm(y ~ . - 1, data = dat_naive, family = family, weights = wt)
  psi_naive <- unname(coef(naive_refit))

  nms <- .make_param_names(z_levels, x_mat, K)
  names(psi_naive) <- nms

  # =====================================================================
  # IMPROVED MC-SIMEX (Sevilimedu & Yu, 2026)
  # =====================================================================
  if (method == "improved") {
    # --- estimate pi_x = P(X = 1) from observed data ---
    pi_x <- .estimate_pi_x(z_hat, Pi)

    # --- find optimal lambda if requested ---
    if (optimal_lambda) {
      lambda <- .find_optimal_lambda(Pi, pi_x)
    }

    lambda <- sort(as.numeric(lambda))
    n_lambda <- length(lambda)

    # --- compute correction factors c_lambda ---
    c_lam_vec <- .compute_c_lambda(Pi, pi_x, lambda)

    # --- intercept correction factors ---
    intercept_correction <- .compute_intercept_correction(Pi, pi_x, lambda)

    # --- simulation step (same C++ code) ---
    sim_result <- mcsimex_sim_cpp(y, z_hat, x_mat, Pi, K, dist_code,
                                  lambda, B, wt, as.integer(seed))

    # --- apply exact correction ---
    # Column 1 = SIMEXvariable coefficient (dummy for z=1)
    # Column 2 = intercept, columns 3+ = confounders
    s <- K - 1L  # number of dummy columns = 1 for K=2

    # Reshape and correct
    theta_list <- vector("list", n_lambda)
    corrected_all <- matrix(0, nrow = 0, ncol = p)

    for (l in seq_len(n_lambda)) {
      rows <- ((l - 1) * B + 1):(l * B)
      theta_mat <- sim_result[rows, , drop = FALSE]
      colnames(theta_mat) <- nms

      # Correct SIMEXvariable coefficient(s): multiply by c_lambda
      theta_corrected <- theta_mat
      for (k in seq_len(s)) {
        theta_corrected[, k] <- c_lam_vec[l] * theta_mat[, k]
      }

      # Correct intercept: beta_0 = beta_0* - beta_x * correction_factor
      icol <- s + 1  # intercept column index
      for (b_idx in seq_len(nrow(theta_corrected))) {
        beta_x_corrected <- theta_corrected[b_idx, 1]
        theta_corrected[b_idx, icol] <- theta_mat[b_idx, icol] -
          beta_x_corrected * intercept_correction[l]
      }

      # Confounders (columns s+2 to p): leave as-is from simulation

      theta_list[[l]] <- theta_corrected
      corrected_all <- rbind(corrected_all, theta_corrected)
    }
    names(theta_list) <- paste0("lambda_", lambda)

    # Final corrected coefficients: average over all (lambda, B) combinations
    corrected_coefs <- colMeans(corrected_all)
    names(corrected_coefs) <- nms

    # Build SIMEX.estimates matrix (for compatibility)
    avg_estimates <- matrix(0, n_lambda, p)
    for (l in seq_len(n_lambda)) {
      avg_estimates[l, ] <- colMeans(theta_list[[l]])
    }
    all_lambda <- c(0, lambda)
    all_estimates <- rbind(psi_naive, avg_estimates)
    rownames(all_estimates) <- paste0("lambda_", all_lambda)
    colnames(all_estimates) <- nms

    # --- variance estimation (closed-form) ---
    vcov_improved <- NULL
    if (!isFALSE(jackknife)) {
      vcov_improved <- .variance_improved(theta_list, corrected_coefs,
                                          c_lam_vec, n_lambda, B, p,
                                          naive_refit, xi_hat, y, wt, fam)
      rownames(vcov_improved) <- colnames(vcov_improved) <- nms
    }

    # --- fitted values and residuals ---
    eta_corrected <- as.numeric(xi_hat %*% corrected_coefs)
    fitted_vals <- family$linkinv(eta_corrected)
    resids <- y - fitted_vals

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
        vcov               = vcov_improved,
        SIMEX.estimates    = all_estimates,
        theta              = theta_list,
        extrapolation      = NULL,
        lambda             = all_lambda,
        B                  = B,
        n                  = n,
        p                  = p,
        SIMEXvariable      = SIMEXvariable,
        mc.matrix          = mc_list,
        naive.model        = naive_fit,
        extrapolation.method = "exact (improved)",
        method             = "improved",
        c.lambda           = setNames(c_lam_vec, paste0("lambda_", lambda)),
        pi.x               = pi_x
      ),
      class = "mcsimex"
    )
    return(out)
  }

  # =====================================================================
  # STANDARD MC-SIMEX (Kuchenhoff et al., 2006)
  # =====================================================================
  lambda <- sort(as.numeric(lambda))
  n_lambda <- length(lambda)

  # --- simulation step (C++) ---
  sim_result <- mcsimex_sim_cpp(y, z_hat, x_mat, Pi, K, dist_code,
                                lambda, B, wt, as.integer(seed))

  # Reshape: list of B x p matrices (one per lambda)
  theta_list <- vector("list", n_lambda)
  avg_estimates <- matrix(0, n_lambda, p)
  for (l in seq_len(n_lambda)) {
    rows <- ((l - 1) * B + 1):(l * B)
    theta_mat <- sim_result[rows, , drop = FALSE]
    colnames(theta_mat) <- nms
    theta_list[[l]] <- theta_mat
    avg_estimates[l, ] <- colMeans(theta_mat)
  }
  names(theta_list) <- paste0("lambda_", lambda)

  # Stack: lambda=0 (naive) + simulation averages
  all_lambda <- c(0, lambda)
  all_estimates <- rbind(psi_naive, avg_estimates)
  rownames(all_estimates) <- paste0("lambda_", all_lambda)
  colnames(all_estimates) <- nms

  # --- extrapolation step ---
  extrap_result <- extrapolate_simex(all_lambda, all_estimates, extrapolation)

  # --- jackknife variance ---
  vcov_jk <- NULL
  if (!isFALSE(jackknife)) {
    jk_method <- if (is.character(jackknife)) jackknife else extrapolation
    vcov_jk <- jackknife_variance_mcsimex(theta_list, psi_naive, naive_refit,
                                          lambda, jk_method, fam, xi_hat, y, wt)
    rownames(vcov_jk) <- colnames(vcov_jk) <- nms
  }

  # --- fitted values and residuals ---
  corrected_coefs <- extrap_result$coefficients
  eta_corrected <- as.numeric(xi_hat %*% corrected_coefs)
  fitted_vals <- family$linkinv(eta_corrected)
  resids <- y - fitted_vals

  # --- build output object ---
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
      mc.matrix          = mc_list,
      naive.model        = naive_fit,
      extrapolation.method = extrapolation,
      method             = "standard"
    ),
    class = "mcsimex"
  )
  out
}


# =========================================================================
# Internal helpers
# =========================================================================

#' Build design matrix [dummies(z, K-1), x] in R (mirrors C++ build_design)
#' @keywords internal
.build_xi_hat <- function(z_hat, x_mat, K) {
  n <- length(z_hat)
  s <- K - 1L
  dummies <- matrix(0, n, s)
  for (k in seq_len(s)) {
    dummies[, k] <- as.numeric(z_hat == k)
  }
  cbind(dummies, x_mat)
}

#' Generate parameter names for the SIMEX model
#' @keywords internal
.make_param_names <- function(z_levels, x_mat, K) {
  if (K == 2L) {
    dummy_names <- z_levels[2]
  } else {
    dummy_names <- z_levels[-1]
  }
  c(dummy_names, colnames(x_mat))
}


# =========================================================================
# Improved MC-SIMEX helpers (Sevilimedu & Yu, 2026)
# =========================================================================

#' Estimate marginal probability P(X=1) from observed data and Pi
#'
#' Uses the relationship: P(W=1) = (1-pi_00)(1-pi_x) + pi_11 pi_x
#' Solving: pi_x = (P(W=1) - 1 + pi_00) / (pi_00 + pi_11 - 1)
#' @keywords internal
.estimate_pi_x <- function(z_hat, Pi) {
  p_w1 <- mean(z_hat == 1L)
  pi_00 <- Pi[1, 1]
  pi_11 <- Pi[2, 2]
  delta <- pi_00 + pi_11 - 1

  if (abs(delta) < 1e-10)
    stop("Misclassification matrix is degenerate (pi_00 + pi_11 = 1)")

  pi_x <- (p_w1 - 1 + pi_00) / delta
  # Clamp to (0, 1)
  pi_x <- max(1e-6, min(1 - 1e-6, pi_x))
  pi_x
}


#' Compute matrix power via eigendecomposition (R version)
#' @keywords internal
.mat_power_r <- function(Pi, power) {
  if (power == 0) return(diag(nrow(Pi)))
  if (power == 1) return(Pi)

  ev <- eigen(Pi)
  V <- ev$vectors
  d <- ev$values

  d_pow <- abs(d)^power * sign(d)
  zapsmall(V %*% diag(d_pow) %*% solve(V))
}


#' Compute correction factor c_lambda for the improved MC-SIMEX
#'
#' From Sevilimedu & Yu (2026), equation (3):
#' beta_na(lambda) = attenuation * beta_x, so
#' beta_x = c_lambda * beta_na(lambda) where c_lambda = 1/attenuation.
#'
#' The correction inverts the attenuation:
#' c_lambda = ((pi_00(lambda) - pi_x delta)(1 - pi_00(lambda) + delta pi_x)) /
#'            (pi_x (1 - pi_x) delta)
#'
#' where pi_00(lambda), pi_11(lambda) are from Pi^(1+lambda) and
#' delta = pi_00(lambda) + pi_11(lambda) - 1.
#'
#' @param Pi 2x2 misclassification matrix
#' @param pi_x marginal probability P(X=1)
#' @param lambda vector of lambda values
#' @return numeric vector of correction factors (> 1 for attenuation correction)
#' @keywords internal
.compute_c_lambda <- function(Pi, pi_x, lambda) {
  sapply(lambda, function(lam) {
    Pi_lam <- .mat_power_r(Pi, 1 + lam)
    pi_00_lam <- Pi_lam[1, 1]
    pi_11_lam <- Pi_lam[2, 2]
    delta_lam <- pi_00_lam + pi_11_lam - 1

    if (abs(delta_lam) < 1e-10) return(NA_real_)

    denominator <- pi_x * (1 - pi_x) * delta_lam
    numerator <- (pi_00_lam - pi_x * delta_lam) *
                 (1 - pi_00_lam + delta_lam * pi_x)

    if (abs(denominator) < 1e-10) return(NA_real_)

    numerator / denominator
  })
}


#' Compute intercept correction factor for the improved MC-SIMEX
#'
#' From the derivation: beta_0 = beta_0* - beta_x * (1-pi_11(lambda)) pi_x / (pi_00(lambda) - delta_lambda pi_x)
#' The correction factor is: (1-pi_11(lambda)) pi_x / (pi_00(lambda) - delta_lambda pi_x)
#'
#' @keywords internal
.compute_intercept_correction <- function(Pi, pi_x, lambda) {
  sapply(lambda, function(lam) {
    Pi_lam <- .mat_power_r(Pi, 1 + lam)
    pi_00_lam <- Pi_lam[1, 1]
    pi_11_lam <- Pi_lam[2, 2]
    delta_lam <- pi_00_lam + pi_11_lam - 1

    denom <- pi_00_lam - delta_lam * pi_x
    if (abs(denom) < 1e-10) return(0)

    (1 - pi_11_lam) * pi_x / denom
  })
}


#' Find optimal lambda that minimizes |c_lambda|
#'
#' Searches over a grid to find the lambda value producing the smallest
#' absolute correction factor, which corresponds to minimum variance.
#'
#' @keywords internal
.find_optimal_lambda <- function(Pi, pi_x, grid = seq(0.1, 3, by = 0.01)) {
  c_values <- .compute_c_lambda(Pi, pi_x, grid)
  valid <- !is.na(c_values) & is.finite(c_values)
  if (!any(valid)) stop("Cannot find valid lambda for the improved method")

  grid[valid][which.min(abs(c_values[valid]))]
}


#' Variance estimation for the improved MC-SIMEX
#'
#' From Sevilimedu & Yu (2026), equation (4):
#' Var(beta_hat^sim) = 1/(B^2 K^2) sum_k sum_b c_lambda_k^2 Var(beta_hat_na(W_bk))
#'
#' We estimate Var(beta_hat_na) using the sandwich variance from the naive model
#' at each lambda level, and the empirical covariance of corrected estimates.
#'
#' @keywords internal
.variance_improved <- function(theta_list, corrected_coefs, c_lam_vec,
                               n_lambda, B, p, naive_fit, xi_hat, y, wt, fam) {
  N <- sum(wt)
  n <- length(y)

  if (n_lambda * B == 1) {
    # Single replicate: use scaled naive model variance
    # Var(beta_x^sim) = c_lambda^2 * Var(beta_na)
    V_naive <- vcov(naive_fit)
    # Map naive vcov to the improved estimate layout
    # The naive model's x matrix matches xi_hat
    psi_mean <- corrected_coefs
    eta_m <- as.numeric(xi_hat %*% psi_mean)
    w_m <- fam$mu_dot(eta_m)
    eps_m <- y - fam$mu(eta_m)
    A_m <- crossprod(xi_hat * (wt * w_m), xi_hat) / N
    C_m <- crossprod(xi_hat * (wt * eps_m), xi_hat * eps_m) / N
    A_m_inv <- tryCatch(solve(A_m), error = function(e) MASS::ginv(A_m))
    V_model <- A_m_inv %*% C_m %*% A_m_inv / N

    # Scale by c_lambda^2 for the SIMEXvariable coefficient
    V_scaled <- V_model * c_lam_vec[1]^2
    return(V_scaled)
  }

  # Multiple replicates: use empirical covariance of corrected estimates
  # Stack all corrected estimates
  all_corrected <- do.call(rbind, theta_list)

  # Empirical variance of the mean
  # Var(mean) = Var(estimates) / (B * K)
  V_emp <- cov(all_corrected) / (n_lambda * B)

  V_emp
}


# =========================================================================
# S3 methods for class "mcsimex"
# =========================================================================

#' @export
print.mcsimex <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall:\n")
  print(x$call)
  method_label <- if (!is.null(x$method) && x$method == "improved")
    "Improved MC-SIMEX" else "MC-SIMEX"
  cat(sprintf("\n%s corrected coefficients:\n", method_label))
  print.default(format(coef(x), digits = digits), print.gap = 2,
                quote = FALSE)
  cat("\n")
  invisible(x)
}


#' @export
summary.mcsimex <- function(object, ...) {
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
  ans$mc.matrix <- object$mc.matrix
  ans$lambda <- object$lambda
  ans$B <- object$B
  ans$extrapolation.method <- object$extrapolation.method
  ans$method <- object$method
  ans$SIMEX.estimates <- object$SIMEX.estimates
  ans$naive.coefficients <- object$naive.coefficients
  ans$residuals <- object$residuals
  ans$n <- n
  ans$df.residual <- rdf
  ans$c.lambda <- object$c.lambda
  ans$pi.x <- object$pi.x

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

  class(ans) <- "summary.mcsimex"
  ans
}


#' @export
print.summary.mcsimex <- function(x,
                                  digits = max(3, getOption("digits") - 3),
                                  ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\nFamily:", x$family$family, "\n")
  cat("SIMEX variable:", x$SIMEXvariable, "\n")

  method_label <- if (!is.null(x$method) && x$method == "improved")
    "improved (Sevilimedu & Yu, 2026)" else "standard"
  cat("Method:", method_label, "\n")
  cat("Extrapolation:", x$extrapolation.method, "\n")
  cat("Lambda grid:", paste(x$lambda, collapse = ", "), "\n")
  cat("B =", x$B, ", n =", x$n, "\n")

  if (!is.null(x$pi.x)) {
    cat("Estimated P(X=1):", round(x$pi.x, 4), "\n")
  }
  if (!is.null(x$c.lambda)) {
    cat("Correction factor(s) c_lambda:", paste(round(x$c.lambda, 4), collapse = ", "), "\n")
  }

  if (!is.null(x$residuals)) {
    cat("\nResiduals:\n")
    rq <- quantile(x$residuals, na.rm = TRUE)
    names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
    print(rq, digits = digits)
  }

  cat("\nNaive coefficients:\n")
  print(round(x$naive.coefficients, digits))

  cat("\nMC-SIMEX corrected coefficients:\n")
  if (ncol(x$coefficients) > 1) {
    printCoefmat(x$coefficients, digits = digits)
  } else {
    print(round(x$coefficients, digits))
  }
  cat("\n")
  invisible(x)
}


#' @export
plot.mcsimex <- function(x,
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

    if (!is.null(x$extrapolation) && !is.null(x$extrapolation[[j]])) {
      # Standard method: plot extrapolation curve
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
    } else {
      # Improved method: plot simulation points + corrected value
      plot(lam + 1, b[, j], xlab = xlab, ylab = ylab[j], type = "n",
           main = paste(p.names[j], "(improved)"), ...)
      points(lam[-1] + 1, b[-1, j], pch = 19)
      points(lam[1] + 1, b[1, j], pch = 1)
      abline(h = x$coefficients[j], lty = 2, col = "red")
      points(0, x$coefficients[j], pch = 4, cex = 1.5, col = "red")
    }
  }
}


#' @export
predict.mcsimex <- function(object, newdata = NULL, type = "response", ...) {
  if (is.null(newdata)) {
    if (type == "link") {
      return(object$family$linkfun(object$fitted.values))
    }
    return(object$fitted.values)
  }

  if (!is.factor(newdata[[object$SIMEXvariable]]))
    newdata[[object$SIMEXvariable]] <- factor(newdata[[object$SIMEXvariable]])

  z_new <- as.integer(newdata[[object$SIMEXvariable]]) - 1L
  K <- object$mc.matrix[[object$SIMEXvariable]] |> nrow()
  other_terms <- setdiff(all.vars(object$formula)[-1], object$SIMEXvariable)
  if (length(other_terms) > 0) {
    other_formula <- reformulate(other_terms, intercept = TRUE)
    x_new <- model.matrix(other_formula, data = newdata)
  } else {
    x_new <- matrix(1, nrow = nrow(newdata), ncol = 1)
  }
  xi_new <- .build_xi_hat(z_new, x_new, K)

  eta <- as.numeric(xi_new %*% object$coefficients)
  if (type == "link") return(eta)
  object$family$linkinv(eta)
}


#' @export
coef.mcsimex <- function(object, ...) object$coefficients

#' @export
vcov.mcsimex <- function(object, ...) object$vcov

#' @export
residuals.mcsimex <- function(object, ...) object$residuals

#' @export
fitted.mcsimex <- function(object, ...) object$fitted.values

#' @export
nobs.mcsimex <- function(object, ...) object$n

#' @export
family.mcsimex <- function(object, ...) object$family

#' @export
formula.mcsimex <- function(x, ...) x$formula

#' @export
confint.mcsimex <- function(object, parm, level = 0.95, ...) {
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
