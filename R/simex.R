# =========================================================================
# Unified SIMEX / MC-SIMEX estimator
# =========================================================================

#' SIMEX estimator for GLMs with measurement error or misclassification
#'
#' Fits a generalized linear model corrected for measurement error or
#' misclassification using the Simulation-Extrapolation (SIMEX) method.
#' Error-prone variables are specified directly in the formula using
#' \code{\link{me}} (continuous measurement error) or \code{\link{mc}}
#' (discrete misclassification).
#'
#' @param formula a formula with \code{\link{me}} and/or \code{\link{mc}} terms.
#'   See Examples.
#' @param family a GLM family (default: \code{gaussian()}).
#' @param data a data frame.
#' @param method character: correction method for \code{mc()} terms.
#'   \code{NULL} (default) auto-selects \code{"improved"}.
#'   Ignored for \code{me()} terms. See Details.
#' @param lambda numeric vector of SIMEX exponents, or \code{"optimal"}
#'   (improved MC-SIMEX only). Default: auto-set based on error type.
#' @param B integer number of simulation replicates. Default: auto-set.
#' @param extrapolation extrapolation method: \code{"quadratic"} (default),
#'   \code{"linear"}, or \code{"loglinear"}. Ignored for improved MC-SIMEX.
#' @param jackknife logical: compute variance estimates? Default \code{TRUE}.
#' @param weights optional prior weights.
#' @param seed random seed (default: 42).
#'
#' @return An object of class \code{"simex"}.
#'
#' @details
#' The function auto-detects the error type from the formula:
#'
#' \itemize{
#'   \item \strong{Continuous measurement error} (\code{me()} terms): Uses the
#'     SIMEX algorithm of Cook and Stefanski (1994). Extra noise is added to
#'     the error-prone variable, the model is refitted B times per lambda level,
#'     and coefficients are extrapolated to lambda = -1.
#'   \item \strong{Discrete misclassification} (\code{mc()} terms): Uses the
#'     MC-SIMEX algorithm. With \code{method = "improved"} (default), the exact
#'     fixed-matrix correction of Sevilimedu and Yu (2026) and its K-level
#'     dummy-vector extension are applied for one misclassified covariate,
#'     requiring only B = 1 replicate. With \code{method = "standard"}, the original
#'     Kuchenhoff et al. (2006) extrapolation-based approach is used.
#' }
#'
#' @references
#' Cook, J.R. and Stefanski, L.A. (1994). Simulation-extrapolation estimation
#' in parametric measurement error models. \emph{JASA}, 89, 1314--1328.
#'
#' Kuchenhoff, H., Mwalili, S.M. and Lesaffre, E. (2006). A general method for
#' dealing with misclassification in regression: The misclassification SIMEX.
#' \emph{Biometrics}, 62(1), 85--96.
#'
#' Sevilimedu, V. and Yu, L. (2026). An improved misclassification simulation
#' extrapolation (MC-SIMEX) algorithm. \emph{Statistics in Medicine}, 45,
#' e70418.
#'
#' @examples
#' # --- Continuous measurement error ---
#' set.seed(42)
#' n <- 500
#' x_true <- rnorm(n)
#' y <- 1 + 2 * x_true + rnorm(n, sd = 0.5)
#' x_obs <- x_true + rnorm(n, sd = 0.5)
#' df <- data.frame(y = y, x = x_obs)
#'
#' fit_me <- simex(y ~ me(x, 0.5), data = df, B = 50)
#' summary(fit_me)
#'
#' # --- Misclassification (improved, default) ---
#' Pi <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
#' z <- rbinom(n, 1, 0.4)
#' y2 <- rpois(n, exp(0.5 + 0.8 * z + 0.3 * x_true))
#' z_star <- z
#' z_star[z == 0] <- rbinom(sum(z == 0), 1, 0.1)
#' z_star[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)
#' df2 <- data.frame(y = y2, z = factor(z_star), x = x_true)
#'
#' fit_mc <- simex(y ~ mc(z, Pi) + x, family = poisson(), data = df2)
#' summary(fit_mc)
#'
#' @export
simex <- function(formula, family = gaussian(), data,
                  method = NULL,
                  lambda = NULL, B = NULL,
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

  # --- parse formula ---
  parsed <- parse_simex_formula(formula, data, parent.frame())

  # --- dispatch ---
  if (parsed$error_type == "me") {
    if (!is.null(method))
      warning("'method' is ignored for continuous measurement error (me() terms).",
              call. = FALSE)
    # defaults for me
    if (is.null(lambda)) lambda <- c(0.5, 1, 1.5, 2)
    if (is.null(B)) B <- 200L
    result <- .simex_continuous(parsed, family, data, lambda, B,
                                extrapolation, jackknife, weights, seed, cl)
  } else {
    # mc terms
    n_mc <- length(parsed$mc_terms)
    has_response_mc <- !is.null(parsed$response_mc)

    if (is.null(method)) {
      method <- if (n_mc == 1L && !has_response_mc) "improved" else "standard"
    }
    method <- match.arg(method, c("standard", "improved"))

    if (method == "improved" && (n_mc > 1L || has_response_mc))
      stop("method = 'improved' is only supported for a single mc() covariate ",
           "term without response misclassification.", call. = FALSE)

    # defaults depend on method
    if (is.null(lambda)) {
      lambda <- if (method == "improved") 1 else c(0.5, 1, 1.5, 2)
    }
    if (is.null(B)) {
      B <- if (method == "improved") 1L else 200L
    }
    result <- .simex_discrete(parsed, family, data, method, lambda, B,
                              extrapolation, jackknife, weights, seed, cl)
  }

  # --- attach common fields ---
  result$formula <- formula
  result$clean.formula <- parsed$clean_formula
  result$call <- cl
  result$error.type <- parsed$error_type
  result$me.terms <- parsed$me_terms
  result$mc.terms <- parsed$mc_terms
  result$response.mc <- parsed$response_mc

  structure(result, class = "simex")
}


# =========================================================================
# Internal: continuous measurement error fitting
# =========================================================================

.simex_continuous <- function(parsed, family, data, lambda, B,
                              extrapolation, jackknife, weights, seed, cl) {
  dist_code <- family_to_dist_code(family)
  me_terms <- parsed$me_terms
  clean_formula <- parsed$clean_formula

  SIMEXvariable <- vapply(me_terms, `[[`, character(1), "variable")

  # Fit naive model — use .wt in data to avoid glm() NSE issues with weights
  if (!is.null(weights)) {
    data$.wt <- as.numeric(weights)
    naive_fit <- glm(clean_formula, family = family, data = data,
                     weights = .wt, x = TRUE, y = TRUE)
  } else {
    naive_fit <- glm(clean_formula, family = family, data = data,
                     x = TRUE, y = TRUE)
  }

  y <- naive_fit$y
  X <- naive_fit$x
  n <- length(y)
  p <- ncol(X)
  p_names <- colnames(X)
  psi_naive <- unname(coef(naive_fit))
  names(psi_naive) <- p_names
  wt <- if (is.null(weights)) rep(1.0, n) else as.numeric(weights)

  lambda <- sort(as.numeric(lambda))
  B <- as.integer(B)
  n_lambda <- length(lambda)
  n_simex <- length(SIMEXvariable)

  # Build measurement error matrix and per-variable parameters from me_terms
  measurement_error <- matrix(NA_real_, nrow = n, ncol = n_simex)
  me_mean <- numeric(n_simex)
  me_error_type <- integer(n_simex)  # 0 = classical, 1 = berkson
  simex_cols <- integer(n_simex)
  for (j in seq_len(n_simex)) {
    sv <- SIMEXvariable[j]
    sd_val <- me_terms[[j]]$sd
    if (length(sd_val) == 1) sd_val <- rep(sd_val, n)
    if (length(sd_val) != n)
      stop("me() sd for '", sv, "' must be scalar or length n.", call. = FALSE)
    measurement_error[, j] <- sd_val

    me_mean[j] <- if (!is.null(me_terms[[j]]$mean)) me_terms[[j]]$mean else 0
    me_error_type[j] <- if (identical(me_terms[[j]]$error_type, "berkson")) 1L else 0L

    idx <- which(p_names == sv)
    if (length(idx) != 1)
      stop("me() variable '", sv,
           "' must appear as exactly one untransformed main effect.",
           call. = FALSE)
    simex_cols[j] <- idx - 1L
  }

  # C++ simulation (returns list with theta matrix and per-lambda avg vcov)
  sim_result <- simex_sim_cpp(y, X, simex_cols, measurement_error,
                              dist_code, lambda, B, wt, as.integer(seed),
                              me_error_type, me_mean)
  sim_theta <- sim_result$theta
  sim_vcov_model <- sim_result$vcov_model

  # Reshape
  theta_list <- vector("list", n_lambda)
  avg_estimates <- matrix(0, n_lambda, p)
  for (l in seq_len(n_lambda)) {
    rows <- ((l - 1) * B + 1):(l * B)
    theta_mat <- sim_theta[rows, , drop = FALSE]
    colnames(theta_mat) <- p_names
    theta_list[[l]] <- theta_mat
    avg_estimates[l, ] <- colMeans(theta_mat)
  }
  names(theta_list) <- paste0("lambda_", lambda)

  # Reshape per-lambda model vcov averages
  vcov_model_list <- vector("list", n_lambda)
  for (l in seq_len(n_lambda)) {
    vcov_model_list[[l]] <- matrix(sim_vcov_model[l, ], p, p)
  }

  all_lambda <- c(0, lambda)
  all_estimates <- rbind(psi_naive, avg_estimates)
  rownames(all_estimates) <- paste0("lambda_", all_lambda)
  colnames(all_estimates) <- p_names

  # Extrapolation
  extrap_result <- extrapolate_simex(all_lambda, all_estimates, extrapolation)

  # Jackknife variance
  vcov_jk <- NULL
  if (!isFALSE(jackknife)) {
    jk_method <- if (is.character(jackknife)) jackknife else extrapolation
    vcov_jk <- jackknife_variance_simex(theta_list, naive_fit, lambda,
                                        jk_method, vcov_model_list)
    rownames(vcov_jk) <- colnames(vcov_jk) <- p_names
  }

  corrected_coefs <- extrap_result$coefficients
  eta <- as.numeric(X %*% corrected_coefs)
  fitted_vals <- family$linkinv(eta)

  list(
    coefficients        = corrected_coefs,
    naive.coefficients  = psi_naive,
    residuals           = y - fitted_vals,
    fitted.values       = fitted_vals,
    family              = family,
    data                = data,
    vcov                = vcov_jk,
    SIMEX.estimates     = all_estimates,
    theta               = theta_list,
    vcov.model          = vcov_model_list,
    extrapolation       = extrap_result$model,
    lambda              = all_lambda,
    B                   = B,
    n                   = n,
    p                   = p,
    naive.model         = naive_fit,
    extrapolation.method = extrapolation,
    method              = NULL
  )
}


# =========================================================================
# Internal: discrete misclassification fitting
# =========================================================================

.simex_discrete <- function(parsed, family, data, method, lambda, B,
                            extrapolation, jackknife, weights, seed, cl) {
  dist_code <- family_to_dist_code(family)
  mc_terms <- parsed$mc_terms
  n_mc <- length(mc_terms)
  clean_formula <- parsed$clean_formula
  response_mc <- parsed$response_mc
  has_response_mc <- !is.null(response_mc)

  SIMEXvariables <- vapply(mc_terms, `[[`, character(1), "variable")

  # Build named list of misclassification matrices
  mc_list <- setNames(
    lapply(mc_terms, `[[`, "mc_matrix"),
    SIMEXvariables
  )

  # Ensure all mc variables are factors

  for (sv in SIMEXvariables) {
    if (!is.factor(data[[sv]]))
      data[[sv]] <- factor(data[[sv]])
  }

  # Single-mc: extract for C++ fast path
  if (n_mc == 1L) {
    Pi <- mc_list[[1]]
    K <- nrow(Pi)
    SIMEXvariable <- SIMEXvariables[1]
  }

  # Improved method currently handles one mc covariate. The correction itself
  # supports K-level factors; response and multi-covariate extensions still use
  # standard MC-SIMEX.
  if (method == "improved") {
    if (n_mc > 1L)
      stop("method = 'improved' is only supported for a single mc() term.",
           call. = FALSE)
  }

  # Handle lambda = "optimal"
  optimal_lambda <- FALSE
  if (is.character(lambda) && lambda == "optimal") {
    optimal_lambda <- TRUE
    lambda <- 1  # placeholder
  }

  # Fit naive model — use .wt in data to avoid glm() NSE issues with weights
  if (!is.null(weights)) {
    data$.wt <- as.numeric(weights)
    naive_fit <- glm(clean_formula, family = family, data = data,
                     weights = .wt, x = TRUE, y = TRUE)
  } else {
    naive_fit <- glm(clean_formula, family = family, data = data,
                     x = TRUE, y = TRUE)
  }

  y <- naive_fit$y
  n <- length(y)
  wt <- if (is.null(weights)) rep(1.0, n) else as.numeric(weights)
  B <- as.integer(B)
  fam <- get_link_funs(family)

  # --- Build design matrix and naive refit ---
  if (has_response_mc && n_mc == 0L) {
    # Response-only mc: use naive GLM's design matrix
    psi_naive <- unname(coef(naive_fit))
    nms <- names(coef(naive_fit))
    names(psi_naive) <- nms
    p <- length(psi_naive)
    xi_hat <- naive_fit$x
    naive_refit <- naive_fit
  } else if (has_response_mc && n_mc >= 1L) {
    # Response mc + covariate mc(s): build covariate design matrix for C++
    z_hats <- list()
    z_levels_list <- list()
    K_vec <- integer(n_mc)
    for (j in seq_len(n_mc)) {
      sv <- SIMEXvariables[j]
      z_factor_j <- data[[sv]]
      z_levels_list[[sv]] <- levels(z_factor_j)
      z_hats[[j]] <- as.integer(z_factor_j) - 1L
      K_vec[j] <- length(levels(z_factor_j))
    }

    xm <- .build_x_mat(clean_formula, SIMEXvariables, data)
    x_mat <- xm$x_mat

    xi_hat <- .build_multi_xi_hat(z_hats, K_vec, x_mat)
    p <- ncol(xi_hat)
    nms <- .make_multi_param_names(z_levels_list, K_vec, x_mat)

    dat_naive <- data.frame(y = y, .wt = wt, xi_hat)
    naive_refit <- glm(y ~ . - .wt - 1, data = dat_naive, family = family,
                       weights = .wt)
    psi_naive <- unname(coef(naive_refit))
    names(psi_naive) <- nms
  } else if (n_mc == 1L) {
    # Single-mc covariate: manual design matrix (compatible with C++ path)
    z_factor <- data[[SIMEXvariable]]
    z_levels <- levels(z_factor)
    z_levels_list <- setNames(list(z_levels), SIMEXvariable)
    z_hat <- as.integer(z_factor) - 1L

    xm <- .build_x_mat(clean_formula, SIMEXvariable, data)
    x_mat <- xm$x_mat

    xi_hat <- .build_xi_hat(z_hat, x_mat, K)
    p <- ncol(xi_hat)
    nms <- .make_param_names(z_levels, x_mat, K)

    dat_naive <- data.frame(y = y, .wt = wt, xi_hat)
    naive_refit <- glm(y ~ . - .wt - 1, data = dat_naive, family = family,
                       weights = .wt)
    psi_naive <- unname(coef(naive_refit))
    names(psi_naive) <- nms
  } else {
    # Multi-mc: build design matrix manually [dummies_z1 | dummies_z2 | ... | x]
    z_hats <- list()
    z_levels_list <- list()
    K_vec <- integer(n_mc)
    for (j in seq_len(n_mc)) {
      sv <- SIMEXvariables[j]
      z_factor_j <- data[[sv]]
      z_levels_list[[sv]] <- levels(z_factor_j)
      z_hats[[j]] <- as.integer(z_factor_j) - 1L
      K_vec[j] <- length(levels(z_factor_j))
    }

    xm <- .build_x_mat(clean_formula, SIMEXvariables, data)
    x_mat <- xm$x_mat

    xi_hat <- .build_multi_xi_hat(z_hats, K_vec, x_mat)
    p <- ncol(xi_hat)
    nms <- .make_multi_param_names(z_levels_list, K_vec, x_mat)

    dat_naive <- data.frame(y = y, .wt = wt, xi_hat)
    naive_refit <- glm(y ~ . - .wt - 1, data = dat_naive, family = family,
                       weights = .wt)
    psi_naive <- unname(coef(naive_refit))
    names(psi_naive) <- nms
  }

  # --- IMPROVED MC-SIMEX (single mc only) ---
  if (method == "improved") {
    pi_vec <- .estimate_pi_vec(z_hat, Pi, wt)
    pi_x <- if (K == 2L) pi_vec[2] else NULL

    if (optimal_lambda) {
      if (K != 2L)
        stop("lambda = 'optimal' is currently only supported for binary ",
             "improved MC-SIMEX.", call. = FALSE)
      lambda <- .find_optimal_lambda(Pi, pi_x)
    }
    lambda <- sort(as.numeric(lambda))
    n_lambda <- length(lambda)

    correction_components <- .compute_k_correction(Pi, pi_vec, lambda)
    c_lam_vec <- if (K == 2L) {
      vapply(correction_components, function(cc) cc$C[1, 1], numeric(1))
    } else {
      NULL
    }

    sim_result <- mcsimex_sim_cpp(y, z_hat, x_mat, Pi, K, dist_code,
                                  lambda, B, wt, as.integer(seed))
    s <- K - 1L
    intercept_idx <- which(nms == "(Intercept)")
    if (length(intercept_idx) != 1L) intercept_idx <- integer(0)
    theta_list <- vector("list", n_lambda)
    corrected_all <- matrix(0, nrow = 0, ncol = p)

    for (l in seq_len(n_lambda)) {
      rows <- ((l - 1) * B + 1):(l * B)
      theta_mat <- sim_result[rows, , drop = FALSE]
      colnames(theta_mat) <- nms
      theta_corrected <- .apply_k_improved_correction(
        theta_mat,
        correction_components[[l]]$C,
        correction_components[[l]]$alpha,
        intercept_idx
      )
      colnames(theta_corrected) <- nms
      theta_list[[l]] <- theta_corrected
      corrected_all <- rbind(corrected_all, theta_corrected)
    }
    names(theta_list) <- paste0("lambda_", lambda)

    corrected_coefs <- colMeans(corrected_all)
    names(corrected_coefs) <- nms

    avg_estimates <- matrix(0, n_lambda, p)
    for (l in seq_len(n_lambda)) {
      avg_estimates[l, ] <- colMeans(theta_list[[l]])
    }
    all_lambda <- c(0, lambda)
    all_estimates <- rbind(psi_naive, avg_estimates)
    rownames(all_estimates) <- paste0("lambda_", all_lambda)
    colnames(all_estimates) <- nms

    vcov_imp <- NULL
    if (!isFALSE(jackknife)) {
      transform <- .k_improved_transform(
        p,
        correction_components[[1]]$C,
        correction_components[[1]]$alpha,
        intercept_idx
      )
      vcov_imp <- .variance_k_improved(theta_list, naive_refit, transform,
                                       lambda, B, p)
      rownames(vcov_imp) <- colnames(vcov_imp) <- nms
    }

    eta <- as.numeric(xi_hat %*% corrected_coefs)
    fitted_vals <- family$linkinv(eta)

    return(list(
      coefficients        = corrected_coefs,
      naive.coefficients  = psi_naive,
      residuals           = y - fitted_vals,
      fitted.values       = fitted_vals,
      family              = family,
      data                = data,
      vcov                = vcov_imp,
      SIMEX.estimates     = all_estimates,
      theta               = theta_list,
      extrapolation       = NULL,
      lambda              = all_lambda,
      B                   = B,
      n                   = n,
      p                   = p,
      naive.model         = naive_fit,
      extrapolation.method = "exact (improved)",
      method              = "improved",
      c.lambda            = if (!is.null(c_lam_vec))
        setNames(c_lam_vec, paste0("lambda_", lambda)) else NULL,
      pi.x                = pi_x,
      pi.vec              = pi_vec,
      correction.matrix   = setNames(lapply(correction_components, `[[`, "C"),
                                     paste0("lambda_", lambda)),
      SIMEXvariable       = SIMEXvariables,
      mc.matrix           = mc_list,
      xi.hat              = xi_hat,
      x.mat               = x_mat,
      z.levels            = z_levels_list
    ))
  }

  # --- STANDARD MC-SIMEX ---
  lambda <- sort(as.numeric(lambda))
  n_lambda <- length(lambda)

  set.seed(seed)

  if (has_response_mc) {
    # Response misclassification: C++ path via mcsimex_multi_sim_cpp
    resp_var <- response_mc$variable
    if (!is.factor(data[[resp_var]]))
      data[[resp_var]] <- factor(data[[resp_var]])
    y_z_hat <- as.integer(data[[resp_var]]) - 1L
    Pi_y <- response_mc$mc_matrix
    K_y <- nrow(Pi_y)

    if (n_mc == 0L) {
      # Response-only mc: no covariate mc terms
      sim_result <- mcsimex_multi_sim_cpp(
        y, list(), list(), integer(0), xi_hat,
        dist_code, lambda, B, wt, as.integer(seed),
        y_z_hat, Pi_y, K_y
      )
    } else {
      # Response mc + covariate mc(s)
      Pi_list_cov <- lapply(mc_list, unname)
      sim_result <- mcsimex_multi_sim_cpp(
        y, z_hats, Pi_list_cov, K_vec, x_mat,
        dist_code, lambda, B, wt, as.integer(seed),
        y_z_hat, Pi_y, K_y
      )
    }

    theta_list <- vector("list", n_lambda)
    avg_estimates <- matrix(0, n_lambda, p)
    for (l in seq_len(n_lambda)) {
      rows <- ((l - 1) * B + 1):(l * B)
      theta_mat <- sim_result[rows, , drop = FALSE]
      colnames(theta_mat) <- nms
      theta_list[[l]] <- theta_mat
      avg_estimates[l, ] <- colMeans(theta_mat)
    }
  } else if (n_mc == 1L) {
    # Single mc covariate: C++ fast path
    sim_result <- mcsimex_sim_cpp(y, z_hat, x_mat, Pi, K, dist_code,
                                  lambda, B, wt, as.integer(seed))

    theta_list <- vector("list", n_lambda)
    avg_estimates <- matrix(0, n_lambda, p)
    for (l in seq_len(n_lambda)) {
      rows <- ((l - 1) * B + 1):(l * B)
      theta_mat <- sim_result[rows, , drop = FALSE]
      colnames(theta_mat) <- nms
      theta_list[[l]] <- theta_mat
      avg_estimates[l, ] <- colMeans(theta_mat)
    }
  } else {
    # Multi-mc covariates: C++ fast path
    Pi_list <- lapply(mc_list, unname)
    sim_result <- mcsimex_multi_sim_cpp(y, z_hats, Pi_list, K_vec, x_mat,
                                        dist_code, lambda, B, wt,
                                        as.integer(seed),
                                        integer(0),
                                        matrix(0, nrow = 0, ncol = 0),
                                        0L)

    theta_list <- vector("list", n_lambda)
    avg_estimates <- matrix(0, n_lambda, p)
    for (l in seq_len(n_lambda)) {
      rows <- ((l - 1) * B + 1):(l * B)
      theta_mat <- sim_result[rows, , drop = FALSE]
      colnames(theta_mat) <- nms
      theta_list[[l]] <- theta_mat
      avg_estimates[l, ] <- colMeans(theta_mat)
    }
  }
  names(theta_list) <- paste0("lambda_", lambda)

  all_lambda <- c(0, lambda)
  all_estimates <- rbind(psi_naive, avg_estimates)
  rownames(all_estimates) <- paste0("lambda_", all_lambda)
  colnames(all_estimates) <- nms

  extrap_result <- extrapolate_simex(all_lambda, all_estimates, extrapolation)

  vcov_jk <- NULL
  if (!isFALSE(jackknife)) {
    jk_method <- if (is.character(jackknife)) jackknife else extrapolation
    vcov_jk <- jackknife_variance_mcsimex(theta_list, psi_naive, naive_refit,
                                          lambda, jk_method, fam, xi_hat, y, wt)
    rownames(vcov_jk) <- colnames(vcov_jk) <- nms
  }

  corrected_coefs <- extrap_result$coefficients
  eta <- as.numeric(xi_hat %*% corrected_coefs)
  fitted_vals <- family$linkinv(eta)

  list(
    coefficients        = corrected_coefs,
    naive.coefficients  = psi_naive,
    residuals           = y - fitted_vals,
    fitted.values       = fitted_vals,
    family              = family,
    data                = data,
    vcov                = vcov_jk,
    SIMEX.estimates     = all_estimates,
    theta               = theta_list,
    extrapolation       = extrap_result$model,
    lambda              = all_lambda,
    B                   = B,
    n                   = n,
    p                   = p,
    naive.model         = naive_fit,
    extrapolation.method = extrapolation,
    method              = "standard",
    SIMEXvariable       = SIMEXvariables,
    mc.matrix           = mc_list,
    xi.hat              = xi_hat,
    x.mat               = if (exists("x_mat")) x_mat else NULL,
    z.levels            = if (exists("z_levels_list")) z_levels_list else NULL
  )
}


# =========================================================================
# S3 methods for class "simex"
# =========================================================================

#' @export
print.simex <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall:\n")
  print(x$call)
  label <- if (x$error.type == "mc") "MC-SIMEX" else "SIMEX"
  cat(sprintf("\n%s corrected coefficients:\n", label))
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
  ans$error.type <- object$error.type
  ans$me.terms <- object$me.terms
  ans$mc.terms <- object$mc.terms
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
  ans$pi.vec <- object$pi.vec

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

  if (x$error.type == "me") {
    me_vars <- vapply(x$me.terms, `[[`, character(1), "variable")
    cat("SIMEX variable(s):", paste(me_vars, collapse = ", "), "\n")
  } else {
    mc_vars <- vapply(x$mc.terms, `[[`, character(1), "variable")
    label_mc <- if (length(mc_vars) > 1) "MC-SIMEX variables:" else "MC-SIMEX variable:"
    cat(label_mc, paste(mc_vars, collapse = ", "), "\n")
    if (!is.null(x$method))
      cat("Method:", x$method, "\n")
  }

  cat("Extrapolation:", x$extrapolation.method, "\n")
  cat("Lambda grid:", paste(x$lambda, collapse = ", "), "\n")
  cat("B =", x$B, ", n =", x$n, "\n")

  if (!is.null(x$pi.x))
    cat("Estimated P(X=1):", round(x$pi.x, 4), "\n")
  else if (!is.null(x$pi.vec))
    cat("Estimated latent probabilities:",
        paste(round(x$pi.vec, 4), collapse = ", "), "\n")
  if (!is.null(x$c.lambda))
    cat("Correction factor(s):", paste(round(x$c.lambda, 4), collapse = ", "), "\n")

  if (!is.null(x$residuals)) {
    cat("\nResiduals:\n")
    rq <- quantile(x$residuals, na.rm = TRUE)
    names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
    print(rq, digits = digits)
  }

  cat("\nNaive coefficients:\n")
  print(round(x$naive.coefficients, digits))

  label <- if (x$error.type == "mc") "MC-SIMEX" else "SIMEX"
  cat(sprintf("\n%s corrected coefficients:\n", label))
  if (ncol(x$coefficients) > 1) {
    printCoefmat(x$coefficients, digits = digits)
  } else {
    print(round(x$coefficients, digits))
  }
  cat("\n")
  invisible(x)
}

#' @export
plot.simex <- function(x, xlab = expression(1 + lambda),
                       ylab = NULL, ask = FALSE, show = NULL, ...) {
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
      plot(lam + 1, b[, j], xlab = xlab, ylab = ylab[j], type = "n",
           main = p.names[j], ...)
      points(lam[-1] + 1, b[-1, j], pch = 19)
      points(lam[1] + 1, b[1, j], pch = 1)
      abline(h = x$coefficients[j], lty = 2, col = "red")
      points(0, x$coefficients[j], pch = 4, cex = 1.5, col = "red")
    }
  }
}

#' @export
predict.simex <- function(object, newdata = NULL, type = "response", ...) {
  if (is.null(newdata)) {
    if (type == "link") return(object$family$linkfun(object$fitted.values))
    return(object$fitted.values)
  }

  mc_vars <- object$SIMEXvariable
  if (object$error.type == "mc" && length(mc_vars) >= 1L) {
    z_levels <- object$z.levels

    # Enforce training factor levels for all mc variables
    for (sv in mc_vars) {
      if (!is.null(z_levels) && sv %in% names(z_levels)) {
        newdata[[sv]] <- factor(newdata[[sv]], levels = z_levels[[sv]])
      } else if (!is.factor(newdata[[sv]])) {
        newdata[[sv]] <- factor(newdata[[sv]])
      }
    }

    if (length(mc_vars) == 1L) {
      sv <- mc_vars[1]
      z_new <- as.integer(newdata[[sv]]) - 1L
      K <- nrow(object$mc.matrix[[sv]])
      xm <- .build_x_mat(object$clean.formula, sv, newdata)
      X_new <- .build_xi_hat(z_new, xm$x_mat, K)
    } else {
      # Multi-mc: use same manual design matrix construction as fitting
      z_hats_new <- list()
      K_vec <- integer(length(mc_vars))
      for (j in seq_along(mc_vars)) {
        sv <- mc_vars[j]
        z_hats_new[[j]] <- as.integer(newdata[[sv]]) - 1L
        K_vec[j] <- nrow(object$mc.matrix[[sv]])
      }
      xm <- .build_x_mat(object$clean.formula, mc_vars, newdata)
      X_new <- .build_multi_xi_hat(z_hats_new, K_vec, xm$x_mat)
    }
  } else {
    # ME predict: use model.matrix from clean formula
    tt <- delete.response(terms(object$clean.formula))
    mf <- model.frame(tt, data = newdata)
    X_new <- model.matrix(tt, mf)
  }

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
  if (missing(parm)) parm <- pnames
  else if (is.numeric(parm)) parm <- pnames[parm]

  ci <- cbind(cf[parm] + fac[1] * ses[parm],
              cf[parm] + fac[2] * ses[parm])
  colnames(ci) <- paste0(format(100 * c(a, 1 - a), trim = TRUE), " %")
  ci
}
