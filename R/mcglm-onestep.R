# ---------------------------------------------------------------------------
# One-step joint estimation via RTMB mixture likelihood
# ---------------------------------------------------------------------------

#' Compute mixture weights from known misclassification rates (binary)
#' @keywords internal
.mcglm_compute_omega_bin <- function(p01, p10, pi_z) {
  c(
    (1 - p01) * (1 - pi_z),
    p10 * pi_z,
    p01 * (1 - pi_z),
    (1 - p10) * pi_z
  )
}

#' Compute mixture weights from known misclassification matrix (multicategory)
#' @keywords internal
.mcglm_compute_omega_multi <- function(Pi, pi_z, K) {
  omega <- numeric(K * K)
  for (j in seq_len(K)) {
    for (l in seq_len(K)) {
      omega[(j - 1) * K + l] <- Pi[j, l] * pi_z[l]
    }
  }
  omega
}


#' One-step joint estimator for binary misclassification (GLM)
#' @keywords internal
.mcglm_fit_onestep_bin <- function(y, z_hat, x, family,
                                   p01 = NULL, p10 = NULL, pi_z = NULL,
                                   fix_omega = FALSE, homoskedastic = TRUE,
                                   optim_control = list(), wt = NULL) {
  if (!requireNamespace("RTMB", quietly = TRUE))
    stop("Package 'RTMB' is required for the one-step method. ",
         "Install with: install.packages('RTMB')")

  fam <- .mcglm_get_link_funs(family)
  n   <- length(y)
  r   <- ncol(x)
  K   <- 2L
  s   <- 1L
  d   <- s + r

  dist_code <- switch(fam$family$family,
    gaussian = 1L, poisson = 2L, binomial = 3L,
    stop("Unsupported family for one-step: ", fam$family$family)
  )

  weights_fixed <- as.integer(fix_omega)
  if (weights_fixed) {
    stopifnot(!is.null(p01), !is.null(p10), !is.null(pi_z))
    omega_data <- .mcglm_compute_omega_bin(p01, p10, pi_z)
  } else {
    omega_data <- rep(0.25, 4)
  }

  if (is.null(wt)) wt_data <- rep(1.0, n) else wt_data <- as.numeric(wt)

  xi_hat <- cbind(z_hat, x)
  naive_fit <- stats::glm(y ~ . - 1,
    data = data.frame(y = y, xi_hat),
    family = fam$family, weights = wt_data
  )
  b_init <- unname(stats::coef(naive_fit))

  theta_init <- b_init
  if (!weights_fixed) theta_init <- c(theta_init, rep(0, K * K - 1))
  if (dist_code == 1L) {
    resid_sd <- stats::sd(stats::residuals(naive_fit))
    if (homoskedastic) {
      theta_init <- c(theta_init, log(resid_sd))
    } else {
      theta_init <- c(theta_init, log(resid_sd), log(resid_sd))
    }
  }

  nll_fn <- .mcglm_make_nll_glm(
    Y = y, X = x, z_hat = as.integer(z_hat), K = K,
    dist_code = dist_code, omega_data = omega_data,
    wt = wt_data, weights_fixed = weights_fixed,
    homoskedastic = as.integer(homoskedastic)
  )

  obj <- RTMB::MakeADFun(
    func = nll_fn,
    parameters = list(theta = theta_init),
    silent = TRUE
  )

  opt <- .mcglm_run_nlminb(obj, optim_control)

  b_hat <- opt$par[1:d]
  V     <- .mcglm_vcov_onestep(obj, opt, d)
  nms   <- c("gamma", paste0("alpha", seq_len(r) - 1))
  names(b_hat) <- nms
  colnames(V) <- rownames(V) <- nms

  list(coefficients = b_hat, vcov = V,
       loglik = -opt$objective, convergence = opt$convergence)
}


#' One-step joint estimator for multicategory misclassification (GLM)
#' @keywords internal
.mcglm_fit_onestep_multi <- function(y, z_hat, x, K, family,
                                     Pi = NULL, pi_z = NULL,
                                     fix_omega = FALSE, homoskedastic = TRUE,
                                     optim_control = list(), wt = NULL) {
  if (!requireNamespace("RTMB", quietly = TRUE))
    stop("Package 'RTMB' is required for the one-step method. ",
         "Install with: install.packages('RTMB')")

  fam <- .mcglm_get_link_funs(family)
  n   <- length(y)
  r   <- ncol(x)
  s   <- K - 1
  d   <- s + r

  dist_code <- switch(fam$family$family,
    gaussian = 1L, poisson = 2L, binomial = 3L,
    stop("Unsupported family for one-step: ", fam$family$family)
  )

  weights_fixed <- as.integer(fix_omega)
  if (weights_fixed) {
    stopifnot(!is.null(Pi), !is.null(pi_z))
    omega_data <- .mcglm_compute_omega_multi(Pi, pi_z, K)
  } else {
    omega_data <- rep(1 / (K * K), K * K)
  }

  if (is.null(wt)) wt_data <- rep(1.0, n) else wt_data <- as.numeric(wt)

  d_hat <- matrix(0, n, s)
  for (k in seq_len(s)) d_hat[, k] <- as.numeric(z_hat == k)
  xi_hat <- cbind(d_hat, x)
  naive_fit <- stats::glm(y ~ . - 1,
    data = data.frame(y = y, xi_hat),
    family = fam$family, weights = wt_data
  )
  b_init <- unname(stats::coef(naive_fit))

  theta_init <- b_init
  if (!weights_fixed) theta_init <- c(theta_init, rep(0, K * K - 1))
  if (dist_code == 1L) {
    resid_sd <- stats::sd(stats::residuals(naive_fit))
    if (homoskedastic) {
      theta_init <- c(theta_init, log(resid_sd))
    } else {
      theta_init <- c(theta_init, log(resid_sd), log(resid_sd))
    }
  }

  nll_fn <- .mcglm_make_nll_glm(
    Y = y, X = x, z_hat = as.integer(z_hat), K = K,
    dist_code = dist_code, omega_data = omega_data,
    wt = wt_data, weights_fixed = weights_fixed,
    homoskedastic = as.integer(homoskedastic)
  )

  obj <- RTMB::MakeADFun(
    func = nll_fn,
    parameters = list(theta = theta_init),
    silent = TRUE
  )

  opt <- .mcglm_run_nlminb(obj, optim_control)

  b_hat <- opt$par[1:d]
  V     <- .mcglm_vcov_onestep(obj, opt, d)
  nms   <- c(paste0("gamma", seq_len(s)), paste0("alpha", seq_len(r) - 1))
  names(b_hat) <- nms
  colnames(V) <- rownames(V) <- nms

  list(coefficients = b_hat, vcov = V,
       loglik = -opt$objective, convergence = opt$convergence)
}


# ======================== MULTINOMIAL ONE-STEP =============================

#' One-step joint estimator for multinomial logistic with misclassified covariate
#' @keywords internal
.mcglm_fit_onestep_multinomial <- function(y, z_hat, x, J, K,
                                           p01 = NULL, p10 = NULL, pi_z = NULL,
                                           Pi = NULL,
                                           fix_omega = FALSE,
                                           optim_control = list(), wt = NULL) {
  if (!requireNamespace("RTMB", quietly = TRUE))
    stop("Package 'RTMB' is required for the one-step method. ",
         "Install with: install.packages('RTMB')")

  n <- length(y)
  r <- ncol(x)
  s <- K - 1
  block_size <- s + r
  d <- (J - 1) * block_size

  is_binary_z <- (K == 2L)

  weights_fixed <- as.integer(fix_omega)
  if (weights_fixed) {
    if (is_binary_z) {
      stopifnot(!is.null(p01), !is.null(p10), !is.null(pi_z))
      omega_data <- .mcglm_compute_omega_bin(p01, p10, pi_z)
    } else {
      stopifnot(!is.null(Pi), !is.null(pi_z))
      omega_data <- .mcglm_compute_omega_multi(Pi, pi_z, K)
    }
  } else {
    omega_data <- rep(1 / (K * K), K * K)
  }

  if (is.null(wt)) wt_data <- rep(1.0, n) else wt_data <- as.numeric(wt)

  # Starting values via category-specific binomial regressions
  theta_init <- numeric(d)
  d_hat <- matrix(0, n, s)
  for (k in seq_len(s)) d_hat[, k] <- as.numeric(z_hat == k)
  xi_naive <- cbind(d_hat, x)

  for (jj in seq_len(J - 1)) {
    y_bin <- as.numeric(y == jj)
    tryCatch({
      fit_j <- stats::glm(y_bin ~ . - 1,
        data = data.frame(y_bin = y_bin, xi_naive),
        family = stats::binomial(), weights = wt_data
      )
      theta_init[((jj - 1) * block_size + 1):(jj * block_size)] <-
        unname(stats::coef(fit_j))
    }, error = function(e) NULL)
  }

  if (!weights_fixed) theta_init <- c(theta_init, rep(0, K * K - 1))

  nll_fn <- .mcglm_make_nll_multinomial(
    Y_int = as.integer(y), X = x, z_hat = as.integer(z_hat),
    J = as.integer(J), K = K, omega_data = omega_data,
    wt = wt_data, weights_fixed = weights_fixed
  )

  obj <- RTMB::MakeADFun(
    func = nll_fn,
    parameters = list(theta = theta_init),
    silent = TRUE
  )

  opt <- .mcglm_run_nlminb(obj, optim_control)

  b_hat <- opt$par[1:d]
  V     <- .mcglm_vcov_onestep(obj, opt, d)

  nms <- .mcglm_make_multinomial_names(J, K, ncol(x))
  names(b_hat) <- nms
  colnames(V) <- rownames(V) <- nms

  list(coefficients = b_hat, vcov = V,
       loglik = -opt$objective, convergence = opt$convergence)
}


#' Fit naive multinomial logistic regression with proxy covariate
#' @keywords internal
.mcglm_fit_naive_multinomial <- function(y, z_hat, x, J, K, wt = NULL) {
  n <- length(y)
  r <- ncol(x)
  s <- K - 1
  block_size <- s + r
  d <- (J - 1) * block_size

  d_hat <- matrix(0, n, s)
  for (k in seq_len(s)) d_hat[, k] <- as.numeric(z_hat == k)
  xi <- cbind(d_hat, x)

  if (is.null(wt)) wt_data <- rep(1, n) else wt_data <- wt

  if (requireNamespace("nnet", quietly = TRUE)) {
    dat <- data.frame(y = factor(y, levels = 0:(J - 1)), xi)
    fit <- nnet::multinom(y ~ . - 1, data = dat, weights = wt_data, trace = FALSE)
    coef_mat <- stats::coef(fit)
    if (!is.matrix(coef_mat)) coef_mat <- matrix(coef_mat, nrow = 1)
    psi <- as.numeric(t(coef_mat))
    return(list(coefficients = psi, multinom_fit = fit))
  }

  psi <- numeric(d)
  for (jj in seq_len(J - 1)) {
    y_bin <- as.numeric(y == jj)
    tryCatch({
      fit_j <- stats::glm(y_bin ~ . - 1,
        data = data.frame(y_bin = y_bin, xi),
        family = stats::binomial(), weights = wt_data
      )
      psi[((jj - 1) * block_size + 1):(jj * block_size)] <-
        unname(stats::coef(fit_j))
    }, error = function(e) NULL)
  }
  list(coefficients = psi, multinom_fit = NULL)
}


# ======================== SHARED UTILITIES ================================

#' Generate parameter names for multinomial model
#' @keywords internal
.mcglm_make_multinomial_names <- function(J, K, r) {
  s <- K - 1
  block_size <- s + r
  d <- (J - 1) * block_size
  nms <- character(d)
  for (jj in seq_len(J - 1)) {
    offset <- (jj - 1) * block_size
    if (s == 1) {
      nms[offset + 1] <- paste0("y", jj, ":gamma")
    } else if (s > 1) {
      nms[offset + seq_len(s)] <- paste0("y", jj, ":gamma", seq_len(s))
    }
    nms[offset + s + seq_len(r)] <- paste0("y", jj, ":alpha", seq_len(r) - 1)
  }
  nms
}


#' Run nlminb with BFGS fallback
#' @keywords internal
.mcglm_run_nlminb <- function(obj, optim_control = list()) {
  ctrl <- list(iter.max = 500, abs.tol = 1e-12, rel.tol = 1e-10)
  ctrl[names(optim_control)] <- optim_control

  opt <- tryCatch({
    stats::nlminb(
      start = obj$par, objective = obj$fn, gradient = obj$gr,
      control = ctrl
    )
  }, error = function(e) {
    warning("nlminb failed, trying BFGS fallback")
    stats::optim(
      par = obj$par, fn = obj$fn, gr = obj$gr,
      method = "BFGS",
      control = list(maxit = 500, reltol = 1e-10)
    )
  })

  if (!is.null(opt$convergence) && opt$convergence != 0) {
    warning("Optimization may not have converged (code: ", opt$convergence, ")")
  }
  opt
}


#' Compute variance-covariance matrix from Hessian
#' @keywords internal
.mcglm_vcov_onestep <- function(obj, opt, d) {
  H <- tryCatch({
    obj$he(opt$par)
  }, error = function(e) {
    warning("Hessian computation failed, using numerical approximation")
    if (requireNamespace("numDeriv", quietly = TRUE)) {
      numDeriv::hessian(obj$fn, opt$par)
    } else {
      warning("numDeriv not available; returning diagonal variance")
      return(diag(1e-4, length(opt$par)))
    }
  })

  if (any(!is.finite(H))) {
    warning("Hessian contains non-finite values")
    H[!is.finite(H)] <- 0
  }

  eig_vals <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  min_eig  <- min(eig_vals)
  if (min_eig <= 1e-12) {
    warning("Hessian not positive definite (min eigenvalue: ", min_eig,
            "), adding regularization")
    H <- H + diag(abs(min_eig) + 1e-8, nrow(H))
  }

  V_full <- tryCatch({
    chol2inv(chol(H))
  }, error = function(e) {
    tryCatch(solve(H), error = function(e2) {
      if (requireNamespace("MASS", quietly = TRUE)) {
        MASS::ginv(H)
      } else {
        warning("All covariance computations failed, returning diagonal")
        diag(1e-4, nrow(H))
      }
    })
  })

  V <- V_full[1:d, 1:d, drop = FALSE]

  if (any(!is.finite(V))) {
    warning("Covariance matrix contains non-finite values")
    V[!is.finite(V)] <- 0
    diag(V) <- pmax(diag(V), 1e-8)
  }

  V
}
