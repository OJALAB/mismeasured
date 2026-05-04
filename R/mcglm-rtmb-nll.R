# ---------------------------------------------------------------------------
# RTMB negative log-likelihood functions for one-step estimation
# ---------------------------------------------------------------------------

#' Build RTMB NLL function for GLM with misclassified covariate
#' @keywords internal
.mcglm_make_nll_glm <- function(Y, X, z_hat, K, dist_code, omega_data,
                                wt, weights_fixed, homoskedastic) {
  n <- length(Y)
  r <- ncol(X)
  s <- K - 1L
  d <- s + r
  n_omega <- K * K

  omega_idx <- matrix(0L, n, K)
  for (i in seq_len(n)) {
    for (l_r in seq_len(K)) {
      omega_idx[i, l_r] <- z_hat[i] * K + l_r
    }
  }

  if (dist_code == 2L) lfact_Y <- lfactorial(Y)

  function(parms) {
    theta <- parms$theta
    alpha <- theta[(s + 1):(s + r)]

    idx <- d + 1L
    if (weights_fixed == 1L) {
      omega <- omega_data
    } else {
      n_free <- n_omega - 1L
      vraw <- theta[idx:(idx + n_free - 1)]
      idx <- idx + n_free
      expv <- c(exp(vraw), 1)
      omega <- expv / sum(expv)
    }

    sigma0 <- 1; sigma1 <- 1
    if (dist_code == 1L) {
      sigma0 <- exp(theta[idx])
      idx <- idx + 1L
      sigma1 <- if (homoskedastic == 1L) sigma0 else exp(theta[idx])
    }

    eta_base <- X %*% alpha

    log_f_cols <- vector("list", K)

    compute_log_f <- function(eta_col, sig) {
      if (dist_code == 1L) {
        -0.5 * log(2 * pi) - log(sig) - 0.5 * ((Y - eta_col) / sig)^2
      } else if (dist_code == 2L) {
        Y * eta_col - exp(eta_col) - lfact_Y
      } else {
        p_col <- 1 / (1 + exp(-eta_col))
        Y * log(p_col + 1e-20) + (1 - Y) * log(1 - p_col + 1e-20)
      }
    }

    log_f_cols[[1]] <- compute_log_f(eta_base, sigma0)
    for (k in seq_len(s)) {
      sig <- if (dist_code == 1L && homoskedastic != 1L) sigma1 else sigma0
      log_f_cols[[k + 1]] <- compute_log_f(eta_base + theta[k], sig)
    }

    nll <- 0
    for (i in seq_len(n)) {
      w1 <- omega[omega_idx[i, 1]]
      acc <- log(w1 + 1e-300) + log_f_cols[[1]][i]

      for (l_r in 2:K) {
        wl <- omega[omega_idx[i, l_r]]
        if (wl > 1e-20) {
          term <- log(wl) + log_f_cols[[l_r]][i]
          acc <- acc + log(1 + exp(term - acc))
        }
      }
      nll <- nll - wt[i] * acc
    }

    nll
  }
}


#' Build RTMB NLL function for multinomial logistic with misclassified covariate
#' @keywords internal
.mcglm_make_nll_multinomial <- function(Y_int, X, z_hat, J, K, omega_data,
                                        wt, weights_fixed) {
  n <- length(Y_int)
  r <- ncol(X)
  s <- K - 1L
  block_size <- s + r
  d <- (J - 1L) * block_size
  n_omega <- K * K

  omega_idx <- matrix(0L, n, K)
  for (i in seq_len(n)) {
    for (l_r in seq_len(K)) {
      omega_idx[i, l_r] <- z_hat[i] * K + l_r
    }
  }
  yi_r <- Y_int + 1L

  function(parms) {
    theta <- parms$theta

    alpha_x_list <- vector("list", J - 1)
    for (jj in seq_len(J - 1)) {
      offset <- (jj - 1L) * block_size
      alpha_j <- theta[(offset + s + 1):(offset + s + r)]
      alpha_x_list[[jj]] <- X %*% alpha_j
    }

    eta_cols <- vector("list", J)
    eta_cols[[1]] <- rep(list(rep(0, n)), K)

    for (jj in seq_len(J - 1)) {
      offset <- (jj - 1L) * block_size
      ax <- alpha_x_list[[jj]]
      cols <- vector("list", K)
      cols[[1]] <- ax
      for (k in seq_len(s)) {
        cols[[k + 1]] <- ax + theta[offset + k]
      }
      eta_cols[[jj + 1]] <- cols
    }

    if (weights_fixed == 1L) {
      omega <- omega_data
    } else {
      idx <- d + 1L
      n_free <- n_omega - 1L
      vraw <- theta[idx:(idx + n_free - 1)]
      expv <- c(exp(vraw), 1)
      omega <- expv / sum(expv)
    }

    nll <- 0
    for (i in seq_len(n)) {
      yi <- yi_r[i]

      w1 <- omega[omega_idx[i, 1]]

      eta_yi_1 <- eta_cols[[yi]][[1]][i]
      lse_1 <- eta_cols[[1]][[1]][i]
      for (jj in 2:J) {
        eta_jj <- eta_cols[[jj]][[1]][i]
        lse_1 <- lse_1 + log(1 + exp(eta_jj - lse_1))
      }
      log_prob_1 <- eta_yi_1 - lse_1
      acc <- log(w1 + 1e-300) + log_prob_1

      if (K > 1) {
        for (l_r in 2:K) {
          wl <- omega[omega_idx[i, l_r]]
          if (wl > 1e-20) {
            eta_yi_l <- eta_cols[[yi]][[l_r]][i]
            lse_l <- eta_cols[[1]][[l_r]][i]
            for (jj in 2:J) {
              eta_jj <- eta_cols[[jj]][[l_r]][i]
              lse_l <- lse_l + log(1 + exp(eta_jj - lse_l))
            }
            log_prob_l <- eta_yi_l - lse_l
            term <- log(wl) + log_prob_l
            acc <- acc + log(1 + exp(term - acc))
          }
        }
      }

      nll <- nll - wt[i] * acc
    }

    nll
  }
}
