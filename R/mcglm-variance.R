# ---------------------------------------------------------------------------
# Sandwich variance estimators for mcglm methods
# ---------------------------------------------------------------------------

#' Sandwich variance for the naive estimator
#' @keywords internal
.mcglm_vcov_naive <- function(psi, y, xi_hat, family, wt = NULL) {
  fam    <- .normalize_family(family)
  n      <- length(y)
  N      <- if (is.null(wt)) n else sum(wt)
  eta    <- as.numeric(xi_hat %*% psi)
  mu_val <- fam$linkinv(eta)
  w      <- fam$mu.eta(eta)
  eps    <- y - mu_val

  if (is.null(wt)) {
    A <- crossprod(xi_hat * w, xi_hat) / N
    C <- crossprod(xi_hat * eps, xi_hat * eps) / N
  } else {
    A <- crossprod(xi_hat * (wt * w), xi_hat) / N
    C <- crossprod(xi_hat * (wt * eps), xi_hat * eps) / N
  }
  A_inv <- solve(A)
  A_inv %*% C %*% A_inv / N
}

#' Sandwich variance for BCA/BCM estimators (binary)
#' @keywords internal
.mcglm_vcov_bc_bin <- function(psi_bc, y, xi_hat, x, family,
                               p01 = NULL, p10 = NULL, pi_z = NULL,
                               c1 = NULL, c2 = NULL,
                               psi_naive = NULL,
                               type = c("bca", "bcm"), corrected = FALSE,
                               wt = NULL) {
  if (!corrected) return(.mcglm_vcov_naive(psi_bc, y, xi_hat, family, wt = wt))

  type   <- match.arg(type)
  fam    <- .normalize_family(family)
  n      <- length(y)
  N      <- if (is.null(wt)) n else sum(wt)
  p      <- length(psi_bc)

  if (is.null(psi_naive)) psi_naive <- psi_bc
  eta    <- as.numeric(xi_hat %*% psi_naive)
  w      <- fam$mu.eta(eta)
  resid  <- y - fam$linkinv(eta)

  if (is.null(wt)) {
    A_hat <- crossprod(xi_hat * w, xi_hat) / N
  } else {
    A_hat <- crossprod(xi_hat * (wt * w), xi_hat) / N
  }
  A_inv <- solve(A_hat)
  c1_loc <- if (!is.null(c1)) c1 else p01 * (1 - pi_z)
  c2_loc <- if (!is.null(c2)) c2 else p01 * (1 - pi_z) - p10 * pi_z
  if (is.null(c1_loc) || is.null(c2_loc) ||
      length(c1_loc) == 0L || length(c2_loc) == 0L)
    stop("Corrected variance for 'bca'/'bcm' requires either ",
         "(p01, p10, pi_z) or (c1, c2).")
  M_hat <- .mcglm_compute_Mhat_bin(psi_naive, x, fam$linkinv, fam$mu.eta, c1_loc,
                                    c2_loc, wt = wt)
  m_mat <- .mcglm_compute_m_bin(psi_naive, x, fam$linkinv, c1_loc, c2_loc)
  if (is.null(wt)) {
    m_bar <- colMeans(m_mat)
  } else {
    m_bar <- colSums(wt * m_mat) / N
  }

  if (type == "bca") {
    G     <- diag(p) - A_inv %*% M_hat
    H_inv <- A_inv
  } else {
    H     <- A_hat + M_hat
    H_inv <- solve(H)
    G     <- diag(p) - H_inv %*% M_hat
  }

  score_mat  <- xi_hat * resid
  centered_m <- sweep(m_mat, 2, m_bar)
  L1 <- A_inv %*% t(G)
  L2 <- t(H_inv)

  if (is.null(wt)) {
    S_ss <- crossprod(score_mat) / N^2
    S_mm <- crossprod(centered_m) / N^2
    S_sm <- crossprod(score_mat, centered_m) / N^2
  } else {
    S_ss <- crossprod(score_mat * wt, score_mat) / N^2
    S_mm <- crossprod(centered_m * wt, centered_m) / N^2
    S_sm <- crossprod(score_mat * wt, centered_m) / N^2
  }

  t(L1) %*% S_ss %*% L1 + t(L2) %*% S_mm %*% L2 -
    t(L1) %*% S_sm %*% L2 - t(L2) %*% t(S_sm) %*% L1
}

#' Sandwich variance for corrected-score estimator (binary)
#' @keywords internal
.mcglm_vcov_cs_bin <- function(psi, y, xi_hat, x, family, p01, p10, pi_z,
                               c1 = NULL, c2 = NULL,
                               wt = NULL, validation = NULL) {
  fam    <- .normalize_family(family)
  n      <- length(y)
  N      <- if (is.null(wt)) n else sum(wt)

  eta_tilde <- as.numeric(xi_hat %*% psi)
  resid     <- y - fam$linkinv(eta_tilde)
  c1_loc <- if (!is.null(c1)) c1 else p01 * (1 - pi_z)
  c2_loc <- if (!is.null(c2)) c2 else p01 * (1 - pi_z) - p10 * pi_z
  if (is.null(c1_loc) || is.null(c2_loc) ||
      length(c1_loc) == 0L || length(c2_loc) == 0L)
    stop("Variance for 'cs' requires either (p01, p10, pi_z) or (c1, c2).")
  m_mat     <- .mcglm_compute_m_bin(psi, x, fam$linkinv, c1_loc, c2_loc)

  phi_mat <- xi_hat * resid - m_mat

  if (is.null(wt)) {
    S <- crossprod(phi_mat) / N
  } else {
    S <- crossprod(phi_mat * wt, phi_mat) / N
  }

  if (!is.null(validation))
    S <- .mcglm_cs_validation_meat(S, phi_mat, psi, x, 2L,
                                    fam$linkinv, validation)

  I_hat <- .mcglm_compute_Ihat(psi, xi_hat, fam$mu.eta, wt = wt)
  M_hat <- .mcglm_compute_Mhat_bin(psi, x, fam$linkinv, fam$mu.eta, c1_loc, c2_loc,
                                    wt = wt)
  J     <- -(I_hat + M_hat)
  J_inv <- solve(J)

  J_inv %*% S %*% t(J_inv) / N
}


# Same fall-through for the BCA/BCM corrected sandwich.
# (The earlier definition above already exists; here we extend the binary
#  helper signature to accept c1/c2 the same way.)


# ---- Multicategory variance estimators ----

#' Sandwich variance for naive estimator (multicategory)
#' @keywords internal
.mcglm_vcov_naive_multi <- function(psi, y, xi_hat, z_hat, x, K, family,
                                    wt = NULL) {
  fam <- .normalize_family(family)
  n   <- length(y)
  N   <- if (is.null(wt)) n else sum(wt)
  s   <- K - 1
  r   <- ncol(x)

  gamma <- c(0, psi[seq_len(s)])
  alpha <- psi[(s + 1):(s + r)]
  eta_base  <- as.numeric(x %*% alpha)
  eta_tilde <- eta_base + gamma[z_hat + 1]

  w   <- fam$mu.eta(eta_tilde)
  eps <- y - fam$linkinv(eta_tilde)

  if (is.null(wt)) {
    A <- crossprod(xi_hat * w, xi_hat) / N
    C <- crossprod(xi_hat * eps, xi_hat * eps) / N
  } else {
    A <- crossprod(xi_hat * (wt * w), xi_hat) / N
    C <- crossprod(xi_hat * (wt * eps), xi_hat * eps) / N
  }
  A_inv <- solve(A)
  A_inv %*% C %*% A_inv / N
}

#' Sandwich variance for BCA/BCM estimators (multicategory)
#' @keywords internal
.mcglm_vcov_bc_multi <- function(psi_bc, y, xi_hat, z_hat, x, K, family,
                                 Pi = NULL, pi_z = NULL,
                                 psi_naive = NULL,
                                 type = c("bca", "bcm"), corrected = FALSE,
                                 wt = NULL,
                                 jacobian = c("analytical", "numerical")) {
  jacobian <- match.arg(jacobian)
  if (!corrected)
    return(.mcglm_vcov_naive_multi(psi_bc, y, xi_hat, z_hat, x, K, family,
                                    wt = wt))

  type <- match.arg(type)
  fam  <- .normalize_family(family)
  n    <- length(y)
  N    <- if (is.null(wt)) n else sum(wt)
  s    <- K - 1
  r    <- ncol(x)
  p    <- s + r

  if (is.null(psi_naive)) psi_naive <- psi_bc

  gamma <- c(0, psi_naive[seq_len(s)])
  alpha <- psi_naive[(s + 1):(s + r)]
  eta_base  <- as.numeric(x %*% alpha)
  eta_tilde <- eta_base + gamma[z_hat + 1]

  w     <- fam$mu.eta(eta_tilde)
  resid <- y - fam$linkinv(eta_tilde)

  if (is.null(wt)) {
    A_hat <- crossprod(xi_hat * w, xi_hat) / N
  } else {
    A_hat <- crossprod(xi_hat * (wt * w), xi_hat) / N
  }
  A_inv <- solve(A_hat)
  M_hat <- .mcglm_compute_Mhat_multi(psi_naive, x, K, fam$linkinv, Pi, pi_z,
                                      wt = wt, jacobian = jacobian,
                                      mu_dot_fun = fam$mu.eta)
  m_mat <- .mcglm_compute_m_multi(psi_naive, x, K, fam$linkinv, Pi, pi_z)
  if (is.null(wt)) {
    m_bar <- colMeans(m_mat)
  } else {
    m_bar <- colSums(wt * m_mat) / N
  }

  if (type == "bca") {
    G     <- diag(p) - A_inv %*% M_hat
    H_inv <- A_inv
  } else {
    H     <- A_hat + M_hat
    H_inv <- solve(H)
    G     <- diag(p) - H_inv %*% M_hat
  }

  score_mat  <- xi_hat * resid
  centered_m <- sweep(m_mat, 2, m_bar)
  L1 <- A_inv %*% t(G)
  L2 <- t(H_inv)

  if (is.null(wt)) {
    S_ss <- crossprod(score_mat) / N^2
    S_mm <- crossprod(centered_m) / N^2
    S_sm <- crossprod(score_mat, centered_m) / N^2
  } else {
    S_ss <- crossprod(score_mat * wt, score_mat) / N^2
    S_mm <- crossprod(centered_m * wt, centered_m) / N^2
    S_sm <- crossprod(score_mat * wt, centered_m) / N^2
  }

  t(L1) %*% S_ss %*% L1 + t(L2) %*% S_mm %*% L2 -
    t(L1) %*% S_sm %*% L2 - t(L2) %*% t(S_sm) %*% L1
}

#' Sandwich variance for corrected-score estimator (multicategory)
#' @keywords internal
.mcglm_vcov_cs_multi <- function(psi, y, xi_hat, z_hat, x, K, family, Pi, pi_z,
                                 wt = NULL,
                                 jacobian = c("analytical", "numerical"),
                                 validation = NULL) {
  jacobian <- match.arg(jacobian)
  fam <- .normalize_family(family)
  n   <- length(y)
  N   <- if (is.null(wt)) n else sum(wt)
  s   <- K - 1
  r   <- ncol(x)

  gamma <- c(0, psi[seq_len(s)])
  alpha <- psi[(s + 1):(s + r)]
  eta_base  <- as.numeric(x %*% alpha)
  eta_tilde <- eta_base + gamma[z_hat + 1]

  resid <- y - fam$linkinv(eta_tilde)
  m_mat <- .mcglm_compute_m_multi(psi, x, K, fam$linkinv, Pi, pi_z)

  phi_mat <- xi_hat * resid - m_mat
  if (is.null(wt)) {
    S <- crossprod(phi_mat) / N
  } else {
    S <- crossprod(phi_mat * wt, phi_mat) / N
  }

  if (!is.null(validation))
    S <- .mcglm_cs_validation_meat(S, phi_mat, psi, x, K,
                                    fam$linkinv, validation)

  I_hat <- .mcglm_compute_Ihat_multi(psi, xi_hat, z_hat, K, fam$mu.eta,
                                      wt = wt)
  M_hat <- .mcglm_compute_Mhat_multi(psi, x, K, fam$linkinv, Pi, pi_z, wt = wt,
                                      jacobian = jacobian,
                                      mu_dot_fun = fam$mu.eta)
  J     <- -(I_hat + M_hat)
  J_inv <- solve(J)

  J_inv %*% S %*% t(J_inv) / N
}

# Validate the sampling information before the variance tryCatch in mcglm.
# Joint cells B[j,l] = P(proxy=j, true=l) are sufficient for both K=2 and
# K>2. For K=2, B[2,1]=c1 and B[1,2]=c1-c2; this is exactly the binary
# (pi,p01,p10) influence-function formula by the chain rule.
.mcglm_prepare_cs_validation <- function(validation, z_hat, K, wt,
                                          c1, c2, Pi, pi_z) {
  if (!is.null(wt) && any(wt != 1))
    stop("validation currently requires unweighted observations or unit weights.",
         call. = FALSE)
  vd <- .mcglm_parse_validation(validation, z_hat, K)
  z <- vd$z
  nv <- vd$n
  index <- vd$index
  proxy <- vd$proxy
  cells <- proxy + 1L + K * z  # column-major vec(B), as in the paper
  b <- tabulate(cells, nbins = K * K) / nv
  B <- matrix(b, K, K)
  if (any(colSums(B) == 0))
    stop("Every true category must occur in validation.", call. = FALSE)
  expected <- if (K == 2L) c(B[2, 1], B[2, 1] - B[1, 2]) else b
  supplied <- if (K == 2L) c(c1, c2) else as.numeric(sweep(Pi, 2, pi_z, "*"))
  if (length(supplied) != length(expected) || any(!is.finite(supplied)) ||
      max(abs(supplied - expected)) > 1e-8)
    stop("CS probabilities must be empirical estimates from the supplied validation sample; ",
         "supply pi_z explicitly (or c1/c2 for K=2).", call. = FALSE)
  list(b = b, cells = cells, index = index, n = nv,
       type = if (is.null(index)) "external" else "internal")
}

# Equations Sigmahat-cs-external/internal and D-multicategory-explicit in
# GLM bias correction.tex. All quantities here are on the sqrt(n) scale;
# the existing sandwich callers divide the resulting covariance by n.
.mcglm_cs_validation_meat <- function(S, phi_mat, psi, x, K, mu_fun,
                                       validation) {
  n <- nrow(x)
  nv <- validation$n
  s <- K - 1L
  gamma <- c(0, psi[seq_len(s)])
  eta <- as.numeric(x %*% psi[-seq_len(s)])
  mu <- vapply(gamma, function(g) mu_fun(eta + g), numeric(n))
  D <- matrix(0, length(psi), K * K)
  for (ell in seq_len(K)) {
    for (j in seq_len(K)) {
      column <- j + K * (ell - 1L)
      delta <- mu[, j] - mu[, ell]
      if (j > 1L) D[j - 1L, column] <- mean(delta)
      D[-seq_len(s), column] <- colMeans(x * delta)
    }
  }

  # Row i equals D s_i. Computing in coefficient space avoids allocating
  # an n_V by K^2 matrix and never inverts the singular cell covariance.
  influence <- sweep(t(D[, validation$cells, drop = FALSE]), 2,
                     as.numeric(D %*% validation$b))
  S <- S + (n / nv) * crossprod(influence) / nv
  if (!is.null(validation$index)) {
    phi_v <- phi_mat[validation$index, , drop = FALSE]
    phi_v <- sweep(phi_v, 2, colMeans(phi_v))
    cross <- crossprod(phi_v, influence) / nv
    S <- S + cross + t(cross)
  }
  S
}
