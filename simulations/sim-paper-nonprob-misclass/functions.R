# =====================================================================
# Mass-imputation estimator under DRIFTING binary misclassification:
# simulation functions.
#
# Setting (theory/ml-glm/paper.tex, Section "Corrected-score mass
# imputation estimator", Theorem thm:mass-imp-cs):
#
#   * finite population U of size N with xi_i = (theta_i, 1, x_i),
#     Y_i ~ GLM(mu(psi0' xi_i)), canonical link;
#   * probability sample S_P (Poisson sampling, known pi_i, w_i = 1/pi_i)
#     carrying the TRUE xi_i but not Y_i;
#   * non-probability sample S_NP carrying Y_i and the PROXY
#     theta_hat_i, misclassified with DRIFTING rates
#         p01 = kappa01 / sqrt(n_NP),   p10 = kappa10 / sqrt(n_NP);
#   * mass-imputation mean estimator
#         Ybar_hat = N^{-1} sum_{i in S_P} w_i mu(psi_hat' xi_i),
#     with psi_hat = naive / BCA / BCM / CS from mcglm().
#
# NOTE: the package does NOT provide a variance estimator for the
# mass-imputation estimator.  The plug-in variance V_hat = V_P + V_NP
# (+ V_V and the phi-s covariance term for CS under estimated
# misclassification probabilities) is implemented HERE, directly from
# the formulas in the paper, using mismeasured's internal helpers for
# phi_i, I_hat and M_hat.
# =====================================================================

suppressPackageStartupMessages({
  library(mismeasured)
})

# ---- link functions (canonical links only, as required by the theory) ----
.link_funs <- function(family) {
  mismeasured:::.mcglm_get_link_funs(family)
}

# ---- population -----------------------------------------------------

#' Generate a finite population.
#'
#' @param N population size.
#' @param family "poisson" or "binomial".
#' @param psi0 c(gamma, alpha0, alpha1).
#' @param pi_theta P(theta = 1).
#' @param dep_theta_x if TRUE, theta depends on x (violates the
#'   theta-independent-of-q assumption; robustness scenario).
make_population <- function(N, family, psi0, pi_theta,
                            dep_theta_x = FALSE) {
  x1 <- rnorm(N)
  if (dep_theta_x) {
    theta <- rbinom(N, 1, plogis(qlogis(pi_theta) + 0.8 * x1))
  } else {
    theta <- rbinom(N, 1, pi_theta)
  }
  eta <- psi0[1] * theta + psi0[2] + psi0[3] * x1
  y <- switch(family,
    poisson  = rpois(N, exp(eta)),
    binomial = rbinom(N, 1, plogis(eta)),
    stop("unsupported family: ", family)
  )
  list(x1 = x1, theta = theta, y = y, N = N,
       ybar_true = mean(y))
}

# ---- samples --------------------------------------------------------

#' Informative non-probability sample: selection depends on x but not
#' on Y (transportability holds).  The intercept is calibrated so that
#' the expected sample size equals n_np.
#'
#' IMPORTANT: by default selection must NOT depend on theta.  Selecting
#' on (x + theta) jointly makes theta and x dependent WITHIN S_NP
#' (collider stratification) and shifts P(theta = 1 | S_NP) away from
#' pi_theta, violating the theta-independent-of-q assumption underlying
#' the c1/c2 drift correction.  b_theta != 0 is reserved for the
#' assumption-violation arm.
draw_np_sample <- function(pop, n_np, b_x = 0.5, b_theta = 0) {
  lin <- b_x * pop$x1 + b_theta * pop$theta
  a0 <- uniroot(function(a) mean(plogis(a + lin)) - n_np / pop$N,
                c(-30, 10))$root
  sel <- rbinom(pop$N, 1, plogis(a0 + lin)) == 1L
  which(sel)
}

#' Probability sample: Poisson sampling with size measure exp(0.3 x).
#' Poisson sampling makes pi_ij = pi_i pi_j, so the HT variance
#' estimator from the paper reduces to a single sum (see mi_variance).
draw_p_sample <- function(pop, n_p) {
  size <- exp(0.3 * pop$x1)
  pi_i <- pmin(n_p * size / sum(size), 1)
  sel  <- rbinom(pop$N, 1, pi_i) == 1L
  list(idx = which(sel), pi = pi_i[sel])
}

#' Nondifferential binary misclassification.
misclassify <- function(theta, p01, p10) {
  theta_hat <- theta
  i0 <- theta == 0
  i1 <- theta == 1
  theta_hat[i0] <- rbinom(sum(i0), 1, p01)
  theta_hat[i1] <- 1 - rbinom(sum(i1), 1, p10)
  theta_hat
}

#' Standard estimators of (p01, p10, pi) from a validation sample in
#' which both theta and theta_hat are observed.
estimate_eta <- function(theta_v, theta_hat_v) {
  n0 <- sum(theta_v == 0)
  n1 <- sum(theta_v == 1)
  list(
    p01 = if (n0 > 0) sum(theta_v == 0 & theta_hat_v == 1) / n0 else 0,
    p10 = if (n1 > 0) sum(theta_v == 1 & theta_hat_v == 0) / n1 else 0,
    pi  = mean(theta_v)
  )
}

# ---- mass-imputation point estimator --------------------------------

#' Mass-imputation mean, true xi in S_P.
#'
#' type = "ht":    Ybar = N^{-1}      sum_{S_P} w_i mu(psi' xi_i)
#'                 (the form in Theorem thm:mass-imp-cs / Kim et al. 2021)
#' type = "hajek": Ybar = N_hat^{-1}  sum_{S_P} w_i mu(psi' xi_i),
#'                 N_hat = sum_{S_P} w_i  (ratio-stabilized; removes the
#'                 random-sample-size variance of Poisson sampling).
mi_point <- function(psi, xi_p, w, N, mu_fun, type = c("ht", "hajek")) {
  type  <- match.arg(type)
  denom <- if (type == "ht") N else sum(w)
  sum(w * mu_fun(as.numeric(xi_p %*% psi))) / denom
}

# ---- mass-imputation variance estimator (NOT in the package) --------
#
# Implements, from the paper:
#   V_P  = N^{-2} sum_{S_P} (1 - pi_i) / pi_i^2 * mu_i^2
#          (HT variance of the paper specialized to Poisson sampling,
#           where pi_ij = pi_i pi_j kills all off-diagonal terms; the
#           UNCENTERED mu_i^2 is correct for the 1/N HT form because it
#           includes the random-sample-size component)
#   V_NP = n_NP^{-1} G_P' A^{-1} S (A^{-1})' G_P
#          with G_P = N^{-1} sum_{S_P} w_i mudot_i xi_i.
#
# For type = "hajek" the design term is the standard ratio
# linearization (influence w_i (mu_i - Ybar) / N_hat):
#   V_P  = N_hat^{-2} sum_{S_P} (1 - pi_i) / pi_i^2 * (mu_i - Ybar)^2
# and G_P uses N_hat^{-1} in place of N^{-1}.  The psi-uncertainty
# structure (V_NP, V_V, V_cov) is unchanged.  NOTE: the paper states
# the linearization only for the 1/N form; the Hajek variant is the
# standard ratio extension and should be flagged as such in the text.
#
# Per Theorems thm:mis-bc / thm:mass-imp-cs:
#   * naive / bca / bcm: A = I_hat, S = mean(score score') -- the naive
#     sandwich, which is the asymptotic variance of BCA/BCM in the
#     drifting regime;
#   * cs: A = I_hat + M_hat, S = mean(phi phi') with the m-corrected phi.
#
# For CS under estimated (p01, p10) (subsection "Corrected-score under
# estimated misclassification probabilities") we additionally return
#   V_V   = n_V^{-1}  G' H^{-1} D V(s) D' (H^{-1})' G
#   V_cov = 2 n_NP^{-1} G' H^{-1} Cov(phi_i, s_i) D' (H^{-1})' G
# (the covariance term applies because S_V is a subsample of S_NP),
# with D = -d mbar / d eta', eta = (p01, p10)':
#   d m_gamma / d p01 = -(1-pi) delta,  d m_gamma / d p10 = 0
#   d m_alpha / d p01 = -(1-pi) delta q, d m_alpha / d p10 = +pi delta q.

mi_variance <- function(psi, method, y_np, z_hat_np, x_np, family,
                        c1, c2, p01, p10, pi_theta,
                        xi_p, w, pi_p, N,
                        val = NULL, type = c("ht", "hajek")) {
  type  <- match.arg(type)
  fam   <- .link_funs(family)
  n_np  <- length(y_np)
  xi_np <- cbind(z_hat_np, x_np)
  p     <- ncol(xi_np)

  # -- V_P: Poisson-sampling design variance ---------------------------
  mu_p <- fam$mu(as.numeric(xi_p %*% psi))
  if (type == "ht") {
    denom <- N
    V_P   <- sum((1 - pi_p) / pi_p^2 * mu_p^2) / N^2
  } else {
    denom <- sum(w)
    ybar  <- sum(w * mu_p) / denom
    V_P   <- sum((1 - pi_p) / pi_p^2 * (mu_p - ybar)^2) / denom^2
  }

  # -- G_P -------------------------------------------------------------
  mudot_p <- fam$mu_dot(as.numeric(xi_p %*% psi))
  G_P <- colSums(w * mudot_p * xi_p) / denom

  # -- A, S on the non-probability sample ------------------------------
  resid <- y_np - fam$mu(as.numeric(xi_np %*% psi))
  score <- xi_np * resid
  I_hat <- mismeasured:::.mcglm_compute_Ihat(psi, xi_np, fam$mu_dot)

  use_cs_form <- method == "cs"
  if (use_cs_form) {
    m_mat <- mismeasured:::.mcglm_compute_m_bin(psi, x_np, fam$mu, c1, c2)
    M_hat <- mismeasured:::.mcglm_compute_Mhat_bin(psi, x_np, fam$mu,
                                                   fam$mu_dot, c1, c2)
    A   <- I_hat + M_hat
    phi <- score - m_mat
  } else {
    A   <- I_hat
    phi <- score
  }
  A_inv <- solve(A)
  S     <- crossprod(phi) / n_np
  hAg   <- t(A_inv) %*% G_P                       # (A^{-1})' G_P
  V_NP  <- as.numeric(t(hAg) %*% S %*% hAg) / n_np

  V_V <- 0; V_cov <- 0
  if (use_cs_form && !is.null(val)) {
    # D = -d mbar / d eta' (p x 2), see header
    delta <- mismeasured:::.mcglm_compute_delta(psi, x_np, fam$mu)
    D <- matrix(0, p, 2)
    D[1, 1]  <- (1 - pi_theta) * mean(delta)
    D[-1, 1] <- (1 - pi_theta) * colMeans(delta * x_np)
    D[-1, 2] <- -pi_theta * colMeans(delta * x_np)

    # influence s_i of (p01_hat, p10_hat), i in S_V
    s_mat <- cbind(
      (val$theta == 0) * ((val$theta_hat == 1) - p01) / (1 - pi_theta),
      (val$theta == 1) * ((val$theta_hat == 0) - p10) / pi_theta
    )
    V_s  <- stats::cov(s_mat)
    DVsD <- D %*% V_s %*% t(D)
    V_V  <- as.numeric(t(hAg) %*% DVsD %*% hAg) / nrow(s_mat)

    # covariance term: S_V is a subsample of S_NP, so phi_i and s_i are
    # observed jointly for i in S_V
    phi_v  <- phi[val$np_pos, , drop = FALSE]
    C_ps   <- stats::cov(phi_v, s_mat)            # p x 2
    V_cov  <- 2 * as.numeric(t(hAg) %*% C_ps %*% t(D) %*% hAg) / n_np
  }

  list(V_P = V_P, V_NP = V_NP, V_V = V_V, V_cov = V_cov,
       V_total = V_P + V_NP + V_V + V_cov)
}

# ---- one Monte Carlo replicate --------------------------------------

#' @param cfg list with: family, psi0, pi_theta, kappa01, kappa10,
#'   n_np, ratio_p (n_p / n_np), eta ("known" | "estimated"),
#'   frac_v (validation fraction of S_NP), N_mult (N / n_np),
#'   dep_theta_x.
#' @return data.frame, one row per method, or NULL on failure.
run_one_rep <- function(cfg, rep_id) {
  fam <- .link_funs(cfg$family)

  N   <- cfg$N_mult * cfg$n_np
  pop <- make_population(N, cfg$family, cfg$psi0, cfg$pi_theta,
                         dep_theta_x = isTRUE(cfg$dep_theta_x))

  # drifting misclassification rates
  p01 <- cfg$kappa01 / sqrt(cfg$n_np)
  p10 <- cfg$kappa10 / sqrt(cfg$n_np)
  stopifnot(p01 < 1, p10 < 1)

  # samples
  idx_np <- draw_np_sample(pop, cfg$n_np,
                           b_theta = if (is.null(cfg$sel_theta)) 0
                                     else cfg$sel_theta)
  sp     <- draw_p_sample(pop, round(cfg$ratio_p * cfg$n_np))
  n_np   <- length(idx_np)

  y_np      <- pop$y[idx_np]
  x_np      <- cbind(1, pop$x1[idx_np])
  theta_np  <- pop$theta[idx_np]
  z_hat_np  <- misclassify(theta_np, p01, p10)

  xi_p <- cbind(pop$theta[sp$idx], 1, pop$x1[sp$idx])
  w    <- 1 / sp$pi

  # misclassification parameters: known truth or validation estimates
  val <- NULL
  if (cfg$eta == "known") {
    p01_use <- p01; p10_use <- p10; pi_use <- cfg$pi_theta
  } else {
    n_v    <- max(100L, round(cfg$frac_v * n_np))
    np_pos <- sample.int(n_np, min(n_v, n_np))
    val    <- list(theta     = theta_np[np_pos],
                   theta_hat = z_hat_np[np_pos],
                   np_pos    = np_pos)
    eta_hat <- estimate_eta(val$theta, val$theta_hat)
    p01_use <- eta_hat$p01; p10_use <- eta_hat$p10; pi_use <- eta_hat$pi
    if (pi_use <= 0 || pi_use >= 1) return(NULL)
  }
  c1 <- p01_use * (1 - pi_use)
  c2 <- p01_use * (1 - pi_use) - p10_use * pi_use

  fit <- tryCatch(
    mcglm(y_np, z_hat = z_hat_np, x = x_np, family = cfg$family,
          method = c("naive", "bca", "bcm", "cs"),
          p01 = p01_use, p10 = p10_use, pi_z = pi_use),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  # oracle: true theta in S_NP, no correction needed
  psi_list <- fit$coefficients
  psi_list$oracle <- tryCatch(
    unname(coef(glm(y_np ~ 0 + theta_np + x_np,
                    family = cfg$family))),
    error = function(e) NULL
  )

  combos <- expand.grid(m = names(psi_list), est = c("ht", "hajek"),
                        stringsAsFactors = FALSE)
  rows <- lapply(seq_len(nrow(combos)), function(k) {
    m   <- combos$m[k]
    est <- combos$est[k]
    psi <- psi_list[[m]]
    if (is.null(psi) || anyNA(psi)) return(NULL)
    ybar <- mi_point(psi, xi_p, w, N, fam$mu, type = est)
    if (m == "oracle") {
      vc <- mi_variance(psi, "naive", y_np, theta_np, x_np, cfg$family,
                        0, 0, 0, 0, cfg$pi_theta,
                        xi_p, w, sp$pi, N, type = est)
    } else {
      vc <- mi_variance(psi, m, y_np, z_hat_np, x_np, cfg$family,
                        c1, c2, p01_use, p10_use, pi_use,
                        xi_p, w, sp$pi, N, val = val, type = est)
    }
    data.frame(
      rep = rep_id, method = m, estimator = est,
      ybar_hat = ybar, ybar_true = pop$ybar_true,
      psi_gamma = psi[1], psi_alpha0 = psi[2], psi_alpha1 = psi[3],
      V_P = vc$V_P, V_NP = vc$V_NP, V_V = vc$V_V, V_cov = vc$V_cov,
      se_hat = sqrt(max(vc$V_total, 0)),
      p01_used = p01_use, p10_used = p10_use,
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)

  for (nm in c("family", "kappa01", "kappa10", "n_np", "ratio_p", "eta"))
    out[[nm]] <- cfg[[nm]]
  out
}
