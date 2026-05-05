# =====================================================================
# Simulation Study --- replicates Section "Simulation Study" of
# theory/ml-glm/glm-theory.tex (Battaglia, Christensen, Hansen and
# Sacher, 2025).
#
# Notation matches the paper:
#   Y_i      scalar response
#   Z_i      latent K-category regressor (Z_i = 0,...,K-1)
#   Z_hat_i  observed proxy
#   x_i      cleanly observed covariates (intercept + q_i)
#   psi      = (gamma, alpha) -- gamma on latent dummies, alpha on x
#   p01,p10  binary misclassification rates  (Scenarios I & II)
#   pi       prevalence P(Z = 1)             (Scenarios I & II)
#   Pi, pi_z K x K misclassification matrix and class prevalences
#                                            (Scenario III, K = 4)
#
# Misclassification probabilities are estimated from a validation
# subsample, as in the paper. Each scenario reports per-estimator
# bias, empirical SD, and MSE for every parameter.
#
# Estimators compared:
#   - naive, BCA, BCM, CS, onestep        (mcglm)
#   - MC-SIMEX standard (Kuechenhoff)     (simex)
#   - MC-SIMEX improved (Sevilimedu-Yu)   (simex)
# =====================================================================

suppressPackageStartupMessages({
  library(mismeasured)
})

# ---- Configuration -------------------------------------------------------
# Defaults below give a ~5-10 minute demo on a laptop. To exactly replicate
# the paper's tables (B = 100, n = 10000 / 20000, full simex internal B),
# run as:
#
#   MCGLM_SIM_B=100 MCGLM_SIM_N_BIN=10000 MCGLM_SIM_N_BIN_VALID=1000 \
#     MCGLM_SIM_N_MULT=20000 Rscript simulations/simulation_study.R
#
# Expect ~4 hours at paper scale because simex (standard, internal B = 200)
# fits 200 * length(lambda) GLMs per outer rep on n = 10000-20000.
B            <- as.integer(Sys.getenv("MCGLM_SIM_B",            unset = "30"))
N_BIN        <- as.integer(Sys.getenv("MCGLM_SIM_N_BIN",        unset = "3000"))
N_BIN_VALID  <- as.integer(Sys.getenv("MCGLM_SIM_N_BIN_VALID",  unset = "300"))
N_MULT       <- as.integer(Sys.getenv("MCGLM_SIM_N_MULT",       unset = "5000"))
N_MULT_VALID_FRAC <- 0.10
USE_SIMEX_STD <- as.logical(Sys.getenv("MCGLM_SIM_SIMEX_STD",     unset = "TRUE"))
USE_SIMEX_IMP <- as.logical(Sys.getenv("MCGLM_SIM_SIMEX_IMP",     unset = "TRUE"))
USE_ONESTEP   <- as.logical(Sys.getenv("MCGLM_SIM_ONESTEP",       unset = "TRUE")) &&
                 requireNamespace("RTMB", quietly = TRUE)

OUT_DIR <- "simulations"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Validation-sample estimators of misclassification --------------------

#' Estimate (p01, p10) from a labelled validation sample (binary case).
#' Returns the empirical row-conditional flip rates with a small floor to
#' keep them strictly inside (0, 1).
.estimate_p01_p10 <- function(z_true, z_hat, eps = 1e-3) {
  n0 <- sum(z_true == 0); n1 <- sum(z_true == 1)
  p01 <- if (n0 > 0) sum(z_true == 0 & z_hat == 1) / n0 else 0
  p10 <- if (n1 > 0) sum(z_true == 1 & z_hat == 0) / n1 else 0
  list(p01 = max(min(p01, 1 - eps), eps),
       p10 = max(min(p10, 1 - eps), eps))
}

#' Estimate the K x K misclassification matrix from a validation sample.
#' Pi_hat[j, l] = P_hat(Z_hat = j-1 | Z = l-1).
.estimate_Pi <- function(z_true, z_hat, K, eps = 1e-3) {
  Pi <- matrix(0, K, K)
  for (l in seq_len(K)) {
    sel <- z_true == (l - 1L)
    if (any(sel)) {
      tab <- tabulate(z_hat[sel] + 1L, nbins = K)
      Pi[, l] <- tab / sum(sel)
    } else {
      Pi[, l] <- 1 / K
    }
  }
  Pi <- pmax(Pi, eps)
  Pi <- sweep(Pi, 2, colSums(Pi), "/")
  Pi
}

# ---- Coefficient extraction shims ----------------------------------------

#' Pull (gamma, alpha) from a simex fit and align to (gamma, intercept, q1)
#' ordering (matching mcglm's ("gamma", "alpha0", "alpha1")).
#' For K > 2 the simex fit reports K-1 dummies named "z1, z2, ...".
.extract_simex_bin <- function(fit) {
  cf <- coef(fit)
  c(gamma = unname(cf["1"]),                     # mc-dummy column
    alpha0 = unname(cf["(Intercept)"]),
    alpha1 = unname(cf["q1"]))
}
.extract_simex_multi <- function(fit, K) {
  cf <- coef(fit)
  gammas <- vapply(seq_len(K - 1L),
                   function(k) unname(cf[as.character(k)]),
                   numeric(1L))
  c(setNames(gammas, paste0("gamma", seq_len(K - 1L))),
    alpha0 = unname(cf["(Intercept)"]),
    alpha1 = unname(cf["q1"]))
}

# ---- Method runners ------------------------------------------------------

run_mcglm <- function(y, z_hat, x, family, p01, p10, pi_z,
                      methods = c("naive", "bca", "bcm", "cs"),
                      Pi_hat = NULL, K = 2L) {
  if (K == 2L) {
    fit <- mcglm(y, z_hat = z_hat, x = x, family = family,
                 method  = methods,
                 p01 = p01, p10 = p10, pi_z = pi_z,
                 fix_omega = TRUE)
  } else {
    fit <- mcglm(y, z_hat = z_hat, x = x, family = family,
                 method  = methods,
                 Pi = Pi_hat, pi_z = pi_z,
                 fix_omega = TRUE)
  }
  fit$coefficients
}

run_simex <- function(y, z_hat, x, family, Pi_hat, simex_method = "standard") {
  df <- data.frame(y = y, z = factor(z_hat, levels = seq_len(nrow(Pi_hat)) - 1L),
                   q1 = x[, 2])
  fit <- simex(y ~ mc(z, Pi_hat) + q1, family = family, data = df,
               method = simex_method, jackknife = FALSE)
  fit
}

# ---- Bias / SD / MSE summary --------------------------------------------

summarise_estimates <- function(est, psi0, scenario, n) {
  # est: B x methods x params  array
  rows <- list()
  for (m in dimnames(est)[[2]]) {
    for (k in seq_along(psi0)) {
      v   <- est[, m, k]
      v   <- v[is.finite(v)]
      tru <- psi0[k]
      bias <- mean(v) - tru
      sdv  <- sd(v)
      mse  <- mean((v - tru)^2)
      rows[[length(rows) + 1L]] <- data.frame(
        scenario  = scenario,
        n         = n,
        estimator = m,
        parameter = names(psi0)[k],
        truth     = unname(tru),
        bias      = bias,
        sd_emp    = sdv,
        mse       = mse,
        n_succ    = length(v),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

# ==========================================================================
# Scenario I --- Poisson regression with binary misclassification
# n = 10000, pi = 0.4, p01 = 0.10, p10 = 0.15, q ~ N(0, 1)
# Truth: gamma = 0.8, alpha0 = -0.5, alpha1 = 0.7
# Validation subsample of size 1000 to estimate (p01, p10).
# ==========================================================================

run_scenario_I <- function(reps = B, n = N_BIN, n_val = N_BIN_VALID,
                            seed = 20260505) {
  B <- reps
  message("\n== Scenario I: Poisson, binary misclassification ==")
  psi0 <- c(gamma = 0.8, alpha0 = -0.5, alpha1 = 0.7)
  pi_true <- 0.4
  p01_tr  <- 0.10
  p10_tr  <- 0.15
  family  <- poisson()

  methods_mcglm <- c("naive", "bca", "bcm", "cs")
  if (USE_ONESTEP) methods_mcglm <- c(methods_mcglm, "onestep")
  estimators <- c(methods_mcglm,
                  if (USE_SIMEX_STD) "simex_std",
                  if (USE_SIMEX_IMP) "simex_imp")

  est <- array(NA_real_, dim = c(B, length(estimators), 3),
               dimnames = list(NULL, estimators, names(psi0)))

  set.seed(seed); t0 <- Sys.time()
  for (b in seq_len(B)) {
    q     <- rnorm(n)
    x_mat <- cbind(1, q)
    z     <- rbinom(n, 1, pi_true)
    eta   <- psi0[1] * z + psi0[2] + psi0[3] * q
    y     <- rpois(n, exp(eta))
    z_hat <- z
    z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01_tr)
    z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10_tr)

    val_idx <- sample.int(n, n_val)
    p_hat   <- .estimate_p01_p10(z[val_idx], z_hat[val_idx])
    Pi_hat  <- matrix(c(1 - p_hat$p01, p_hat$p01,
                        p_hat$p10,     1 - p_hat$p10),
                      nrow = 2, byrow = FALSE)
    pi_hat_z <- mean(z[val_idx])  # estimate pi from validation data

    cf <- tryCatch(
      run_mcglm(y, z_hat, x_mat, family,
                p01 = p_hat$p01, p10 = p_hat$p10, pi_z = pi_hat_z,
                methods = methods_mcglm),
      error = function(e) NULL
    )
    if (!is.null(cf)) for (m in methods_mcglm) est[b, m, ] <- cf[[m]]

    if (USE_SIMEX_STD) {
      f1 <- tryCatch(run_simex(y, z_hat, x_mat, family, Pi_hat, "standard"),
                     error = function(e) NULL)
      if (!is.null(f1)) est[b, "simex_std", ] <- .extract_simex_bin(f1)
    }
    if (USE_SIMEX_IMP) {
      f2 <- tryCatch(run_simex(y, z_hat, x_mat, family, Pi_hat, "improved"),
                     error = function(e) NULL)
      if (!is.null(f2)) est[b, "simex_imp", ] <- .extract_simex_bin(f2)
    }
    if (b %% 10 == 0)
      message(sprintf("  rep %3d / %d  (%.1fs)",
                      b, B, as.numeric(Sys.time() - t0, units = "secs")))
  }
  summarise_estimates(est, psi0, "I_poisson_binary", n)
}

# ==========================================================================
# Scenario II --- Logistic regression with binary misclassification
# Same DGP as Scenario I but with binary response.
# ==========================================================================

run_scenario_II <- function(reps = B, n = N_BIN, n_val = N_BIN_VALID,
                             seed = 20260606) {
  B <- reps
  message("\n== Scenario II: Logistic, binary misclassification ==")
  psi0 <- c(gamma = 0.8, alpha0 = -0.5, alpha1 = 0.7)
  pi_true <- 0.4
  p01_tr  <- 0.10
  p10_tr  <- 0.15
  family  <- binomial()

  methods_mcglm <- c("naive", "bca", "bcm", "cs")
  if (USE_ONESTEP) methods_mcglm <- c(methods_mcglm, "onestep")
  estimators <- c(methods_mcglm,
                  if (USE_SIMEX_STD) "simex_std",
                  if (USE_SIMEX_IMP) "simex_imp")

  est <- array(NA_real_, dim = c(B, length(estimators), 3),
               dimnames = list(NULL, estimators, names(psi0)))

  set.seed(seed); t0 <- Sys.time()
  for (b in seq_len(B)) {
    q     <- rnorm(n)
    x_mat <- cbind(1, q)
    z     <- rbinom(n, 1, pi_true)
    eta   <- psi0[1] * z + psi0[2] + psi0[3] * q
    y     <- rbinom(n, 1, plogis(eta))
    z_hat <- z
    z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01_tr)
    z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10_tr)

    val_idx <- sample.int(n, n_val)
    p_hat   <- .estimate_p01_p10(z[val_idx], z_hat[val_idx])
    Pi_hat  <- matrix(c(1 - p_hat$p01, p_hat$p01,
                        p_hat$p10,     1 - p_hat$p10),
                      nrow = 2, byrow = FALSE)
    pi_hat_z <- mean(z[val_idx])

    cf <- tryCatch(
      run_mcglm(y, z_hat, x_mat, family,
                p01 = p_hat$p01, p10 = p_hat$p10, pi_z = pi_hat_z,
                methods = methods_mcglm),
      error = function(e) NULL
    )
    if (!is.null(cf)) for (m in methods_mcglm) est[b, m, ] <- cf[[m]]

    if (USE_SIMEX_STD) {
      f1 <- tryCatch(run_simex(y, z_hat, x_mat, family, Pi_hat, "standard"),
                     error = function(e) NULL)
      if (!is.null(f1)) est[b, "simex_std", ] <- .extract_simex_bin(f1)
    }
    if (USE_SIMEX_IMP) {
      f2 <- tryCatch(run_simex(y, z_hat, x_mat, family, Pi_hat, "improved"),
                     error = function(e) NULL)
      if (!is.null(f2)) est[b, "simex_imp", ] <- .extract_simex_bin(f2)
    }
    if (b %% 10 == 0)
      message(sprintf("  rep %3d / %d  (%.1fs)",
                      b, B, as.numeric(Sys.time() - t0, units = "secs")))
  }
  summarise_estimates(est, psi0, "II_logistic_binary", n)
}

# ==========================================================================
# Scenario III --- Poisson regression with K = 4 categorical
# misclassification.
#  pi_z   = (0.5, 0.25, 0.15, 0.10)
#  Pi     = 0.8 * I_4 + 0.2/3 * (J_4 - I_4)   (off-diagonal flip mass split equally)
# Truth: gamma = (1.0, -0.9, 0.2), alpha = (0.8, -0.7)
# Validation subsample = 10% of n.
# ==========================================================================

run_scenario_III <- function(reps = B, n = N_MULT, val_frac = N_MULT_VALID_FRAC,
                              seed = 20260707) {
  B <- reps
  message("\n== Scenario III: Poisson, K = 4 multicategory misclassification ==")
  K <- 4L
  pi_true <- c(0.50, 0.25, 0.15, 0.10)
  Pi_true <- 0.8 * diag(K) + (0.2 / (K - 1)) * (matrix(1, K, K) - diag(K))
  gamma_true <- c(0, 1.0, -0.9, 0.2)
  alpha_true <- c(0.8, -0.7)
  psi0 <- c(gamma1 = gamma_true[2], gamma2 = gamma_true[3], gamma3 = gamma_true[4],
            alpha0 = alpha_true[1], alpha1 = alpha_true[2])
  family <- poisson()

  methods_mcglm <- c("naive", "bca", "bcm", "cs")
  if (USE_ONESTEP) methods_mcglm <- c(methods_mcglm, "onestep")
  # Improved MC-SIMEX supports K-level mc, but only for a single mc() term.
  estimators <- c(methods_mcglm,
                  if (USE_SIMEX_STD) "simex_std",
                  if (USE_SIMEX_IMP) "simex_imp")

  est <- array(NA_real_, dim = c(B, length(estimators), length(psi0)),
               dimnames = list(NULL, estimators, names(psi0)))

  n_val <- ceiling(val_frac * n)
  set.seed(seed); t0 <- Sys.time()
  for (b in seq_len(B)) {
    q     <- rnorm(n)
    x_mat <- cbind(1, q)
    z     <- sample(0:(K - 1L), n, replace = TRUE, prob = pi_true)
    z_hat <- vapply(z, function(zi) sample(0:(K - 1L), 1, prob = Pi_true[, zi + 1L]),
                    integer(1L))
    eta   <- gamma_true[z + 1L] + alpha_true[1] + alpha_true[2] * q
    y     <- rpois(n, exp(eta))

    val_idx  <- sample.int(n, n_val)
    Pi_hat   <- .estimate_Pi(z[val_idx], z_hat[val_idx], K)
    pi_hat_z <- tabulate(z[val_idx] + 1L, nbins = K) / n_val

    cf <- tryCatch(
      run_mcglm(y, z_hat, x_mat, family,
                p01 = NULL, p10 = NULL, pi_z = pi_hat_z,
                methods = methods_mcglm,
                Pi_hat = Pi_hat, K = K),
      error = function(e) NULL
    )
    if (!is.null(cf)) for (m in methods_mcglm) est[b, m, ] <- cf[[m]]

    if (USE_SIMEX_STD) {
      f1 <- tryCatch(run_simex(y, z_hat, x_mat, family, Pi_hat, "standard"),
                     error = function(e) NULL)
      if (!is.null(f1)) est[b, "simex_std", ] <- .extract_simex_multi(f1, K)
    }
    if (USE_SIMEX_IMP) {
      f2 <- tryCatch(run_simex(y, z_hat, x_mat, family, Pi_hat, "improved"),
                     error = function(e) NULL)
      if (!is.null(f2)) est[b, "simex_imp", ] <- .extract_simex_multi(f2, K)
    }
    if (b %% 10 == 0)
      message(sprintf("  rep %3d / %d  (%.1fs)",
                      b, B, as.numeric(Sys.time() - t0, units = "secs")))
  }
  summarise_estimates(est, psi0, "III_poisson_multi_K4", n)
}

# ==========================================================================
# Driver
# ==========================================================================

main <- function() {
  message(sprintf(
    "Simulation study (B = %d, n_bin = %d, n_mult = %d, onestep = %s, simex_std = %s, simex_imp = %s)",
    B, N_BIN, N_MULT, USE_ONESTEP, USE_SIMEX_STD, USE_SIMEX_IMP
  ))
  res <- rbind(
    run_scenario_I(),
    run_scenario_II(),
    run_scenario_III()
  )
  out_csv <- file.path(OUT_DIR, "simulation_study.csv")
  write.csv(res, out_csv, row.names = FALSE)
  message(sprintf("Results written to %s", out_csv))

  # Pretty per-scenario summary
  for (sc in unique(res$scenario)) {
    cat(sprintf("\n=== %s (n = %d) ===\n",
                sc, unique(res$n[res$scenario == sc])))
    sub <- res[res$scenario == sc, c("estimator", "parameter", "truth",
                                     "bias", "sd_emp", "mse", "n_succ")]
    print(format(sub, digits = 4), row.names = FALSE)
  }
  invisible(res)
}

if (sys.nframe() == 0L) {
  main()
}
