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
# subsample, as in the paper. Per-rep results carry both estimates
# and asymptotic SEs; the summary reports bias, empirical SD, mean
# SE, MSE, and Wald 95% coverage.
#
# Estimators compared:
#   - naive, BCA, BCM, CS, onestep        (mcglm with the correct family)
#   - <method>_g                          (mcglm with family = gaussian())
#   - <method>_b                          (mcglm with family = binomial());
#                                         only fitted on Sim I and Sim III
#                                         after binarising y to 1{y > 0}
#                                         (extensive-margin / incidence
#                                         model). Skipped for Sim II
#                                         where binomial is the main fit.
#   - MC-SIMEX standard (Kuechenhoff)     (simex)
#   - MC-SIMEX improved (Sevilimedu-Yu)   (simex)
#   - ols: lm(y ~ z_hat + q) on the raw response (count y for Poisson,
#     0/1 y for logistic = linear probability model). The OLS coefficient
#     on z is *not* the GLM gamma -- it is a linear-projection slope --
#     so its "bias" against psi_0 conflates link-scale mismatch with the
#     classical errors-in-variables attenuation.  naive_g should equal
#     ols (modulo numerical tolerance) and serves as a sanity check.
#
# Parallelism: reps are dispatched via future.apply::future_lapply.
# Set MCGLM_SIM_WORKERS to override the worker count (default:
# parallel::detectCores() - 1). Set to 1 for serial.
# =====================================================================

suppressPackageStartupMessages({
  library(mismeasured)
  library(future)
  library(future.apply)
  library(progressr)
})

# ---- Configuration -------------------------------------------------------
# Defaults below give a ~5-10 minute demo on a laptop. To exactly replicate
# the paper's tables (B = 100, n = 10000 / 20000, full simex internal B),
# run as:
#
#   MCGLM_SIM_B=100 MCGLM_SIM_N_BIN=10000 MCGLM_SIM_N_BIN_VALID=1000 \
#     MCGLM_SIM_N_MULT=20000 Rscript simulations/simulation_study.R
#
B            <- as.integer(Sys.getenv("MCGLM_SIM_B",            unset = "100"))
N_BIN        <- as.integer(Sys.getenv("MCGLM_SIM_N_BIN",        unset = "10000"))
N_BIN_VALID  <- as.integer(Sys.getenv("MCGLM_SIM_N_BIN_VALID",  unset = "1000"))
N_MULT       <- as.integer(Sys.getenv("MCGLM_SIM_N_MULT",       unset = "20000"))
N_MULT_VALID_FRAC <- 0.10
USE_SIMEX_STD <- as.logical(Sys.getenv("MCGLM_SIM_SIMEX_STD",     unset = "TRUE"))
USE_SIMEX_IMP <- as.logical(Sys.getenv("MCGLM_SIM_SIMEX_IMP",     unset = "TRUE"))
USE_ONESTEP   <- as.logical(Sys.getenv("MCGLM_SIM_ONESTEP",       unset = "TRUE")) &&
                 requireNamespace("RTMB", quietly = TRUE)

WORKERS <- as.integer(Sys.getenv("MCGLM_SIM_WORKERS",
                                 unset = max(1L, parallel::detectCores() - 1L)))

OUT_DIR <- "simulations"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Future plan ---------------------------------------------------------
if (WORKERS > 1L) {
  plan(multisession, workers = WORKERS)
} else {
  plan(sequential)
}
on.exit(plan(sequential), add = TRUE)

# ---- Validation-sample estimators of misclassification --------------------

#' Estimate (p01, p10) from a labelled validation sample (binary case).
.estimate_p01_p10 <- function(z_true, z_hat, eps = 1e-3) {
  n0 <- sum(z_true == 0); n1 <- sum(z_true == 1)
  p01 <- if (n0 > 0) sum(z_true == 0 & z_hat == 1) / n0 else 0
  p10 <- if (n1 > 0) sum(z_true == 1 & z_hat == 0) / n1 else 0
  list(p01 = max(min(p01, 1 - eps), eps),
       p10 = max(min(p10, 1 - eps), eps))
}

#' Estimate the K x K misclassification matrix from a validation sample.
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
  sweep(Pi, 2, colSums(Pi), "/")
}

# ---- Coefficient extraction shims ----------------------------------------

#' Pull (est, se) from a simex fit and align to mcglm's parameter ordering.
.extract_simex_bin <- function(fit) {
  cf <- coef(fit)
  V  <- fit$vcov
  se <- if (!is.null(V)) sqrt(pmax(diag(V), 0)) else
        setNames(rep(NA_real_, length(cf)), names(cf))
  list(
    est = c(gamma  = unname(cf["1"]),
            alpha0 = unname(cf["(Intercept)"]),
            alpha1 = unname(cf["q1"])),
    se  = c(gamma  = unname(se["1"]),
            alpha0 = unname(se["(Intercept)"]),
            alpha1 = unname(se["q1"]))
  )
}

.extract_simex_multi <- function(fit, K) {
  cf <- coef(fit)
  V  <- fit$vcov
  se <- if (!is.null(V)) sqrt(pmax(diag(V), 0)) else
        setNames(rep(NA_real_, length(cf)), names(cf))
  gnames <- as.character(seq_len(K - 1L))
  list(
    est = c(setNames(cf[gnames], paste0("gamma", seq_len(K - 1L))),
            alpha0 = unname(cf["(Intercept)"]),
            alpha1 = unname(cf["q1"])),
    se  = c(setNames(se[gnames], paste0("gamma", seq_len(K - 1L))),
            alpha0 = unname(se["(Intercept)"]),
            alpha1 = unname(se["q1"]))
  )
}

# ---- Method runners ------------------------------------------------------

run_mcglm <- function(y, z_hat, x, family, p01, p10, pi_z,
                      methods = c("naive", "bca", "bcm", "cs"),
                      Pi_hat = NULL, K = 2L) {
  fit <- if (K == 2L) {
    mcglm(y, z_hat = z_hat, x = x, family = family, method = methods,
          p01 = p01, p10 = p10, pi_z = pi_z, fix_omega = TRUE)
  } else {
    mcglm(y, z_hat = z_hat, x = x, family = family, method = methods,
          Pi = Pi_hat, pi_z = pi_z, fix_omega = TRUE)
  }
  list(coef = fit$coefficients, se = fit$se)
}

run_simex <- function(y, z_hat, x, family, Pi_hat, simex_method = "standard") {
  df <- data.frame(y = y,
                   z = factor(z_hat, levels = seq_len(nrow(Pi_hat)) - 1L),
                   q1 = x[, 2])
  simex(y ~ mc(z, Pi_hat) + q1, family = family, data = df,
        method = simex_method, jackknife = TRUE)
}

#' OLS on the raw response: count regression run as `lm(y ~ ...)` and the
#' linear probability model run as `lm(y ~ ...)`. The reported coefficients
#' are aligned to the GLM parameter names so they can be compared to the
#' GLM truth ψ_0 even though OLS estimates a linear-projection slope, not
#' the link-scale γ. The gap between the two captures both the link-scale
#' mismatch *and* the classical errors-in-variables attenuation.
fit_ols_bin <- function(y, z_hat, q) {
  fit <- lm(y ~ z_hat + q)
  cf  <- coef(fit)
  se  <- sqrt(pmax(diag(vcov(fit)), 0))
  list(
    est = c(gamma  = unname(cf["z_hat"]),
            alpha0 = unname(cf["(Intercept)"]),
            alpha1 = unname(cf["q"])),
    se  = c(gamma  = unname(se["z_hat"]),
            alpha0 = unname(se["(Intercept)"]),
            alpha1 = unname(se["q"]))
  )
}

fit_ols_multi <- function(y, z_hat, q, K) {
  z_factor <- factor(z_hat, levels = 0:(K - 1L))
  fit <- lm(y ~ z_factor + q)
  cf  <- coef(fit)
  se  <- sqrt(pmax(diag(vcov(fit)), 0))
  gnames_in <- paste0("z_factor", seq_len(K - 1L))
  list(
    est = c(setNames(cf[gnames_in], paste0("gamma", seq_len(K - 1L))),
            alpha0 = unname(cf["(Intercept)"]),
            alpha1 = unname(cf["q"])),
    se  = c(setNames(se[gnames_in], paste0("gamma", seq_len(K - 1L))),
            alpha0 = unname(se["(Intercept)"]),
            alpha1 = unname(se["q"]))
  )
}

# ---- Tidy collection -----------------------------------------------------

#' Merge an auxiliary mcglm fit into the main fit, suffixing method names
#' with `suffix` (e.g. "_g" for gaussian, "_b" for binomial).
.add_suffixed_fit <- function(fits, fit_aux, suffix) {
  if (is.null(fit_aux)) return(fits)
  names(fit_aux$coef) <- paste0(names(fit_aux$coef), suffix)
  names(fit_aux$se)   <- paste0(names(fit_aux$se),   suffix)
  if (is.null(fits)) return(fit_aux)
  list(coef = c(fits$coef, fit_aux$coef),
       se   = c(fits$se,   fit_aux$se))
}

#' Build a long data.frame with one row per (estimator, parameter) for one rep.
.collect_long <- function(scenario, n, rep_id, mc_fit, simex_fits, param_names) {
  P <- length(param_names)
  parts <- list()
  if (!is.null(mc_fit)) {
    for (m in names(mc_fit$coef)) {
      cf <- mc_fit$coef[[m]]
      if (is.null(cf) || length(cf) != P) next
      sg <- mc_fit$se[[m]]
      if (is.null(sg) || length(sg) != P) sg <- rep(NA_real_, P)
      parts[[length(parts) + 1L]] <- data.frame(
        scenario = scenario, n = n, rep = rep_id,
        estimator = m, parameter = param_names,
        estimate = unname(cf), se = unname(sg),
        stringsAsFactors = FALSE
      )
    }
  }
  for (m in names(simex_fits)) {
    ex <- simex_fits[[m]]
    if (is.null(ex)) next
    parts[[length(parts) + 1L]] <- data.frame(
      scenario = scenario, n = n, rep = rep_id,
      estimator = m, parameter = param_names,
      estimate = unname(ex$est), se = unname(ex$se),
      stringsAsFactors = FALSE
    )
  }
  if (length(parts) == 0L) return(NULL)
  do.call(rbind, parts)
}

# ---- Per-rep simulators --------------------------------------------------
# Each returns a long data.frame for ONE rep, or NULL on total failure.

sim_one_I <- function(rep_id, n = N_BIN, n_val = N_BIN_VALID) {
  psi0 <- c(gamma = 0.8, alpha0 = -0.5, alpha1 = 0.7)
  family <- poisson()
  pi_true <- 0.4; p01_tr <- 0.10; p10_tr <- 0.15

  q     <- rnorm(n); x_mat <- cbind(1, q)
  z     <- rbinom(n, 1, pi_true)
  eta   <- psi0[1] * z + psi0[2] + psi0[3] * q
  y     <- rpois(n, exp(eta))
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01_tr)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10_tr)

  val_idx <- sample.int(n, n_val)
  p_hat   <- .estimate_p01_p10(z[val_idx], z_hat[val_idx])
  Pi_hat  <- matrix(c(1 - p_hat$p01, p_hat$p01,
                      p_hat$p10,     1 - p_hat$p10), nrow = 2)
  pi_hat_z <- mean(z[val_idx])

  methods_mcglm <- c("naive", "bca", "bcm", "cs",
                     if (USE_ONESTEP) "onestep")
  mc_fit <- tryCatch(
    run_mcglm(y, z_hat, x_mat, family,
              p01 = p_hat$p01, p10 = p_hat$p10, pi_z = pi_hat_z,
              methods = methods_mcglm),
    error = function(e) NULL
  )
  # mcglm with the wrong family: gaussian (linear-in-y) and binomial on
  # the binarised count (extensive-margin / incidence model).
  mc_fit_g <- tryCatch(
    run_mcglm(y, z_hat, x_mat, gaussian(),
              p01 = p_hat$p01, p10 = p_hat$p10, pi_z = pi_hat_z,
              methods = methods_mcglm),
    error = function(e) NULL
  )
  y_bin <- as.integer(y > 0)
  mc_fit_b <- tryCatch(
    run_mcglm(y_bin, z_hat, x_mat, binomial(),
              p01 = p_hat$p01, p10 = p_hat$p10, pi_z = pi_hat_z,
              methods = methods_mcglm),
    error = function(e) NULL
  )
  mc_fit <- .add_suffixed_fit(mc_fit,   mc_fit_g, "_g")
  mc_fit <- .add_suffixed_fit(mc_fit,   mc_fit_b, "_b")

  simex_fits <- list()
  if (USE_SIMEX_STD) simex_fits$simex_std <-
    tryCatch(.extract_simex_bin(run_simex(y, z_hat, x_mat, family, Pi_hat, "standard")),
             error = function(e) NULL)
  if (USE_SIMEX_IMP) simex_fits$simex_imp <-
    tryCatch(.extract_simex_bin(run_simex(y, z_hat, x_mat, family, Pi_hat, "improved")),
             error = function(e) NULL)
  simex_fits$ols <- tryCatch(fit_ols_bin(y, z_hat, q), error = function(e) NULL)

  .collect_long("I_poisson_binary", n, rep_id, mc_fit, simex_fits, names(psi0))
}

sim_one_II <- function(rep_id, n = N_BIN, n_val = N_BIN_VALID) {
  psi0 <- c(gamma = 0.8, alpha0 = -0.5, alpha1 = 0.7)
  family <- binomial()
  pi_true <- 0.4; p01_tr <- 0.10; p10_tr <- 0.15

  q     <- rnorm(n); x_mat <- cbind(1, q)
  z     <- rbinom(n, 1, pi_true)
  eta   <- psi0[1] * z + psi0[2] + psi0[3] * q
  y     <- rbinom(n, 1, plogis(eta))
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01_tr)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10_tr)

  val_idx <- sample.int(n, n_val)
  p_hat   <- .estimate_p01_p10(z[val_idx], z_hat[val_idx])
  Pi_hat  <- matrix(c(1 - p_hat$p01, p_hat$p01,
                      p_hat$p10,     1 - p_hat$p10), nrow = 2)
  pi_hat_z <- mean(z[val_idx])

  methods_mcglm <- c("naive", "bca", "bcm", "cs",
                     if (USE_ONESTEP) "onestep")
  mc_fit <- tryCatch(
    run_mcglm(y, z_hat, x_mat, family,
              p01 = p_hat$p01, p10 = p_hat$p10, pi_z = pi_hat_z,
              methods = methods_mcglm),
    error = function(e) NULL
  )
  # Wrong-family comparator: gaussian on 0/1 y (linear probability model
  # with bias correction). No "_b" version: binomial *is* the main family.
  mc_fit_g <- tryCatch(
    run_mcglm(y, z_hat, x_mat, gaussian(),
              p01 = p_hat$p01, p10 = p_hat$p10, pi_z = pi_hat_z,
              methods = methods_mcglm),
    error = function(e) NULL
  )
  mc_fit <- .add_suffixed_fit(mc_fit, mc_fit_g, "_g")

  simex_fits <- list()
  if (USE_SIMEX_STD) simex_fits$simex_std <-
    tryCatch(.extract_simex_bin(run_simex(y, z_hat, x_mat, family, Pi_hat, "standard")),
             error = function(e) NULL)
  if (USE_SIMEX_IMP) simex_fits$simex_imp <-
    tryCatch(.extract_simex_bin(run_simex(y, z_hat, x_mat, family, Pi_hat, "improved")),
             error = function(e) NULL)
  simex_fits$ols <- tryCatch(fit_ols_bin(y, z_hat, q), error = function(e) NULL)

  .collect_long("II_logistic_binary", n, rep_id, mc_fit, simex_fits, names(psi0))
}

sim_one_III <- function(rep_id, n = N_MULT, val_frac = N_MULT_VALID_FRAC) {
  K <- 4L
  pi_true <- c(0.50, 0.25, 0.15, 0.10)
  Pi_true <- 0.8 * diag(K) + (0.2 / (K - 1)) * (matrix(1, K, K) - diag(K))
  gamma_true <- c(0, 1.0, -0.9, 0.2); alpha_true <- c(0.8, -0.7)
  psi0 <- c(gamma1 = gamma_true[2], gamma2 = gamma_true[3], gamma3 = gamma_true[4],
            alpha0 = alpha_true[1], alpha1 = alpha_true[2])
  family <- poisson()

  n_val <- ceiling(val_frac * n)
  q     <- rnorm(n); x_mat <- cbind(1, q)
  z     <- sample(0:(K - 1L), n, replace = TRUE, prob = pi_true)
  z_hat <- vapply(z, function(zi) sample(0:(K - 1L), 1, prob = Pi_true[, zi + 1L]),
                  integer(1L))
  eta   <- gamma_true[z + 1L] + alpha_true[1] + alpha_true[2] * q
  y     <- rpois(n, exp(eta))

  val_idx  <- sample.int(n, n_val)
  Pi_hat   <- .estimate_Pi(z[val_idx], z_hat[val_idx], K)
  pi_hat_z <- tabulate(z[val_idx] + 1L, nbins = K) / n_val

  methods_mcglm <- c("naive", "bca", "bcm", "cs",
                     if (USE_ONESTEP) "onestep")
  mc_fit <- tryCatch(
    run_mcglm(y, z_hat, x_mat, family,
              p01 = NULL, p10 = NULL, pi_z = pi_hat_z,
              methods = methods_mcglm, Pi_hat = Pi_hat, K = K),
    error = function(e) NULL
  )
  mc_fit_g <- tryCatch(
    run_mcglm(y, z_hat, x_mat, gaussian(),
              p01 = NULL, p10 = NULL, pi_z = pi_hat_z,
              methods = methods_mcglm, Pi_hat = Pi_hat, K = K),
    error = function(e) NULL
  )
  y_bin <- as.integer(y > 0)
  mc_fit_b <- tryCatch(
    run_mcglm(y_bin, z_hat, x_mat, binomial(),
              p01 = NULL, p10 = NULL, pi_z = pi_hat_z,
              methods = methods_mcglm, Pi_hat = Pi_hat, K = K),
    error = function(e) NULL
  )
  mc_fit <- .add_suffixed_fit(mc_fit, mc_fit_g, "_g")
  mc_fit <- .add_suffixed_fit(mc_fit, mc_fit_b, "_b")

  simex_fits <- list()
  if (USE_SIMEX_STD) simex_fits$simex_std <-
    tryCatch(.extract_simex_multi(run_simex(y, z_hat, x_mat, family, Pi_hat, "standard"), K),
             error = function(e) NULL)
  if (USE_SIMEX_IMP) simex_fits$simex_imp <-
    tryCatch(.extract_simex_multi(run_simex(y, z_hat, x_mat, family, Pi_hat, "improved"), K),
             error = function(e) NULL)
  simex_fits$ols <- tryCatch(fit_ols_multi(y, z_hat, q, K), error = function(e) NULL)

  .collect_long("III_poisson_multi_K4", n, rep_id, mc_fit, simex_fits, names(psi0))
}

# Drifting-regime variant: Pi shrinks at rate n^{-1/2}.  Disabled by
# default; uncomment the run_scenario("IIIb", ...) call in main() to enable.
sim_one_IIIb <- function(rep_id, n, c_off = 1.0,
                         val_frac = N_MULT_VALID_FRAC) {
  K <- 4L
  pi_true <- c(0.50, 0.25, 0.15, 0.10)
  off_per_col <- c_off / sqrt(n)
  Pi_true <- (1 - off_per_col) * diag(K) +
             (off_per_col / (K - 1L)) * (matrix(1, K, K) - diag(K))
  gamma_true <- c(0, 1.0, -0.9, 0.2); alpha_true <- c(0.8, -0.7)
  psi0 <- c(gamma1 = gamma_true[2], gamma2 = gamma_true[3], gamma3 = gamma_true[4],
            alpha0 = alpha_true[1], alpha1 = alpha_true[2])
  family <- poisson()

  n_val <- ceiling(val_frac * n)
  q     <- rnorm(n); x_mat <- cbind(1, q)
  z     <- sample(0:(K - 1L), n, replace = TRUE, prob = pi_true)
  z_hat <- vapply(z, function(zi) sample(0:(K - 1L), 1, prob = Pi_true[, zi + 1L]),
                  integer(1L))
  eta   <- gamma_true[z + 1L] + alpha_true[1] + alpha_true[2] * q
  y     <- rpois(n, exp(eta))

  val_idx  <- sample.int(n, n_val)
  Pi_hat   <- .estimate_Pi(z[val_idx], z_hat[val_idx], K)
  pi_hat_z <- tabulate(z[val_idx] + 1L, nbins = K) / n_val

  methods_mcglm <- c("naive", "bca", "bcm", "cs",
                     if (USE_ONESTEP) "onestep")
  mc_fit <- tryCatch(
    run_mcglm(y, z_hat, x_mat, family,
              p01 = NULL, p10 = NULL, pi_z = pi_hat_z,
              methods = methods_mcglm, Pi_hat = Pi_hat, K = K),
    error = function(e) NULL
  )
  mc_fit_g <- tryCatch(
    run_mcglm(y, z_hat, x_mat, gaussian(),
              p01 = NULL, p10 = NULL, pi_z = pi_hat_z,
              methods = methods_mcglm, Pi_hat = Pi_hat, K = K),
    error = function(e) NULL
  )
  y_bin <- as.integer(y > 0)
  mc_fit_b <- tryCatch(
    run_mcglm(y_bin, z_hat, x_mat, binomial(),
              p01 = NULL, p10 = NULL, pi_z = pi_hat_z,
              methods = methods_mcglm, Pi_hat = Pi_hat, K = K),
    error = function(e) NULL
  )
  mc_fit <- .add_suffixed_fit(mc_fit, mc_fit_g, "_g")
  mc_fit <- .add_suffixed_fit(mc_fit, mc_fit_b, "_b")

  simex_fits <- list()
  if (USE_SIMEX_STD) simex_fits$simex_std <-
    tryCatch(.extract_simex_multi(run_simex(y, z_hat, x_mat, family, Pi_hat, "standard"), K),
             error = function(e) NULL)
  if (USE_SIMEX_IMP) simex_fits$simex_imp <-
    tryCatch(.extract_simex_multi(run_simex(y, z_hat, x_mat, family, Pi_hat, "improved"), K),
             error = function(e) NULL)
  simex_fits$ols <- tryCatch(fit_ols_multi(y, z_hat, q, K), error = function(e) NULL)

  scenario <- sprintf("IIIb_drifting_c%.2f", c_off)
  .collect_long(scenario, n, rep_id, mc_fit, simex_fits, names(psi0))
}

# ---- Parallel scenario dispatcher ---------------------------------------

#' Run B reps of `sim_fn` in parallel via future_lapply.
#' Returns a long data.frame with all reps stacked.
run_reps <- function(sim_fn, B, ..., seed) {
  message(sprintf("  dispatching %d reps across %s",
                  B, if (WORKERS > 1L) sprintf("%d workers", WORKERS) else "1 process"))
  t0 <- Sys.time()
  with_progress({
    p <- progressor(steps = B)
    out <- future_lapply(seq_len(B), function(b) {
      res <- sim_fn(b, ...)
      p()
      res
    }, future.seed = seed, future.globals = TRUE)
  })
  out <- Filter(Negate(is.null), out)
  message(sprintf("    done in %.1fs (%d / %d reps succeeded)",
                  as.numeric(Sys.time() - t0, units = "secs"), length(out), B))
  do.call(rbind, out)
}

# ---- Truth tables --------------------------------------------------------

.truth_table <- function() {
  rbind(
    data.frame(scenario = "I_poisson_binary",
               parameter = c("gamma", "alpha0", "alpha1"),
               truth = c(0.8, -0.5, 0.7)),
    data.frame(scenario = "II_logistic_binary",
               parameter = c("gamma", "alpha0", "alpha1"),
               truth = c(0.8, -0.5, 0.7)),
    data.frame(scenario = "III_poisson_multi_K4",
               parameter = c("gamma1", "gamma2", "gamma3", "alpha0", "alpha1"),
               truth = c(1.0, -0.9, 0.2, 0.8, -0.7))
  )
}

# Drifting variants share the K=4 truth (gamma1, gamma2, gamma3, alpha0, alpha1).
.truth_table_IIIb <- function(c_off = 1.0) {
  data.frame(scenario  = sprintf("IIIb_drifting_c%.2f", c_off),
             parameter = c("gamma1", "gamma2", "gamma3", "alpha0", "alpha1"),
             truth     = c(1.0, -0.9, 0.2, 0.8, -0.7))
}

# ---- Summariser ----------------------------------------------------------

#' Aggregate a long raw data.frame to per-(scenario, n, estimator, parameter)
#' bias / sd / mse / coverage. Wald intervals at the requested confidence level.
summarise_long <- function(raw, truth_df, level = 0.95) {
  z_crit <- qnorm(0.5 + level / 2)
  d <- merge(raw, truth_df, by = c("scenario", "parameter"), sort = FALSE)
  d$err <- d$estimate - d$truth
  d$in_ci <- with(d,
    ifelse(is.finite(estimate) & is.finite(se) & se > 0,
           (estimate - z_crit * se <= truth) & (truth <= estimate + z_crit * se),
           NA))

  key  <- with(d, interaction(scenario, n, estimator, parameter, drop = TRUE))
  rows <- lapply(split(d, key), function(g) {
    v   <- g$estimate[is.finite(g$estimate)]
    s   <- g$se[is.finite(g$se) & g$se > 0]
    tru <- g$truth[1]
    data.frame(
      scenario  = g$scenario[1],
      n         = g$n[1],
      estimator = g$estimator[1],
      parameter = g$parameter[1],
      truth     = tru,
      bias      = if (length(v)) mean(v) - tru else NA_real_,
      sd_emp    = if (length(v) > 1L) sd(v) else NA_real_,
      mean_se   = if (length(s)) mean(s) else NA_real_,
      mse       = if (length(v)) mean((v - tru)^2) else NA_real_,
      coverage  = mean(g$in_ci, na.rm = TRUE),
      n_succ    = length(v),
      n_se      = length(s),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL
  out[order(out$scenario, out$n, out$estimator, out$parameter), ]
}

# ---- Driver --------------------------------------------------------------

main <- function() {
  message(sprintf(
    "Simulation study (B = %d, n_bin = %d, n_mult = %d, workers = %d, onestep = %s, simex_std = %s, simex_imp = %s)",
    B, N_BIN, N_MULT, WORKERS, USE_ONESTEP, USE_SIMEX_STD, USE_SIMEX_IMP
  ))

  message("\n== Scenario I: Poisson, binary misclassification ==")
  raw_I   <- run_reps(sim_one_I,   B, seed = 20260505L)

  message("\n== Scenario II: Logistic, binary misclassification ==")
  raw_II  <- run_reps(sim_one_II,  B, seed = 20260606L)

  message("\n== Scenario III: Poisson, K = 4 multicategory misclassification ==")
  raw_III <- run_reps(sim_one_III, B, seed = 20260707L)

  raw <- rbind(raw_I, raw_II, raw_III)
  truth_df <- .truth_table()

  # ---- (Optional) drifting-regime sweep --------------------------------
  # Uncomment to add Scenario IIIb at multiple n values.
  # n_grid <- c(2500L, 5000L, 10000L, 20000L)
  # message("\n== Scenario IIIb: Poisson, K = 4, drifting regime ==")
  # raw_IIIb <- do.call(rbind, lapply(n_grid, function(n) {
  #   message(sprintf("  n = %d", n))
  #   run_reps(sim_one_IIIb, B, n = n, seed = 20260808L + n)
  # }))
  # raw      <- rbind(raw, raw_IIIb)
  # truth_df <- rbind(truth_df, .truth_table_IIIb())

  summary <- summarise_long(raw, truth_df)

  out_csv <- file.path(OUT_DIR, "simulation_study.csv")
  out_rds <- file.path(OUT_DIR, "simulation_study.rds")
  write.csv(summary, out_csv, row.names = FALSE)
  saveRDS(list(summary  = summary,
               raw      = raw,
               truth    = truth_df,
               settings = list(B = B, N_BIN = N_BIN, N_BIN_VALID = N_BIN_VALID,
                               N_MULT = N_MULT, N_MULT_VALID_FRAC = N_MULT_VALID_FRAC,
                               WORKERS = WORKERS,
                               USE_ONESTEP = USE_ONESTEP,
                               USE_SIMEX_STD = USE_SIMEX_STD,
                               USE_SIMEX_IMP = USE_SIMEX_IMP),
               session  = sessionInfo()),
          out_rds)
  message(sprintf("Summary written to %s\nRaw data + summary written to %s",
                  out_csv, out_rds))

  for (sc in unique(summary$scenario)) {
    cat(sprintf("\n=== %s ===\n", sc))
    sub <- summary[summary$scenario == sc,
                   c("n", "estimator", "parameter", "truth", "bias",
                     "sd_emp", "mean_se", "mse", "coverage", "n_succ", "n_se")]
    print(format(sub, digits = 4), row.names = FALSE)
  }
  invisible(list(summary = summary, raw = raw))
}

if (sys.nframe() == 0L) {
  main()
}
