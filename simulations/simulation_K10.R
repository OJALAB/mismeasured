# =====================================================================
# K = 10 misclassification simulation study.
#
# Design choice: this script targets the high-K regime motivated by
# applied settings such as occupation (ISCO-1d, K ~= 10) where the
# K x K misclassification matrix Pi has K(K-1) = 90 free off-diagonal
# entries, and validation cells for rare classes are thin.  Estimators:
#
#   - naive     : two-step proxy-score MLE                     (mcglm)
#   - bca       : additive bias-corrected estimator            (mcglm)
#   - bcm       : multiplicative bias-corrected estimator      (mcglm)
#   - cs        : corrected-score, drift formulation           (mcglm)
#   - cs_akn    : corrected-score, AKN98 unbiased surrogate     (mcglm)
#                 (Akazawa, Kinukawa & Nakamura 1998; needs only Pi,
#                 not pi_z)
#   - onestep   : joint mixture-likelihood, Battaglia et al.   (mcglm,
#                 requires RTMB)
#   - simex_std : MC-SIMEX, standard (Kuechenhoff et al. 2006)
#   - simex_imp : MC-SIMEX, improved (Sevilimedu & Yu 2026)
#
# Two scenarios are run, both with the same K = 10 misclassification
# structure and the same coefficient values but different outcome
# families:
#
#   K10_poisson  : log-link Poisson regression (count outcome)
#   K10_logistic : logit-link binomial regression (binary outcome)
#
# Misclassification probabilities and prevalences are estimated from a
# validation subsample exactly as in Scenario III of glm-theory.tex.
#
# Parallelism: future.apply::future_lapply over reps.  Override via env
# vars; defaults give a multi-minute run on a laptop.
#
#   MCGLM_SIM_K10_B          number of Monte Carlo replications
#   MCGLM_SIM_K10_N          primary sample size            (default 10000)
#   MCGLM_SIM_K10_N_VAL      validation subsample size      (default 1000)
#   MCGLM_SIM_K10_WORKERS    parallel workers (default detectCores() - 1)
#   MCGLM_SIM_K10_ONESTEP    "TRUE"/"FALSE", whether to fit onestep
#   MCGLM_SIM_K10_SIMEX_STD  "TRUE"/"FALSE", whether to fit standard simex
#   MCGLM_SIM_K10_SIMEX_IMP  "TRUE"/"FALSE", whether to fit improved simex
#   MCGLM_SIM_K10_JACKKNIFE  "TRUE"/"FALSE", whether simex computes
#                            jackknife SEs (slow; default FALSE at K = 10)
#
# Outputs: simulations/simulation_K10.csv (summary table)
#          simulations/simulation_K10.rds (raw + summary + settings)
# =====================================================================

suppressPackageStartupMessages({
  library(mismeasured)
  library(future)
  library(future.apply)
  library(progressr)
})

# ---- Configuration -------------------------------------------------------

B             <- as.integer(Sys.getenv("MCGLM_SIM_K10_B",         unset = "100"))
N             <- as.integer(Sys.getenv("MCGLM_SIM_K10_N",         unset = "10000"))
N_VAL         <- as.integer(Sys.getenv("MCGLM_SIM_K10_N_VAL",     unset = "1000"))
USE_ONESTEP   <- as.logical(Sys.getenv("MCGLM_SIM_K10_ONESTEP",   unset = "TRUE")) &&
                 requireNamespace("RTMB", quietly = TRUE)
USE_SIMEX_STD <- as.logical(Sys.getenv("MCGLM_SIM_K10_SIMEX_STD", unset = "TRUE"))
USE_SIMEX_IMP <- as.logical(Sys.getenv("MCGLM_SIM_K10_SIMEX_IMP", unset = "TRUE"))
SIMEX_JK      <- as.logical(Sys.getenv("MCGLM_SIM_K10_JACKKNIFE", unset = "FALSE"))

WORKERS <- as.integer(Sys.getenv("MCGLM_SIM_K10_WORKERS",
                                 unset = max(1L, parallel::detectCores() - 1L)))

OUT_DIR <- "simulations"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- DGP -----------------------------------------------------------------
# K = 10, skewed class prevalences typical of occupational data.
# Pi_true: 0.80 on the diagonal, the remaining 0.20 mass split uniformly
# across the K - 1 = 9 off-diagonal cells in each column.  This matches
# the multicategory pattern in Scenario III of glm-theory.tex (extended
# from K = 4 to K = 10) and keeps each off-diagonal entry around 0.022,
# i.e. roughly the regime where validation cell counts start to bite.

K_DGP        <- 10L
PI_TRUE_DIAG <- 0.80

# Skewed prevalences summing to 1 (roughly Zipf-ish but not too extreme).
PI_PREV_TRUE <- c(0.22, 0.16, 0.13, 0.11, 0.09, 0.08, 0.07, 0.06, 0.05, 0.03)
stopifnot(length(PI_PREV_TRUE) == K_DGP, abs(sum(PI_PREV_TRUE) - 1) < 1e-12)

# Latent-dummy coefficients (baseline gamma_0 = 0).  Mix of signs and
# magnitudes; calibrated so eta stays in a sensible range for q ~ N(0,1)
# both for log-link Poisson (exp(eta) bounded) and logit-link binomial
# (probabilities away from 0/1).  With alpha_0 = 0.5, alpha_1 = -0.6 and
# gamma in [-0.8, 0.9], eta ranges roughly in [-2.5, 2.5].
GAMMA_TRUE <- c(0,
                 0.90, -0.80,  0.60, -0.50,  0.40, -0.30,  0.20, -0.10,  0.05)
stopifnot(length(GAMMA_TRUE) == K_DGP)

ALPHA_TRUE <- c(intercept = 0.5, q1 = -0.6)

# Build truth vector matching the parameter ordering used by mcglm():
# (gamma_1, ..., gamma_{K-1}, alpha_0, alpha_1).
PSI0 <- c(setNames(GAMMA_TRUE[-1L], paste0("gamma", seq_len(K_DGP - 1L))),
          alpha0 = unname(ALPHA_TRUE[1]),
          alpha1 = unname(ALPHA_TRUE[2]))
PARAM_NAMES <- names(PSI0)

build_Pi_true <- function(K, diag_p) {
  off <- (1 - diag_p) / (K - 1L)
  diag_p * diag(K) + off * (matrix(1, K, K) - diag(K))
}
PI_TRUE <- build_Pi_true(K_DGP, PI_TRUE_DIAG)

# ---- Future plan ---------------------------------------------------------
if (WORKERS > 1L) {
  plan(multisession, workers = WORKERS)
} else {
  plan(sequential)
}
on.exit(plan(sequential), add = TRUE)

# ---- Validation-sample estimator of Pi -----------------------------------
# Adds floor `eps` to avoid exact zeros, then renormalises columns.  This
# mirrors what an applied user would do when an off-diagonal validation
# cell happens to be empty.
estimate_Pi <- function(z_true, z_hat, K, eps = 1e-3) {
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

# ---- SIMEX runner and extractor -----------------------------------------
# At K = 10 the misclassification matrix is 10 x 10 and each off-diagonal
# entry is roughly 0.022; this is the regime where MC-SIMEX
# extrapolation is noisiest, so SIMEX fits are wrapped in tryCatch and a
# fit failure on a given rep is silently dropped.

run_simex_K10 <- function(y, z_hat, x_mat, family_obj, Pi_hat,
                          simex_method, jackknife = FALSE) {
  df <- data.frame(y  = y,
                   z  = factor(z_hat, levels = seq_len(nrow(Pi_hat)) - 1L),
                   q1 = x_mat[, 2])
  simex(y ~ mc(z, Pi_hat) + q1, family = family_obj, data = df,
        method = simex_method, jackknife = jackknife)
}

#' Pull (est, se) from a multicategory simex fit and align to mcglm's
#' parameter ordering: (gamma_1, ..., gamma_{K-1}, alpha0, alpha1).
extract_simex_K10 <- function(fit, K) {
  cf <- coef(fit)
  V  <- fit$vcov
  se <- if (!is.null(V)) sqrt(pmax(diag(V), 0)) else
        setNames(rep(NA_real_, length(cf)), names(cf))
  gnames_in <- as.character(seq_len(K - 1L))
  list(
    est = c(setNames(unname(cf[gnames_in]), paste0("gamma", seq_len(K - 1L))),
            alpha0 = unname(cf["(Intercept)"]),
            alpha1 = unname(cf["q1"])),
    se  = c(setNames(unname(se[gnames_in]), paste0("gamma", seq_len(K - 1L))),
            alpha0 = unname(se["(Intercept)"]),
            alpha1 = unname(se["q1"]))
  )
}

# ---- Per-rep simulator ---------------------------------------------------
# `family_name` selects the outcome family: "poisson" => log-link counts,
# "logistic" => logit-link 0/1.  Everything else (latent Z, Z_hat,
# covariates, regression coefficients, validation subsample) is shared
# across families so only the outcome model differs.

sim_one <- function(rep_id, family_name,
                    n = N, n_val = N_VAL) {
  K <- K_DGP

  q     <- rnorm(n)
  x_mat <- cbind(1, q)
  z     <- sample(0:(K - 1L), n, replace = TRUE, prob = PI_PREV_TRUE)
  z_hat <- vapply(z,
                  function(zi) sample(0:(K - 1L), 1L,
                                      prob = PI_TRUE[, zi + 1L]),
                  integer(1L))
  eta   <- GAMMA_TRUE[z + 1L] + ALPHA_TRUE[1] + ALPHA_TRUE[2] * q

  family_obj <- switch(family_name,
    poisson  = poisson(),
    logistic = binomial(),
    stop("Unknown family_name: ", family_name)
  )
  y <- switch(family_name,
    poisson  = rpois(n, exp(eta)),
    logistic = rbinom(n, 1L, plogis(eta))
  )

  val_idx  <- sample.int(n, n_val)
  Pi_hat   <- estimate_Pi(z[val_idx], z_hat[val_idx], K)
  pi_hat_z <- tabulate(z[val_idx] + 1L, nbins = K) / n_val

  methods <- c("naive", "bca", "bcm", "cs", "cs_akn",
               if (USE_ONESTEP) "onestep")

  scenario <- sprintf("K10_%s_n%d", family_name, n)
  parts <- list()

  fit <- tryCatch(
    mcglm(y, z_hat = z_hat, x = x_mat, family = family_obj,
          method = methods, Pi = Pi_hat, pi_z = pi_hat_z,
          fix_omega = TRUE, jacobian = "analytical"),
    error = function(e) NULL
  )
  if (!is.null(fit)) {
    for (m in names(fit$coefficients)) {
      cf <- fit$coefficients[[m]]
      if (is.null(cf) || length(cf) != length(PSI0)) next
      sg <- fit$se[[m]]
      if (is.null(sg) || length(sg) != length(PSI0)) sg <- rep(NA_real_, length(PSI0))
      parts[[length(parts) + 1L]] <- data.frame(
        scenario  = scenario,
        family    = family_name,
        n         = n,
        rep       = rep_id,
        estimator = m,
        parameter = PARAM_NAMES,
        estimate  = unname(cf),
        se        = unname(sg),
        stringsAsFactors = FALSE
      )
    }
  }

  add_simex <- function(simex_method, label) {
    ex <- tryCatch(
      extract_simex_K10(
        run_simex_K10(y, z_hat, x_mat, family_obj, Pi_hat,
                      simex_method = simex_method,
                      jackknife = SIMEX_JK),
        K = K),
      error = function(e) NULL
    )
    if (is.null(ex)) return(NULL)
    data.frame(
      scenario  = scenario,
      family    = family_name,
      n         = n,
      rep       = rep_id,
      estimator = label,
      parameter = PARAM_NAMES,
      estimate  = unname(ex$est),
      se        = unname(ex$se),
      stringsAsFactors = FALSE
    )
  }
  if (USE_SIMEX_STD) parts[[length(parts) + 1L]] <- add_simex("standard", "simex_std")
  if (USE_SIMEX_IMP) parts[[length(parts) + 1L]] <- add_simex("improved", "simex_imp")

  parts <- Filter(Negate(is.null), parts)
  if (!length(parts)) return(NULL)
  do.call(rbind, parts)
}

# ---- Parallel dispatcher -------------------------------------------------

run_reps <- function(sim_fn, B, ..., seed) {
  message(sprintf("  dispatching %d reps across %s",
                  B,
                  if (WORKERS > 1L) sprintf("%d workers", WORKERS) else "1 process"))
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
                  as.numeric(Sys.time() - t0, units = "secs"),
                  length(out), B))
  do.call(rbind, out)
}

# ---- Truth and summariser ------------------------------------------------

build_truth_df <- function(scenarios) {
  do.call(rbind, lapply(scenarios, function(sc) {
    data.frame(scenario  = sc,
               parameter = PARAM_NAMES,
               truth     = unname(PSI0),
               stringsAsFactors = FALSE)
  }))
}

summarise_long <- function(raw, truth_df, level = 0.95) {
  z_crit <- qnorm(0.5 + level / 2)
  d <- merge(raw, truth_df, by = c("scenario", "parameter"), sort = FALSE)
  d$err   <- d$estimate - d$truth
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
      bias      = if (length(v))      mean(v) - tru          else NA_real_,
      sd_emp    = if (length(v) > 1L) sd(v)                  else NA_real_,
      mean_se   = if (length(s))      mean(s)                else NA_real_,
      mse       = if (length(v))      mean((v - tru)^2)      else NA_real_,
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
    "K = 10 simulation (B = %d, n = %d, n_val = %d, workers = %d, onestep = %s, simex_std = %s, simex_imp = %s, jackknife = %s)",
    B, N, N_VAL, WORKERS, USE_ONESTEP, USE_SIMEX_STD, USE_SIMEX_IMP, SIMEX_JK
  ))
  message(sprintf("  Pi diag = %.2f, off-diag = %.4f, avg validation per column = %.1f",
                  PI_TRUE_DIAG, (1 - PI_TRUE_DIAG) / (K_DGP - 1L),
                  N_VAL / K_DGP))

  message("\n== Poisson regression ==")
  raw_p <- run_reps(sim_one, B, family_name = "poisson",
                    seed = 20260906L)

  message("\n== Logistic regression ==")
  raw_l <- run_reps(sim_one, B, family_name = "logistic",
                    seed = 20260907L)

  raw <- rbind(raw_p, raw_l)
  truth_df <- build_truth_df(unique(raw$scenario))
  summary  <- summarise_long(raw, truth_df)

  out_csv <- file.path(OUT_DIR, "simulation_K10.csv")
  out_rds <- file.path(OUT_DIR, "simulation_K10.rds")
  write.csv(summary, out_csv, row.names = FALSE)
  saveRDS(list(summary  = summary,
               raw      = raw,
               truth    = truth_df,
               settings = list(B = B, N = N, N_VAL = N_VAL,
                               WORKERS = WORKERS,
                               USE_ONESTEP = USE_ONESTEP,
                               USE_SIMEX_STD = USE_SIMEX_STD,
                               USE_SIMEX_IMP = USE_SIMEX_IMP,
                               SIMEX_JK = SIMEX_JK,
                               K = K_DGP,
                               PI_TRUE = PI_TRUE,
                               PI_PREV_TRUE = PI_PREV_TRUE,
                               GAMMA_TRUE = GAMMA_TRUE,
                               ALPHA_TRUE = ALPHA_TRUE,
                               families = c("poisson", "logistic"),
                               package_version =
                                 as.character(utils::packageVersion("mismeasured"))),
               session  = sessionInfo()),
          out_rds)
  message(sprintf("Summary written to %s\nRaw + summary written to %s",
                  out_csv, out_rds))

  for (sc in unique(summary$scenario)) {
    cat(sprintf("\n=== %s ===\n", sc))
    sub <- summary[summary$scenario == sc,
                   c("estimator", "parameter", "truth", "bias",
                     "sd_emp", "mean_se", "mse", "coverage", "n_succ", "n_se")]
    print(format(sub, digits = 4), row.names = FALSE)
  }
  invisible(list(summary = summary, raw = raw))
}

if (sys.nframe() == 0L) {
  main()
}
