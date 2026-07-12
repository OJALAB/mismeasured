# =====================================================================
# Driver: mass-imputation estimator under drifting misclassification.
#
# Grid (defaults; see README.md):
#   family   : poisson, binomial
#   kappa    : (1,1), (3,2), (6,4)     [p01 = k01/sqrt(n), p10 = k10/sqrt(n)]
#   n_np     : 1000, 4000, 16000, 64000
#   ratio_p  : 1, 0.1                  [n_p = ratio_p * n_np]
#   eta      : known, estimated       [validation subsample, 10% of S_NP]
#
# Environment overrides:
#   SIM_B      number of Monte Carlo replicates per cell (default 500)
#   SIM_CORES  parallel workers (default: all but one)
#   SIM_QUICK  if set (any value), tiny smoke-test grid, B = 4
#
# Usage (from the package root):
#   Rscript simulations/sim-paper-nonprob-misclass/run-simulation.R
#
# Output (in this folder):
#   results-raw.rds   all per-replicate rows
#   results-raw.csv   same, flat CSV
# =====================================================================

suppressPackageStartupMessages({
  library(parallel)
})

sim_dir <- "simulations/sim-paper-nonprob-misclass"
if (!dir.exists(sim_dir)) sim_dir <- "."   # allow running from the folder
source(file.path(sim_dir, "functions.R"))

quick <- nzchar(Sys.getenv("SIM_QUICK"))
B     <- as.integer(Sys.getenv("SIM_B", if (quick) "4" else "500"))
cores <- as.integer(Sys.getenv("SIM_CORES", max(1L, detectCores() - 1L)))

base <- list(
  psi0        = c(gamma = 0.8, alpha0 = -0.5, alpha1 = 0.7),
  pi_theta    = 0.4,
  frac_v      = 0.10,
  N_mult      = 20L,       # N = 20 * n_np (5% sampling fractions)
  dep_theta_x = FALSE
)

if (quick) {
  grid <- expand.grid(
    family  = "poisson",
    kappa   = I(list(c(3, 2))),
    n_np    = 1000L,
    ratio_p = 1,
    eta     = c("known", "estimated"),
    stringsAsFactors = FALSE
  )
} else {
  grid <- expand.grid(
    family  = c("poisson", "binomial"),
    kappa   = I(list(c(1, 1), c(3, 2), c(6, 4))),
    n_np    = c(1000L, 4000L, 16000L, 64000L),
    ratio_p = c(1, 0.1),
    eta     = c("known", "estimated"),
    stringsAsFactors = FALSE
  )
}

message(sprintf("Grid: %d cells x %d reps, %d cores", nrow(grid), B, cores))

set.seed(20260712)
cell_seeds <- sample.int(.Machine$integer.max, nrow(grid))

t0 <- Sys.time()
all_results <- vector("list", nrow(grid))

for (g in seq_len(nrow(grid))) {
  cfg <- c(base, list(
    family  = grid$family[g],
    kappa01 = grid$kappa[[g]][1],
    kappa10 = grid$kappa[[g]][2],
    n_np    = grid$n_np[g],
    ratio_p = grid$ratio_p[g],
    eta     = grid$eta[g]
  ))
  message(sprintf(
    "[%2d/%d] %s kappa=(%g,%g) n_np=%d ratio=%g eta=%s  (elapsed %.1f min)",
    g, nrow(grid), cfg$family, cfg$kappa01, cfg$kappa10,
    cfg$n_np, cfg$ratio_p, cfg$eta,
    as.numeric(Sys.time() - t0, units = "mins")))

  res <- mclapply(seq_len(B), function(b) {
    set.seed(cell_seeds[g] %% 1000000L * 1000L + b)
    tryCatch(run_one_rep(cfg, b), error = function(e) NULL)
  }, mc.cores = cores)

  ok <- Filter(Negate(is.null), res)
  if (length(ok) < B)
    message(sprintf("    %d / %d replicates failed", B - length(ok), B))
  if (length(ok)) all_results[[g]] <- do.call(rbind, ok)
}

results <- do.call(rbind, Filter(Negate(is.null), all_results))
rownames(results) <- NULL

saveRDS(results, file.path(sim_dir, "results-raw.rds"))
write.csv(results, file.path(sim_dir, "results-raw.csv"), row.names = FALSE)
message(sprintf("Done in %.1f min. %d rows written to %s",
                as.numeric(Sys.time() - t0, units = "mins"),
                nrow(results), file.path(sim_dir, "results-raw.rds")))
