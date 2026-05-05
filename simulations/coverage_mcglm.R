# =====================================================================
# Coverage simulation for mcglm() Wald confidence intervals
#
# Replicates Scenario I (Poisson with binary misclassification) from
# Battaglia, Christensen, Hansen and Sacher (2025) and reports, for each
# estimation method, the empirical bias, empirical standard deviation,
# mean estimated standard error and 95% Wald CI coverage.
# =====================================================================

suppressPackageStartupMessages({
  library(mismeasured)
})

# --- DGP ---
n      <- 2000
psi0   <- c(gamma = 0.8, alpha0 = -0.5, alpha1 = 0.7)
p01    <- 0.10
p10    <- 0.15
pi_z   <- 0.4
use_onestep <- requireNamespace("RTMB", quietly = TRUE)
methods_fit <- if (use_onestep) {
  c("naive", "bca", "bcm", "cs", "onestep")
} else {
  c("naive", "bca", "bcm", "cs")
}

# We report BCA/BCM under both variance flavours (default = naive sandwich
# matching Theorem 5; "corrected" = joint score-and-drift sandwich).
methods <- c("naive",
             "bca", "bca_corr",
             "bcm", "bcm_corr",
             "cs",
             if (use_onestep) "onestep")

B      <- 500
level  <- 0.95
zcrit  <- qnorm(1 - (1 - level) / 2)

set.seed(20260505)

est <- array(NA_real_, c(B, length(methods), 3),
             dimnames = list(NULL, methods, names(psi0)))
se  <- array(NA_real_, c(B, length(methods), 3),
             dimnames = list(NULL, methods, names(psi0)))
fail <- 0L

t0 <- Sys.time()
for (b in seq_len(B)) {
  x_mat <- cbind(1, rnorm(n))
  z     <- rbinom(n, 1, pi_z)
  eta   <- psi0[1] * z + psi0[2] * x_mat[, 1] + psi0[3] * x_mat[, 2]
  y     <- rpois(n, exp(eta))
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)

  fit_default <- tryCatch(
    mcglm(y, z_hat = z_hat, x = x_mat, family = "poisson",
          method = methods_fit,
          p01 = p01, p10 = p10, pi_z = pi_z,
          vcov_corrected = FALSE,
          fix_omega = TRUE),
    error = function(e) NULL
  )
  fit_corr <- tryCatch(
    mcglm(y, z_hat = z_hat, x = x_mat, family = "poisson",
          method = c("bca", "bcm"),
          p01 = p01, p10 = p10, pi_z = pi_z,
          vcov_corrected = TRUE),
    error = function(e) NULL
  )
  if (is.null(fit_default) || is.null(fit_corr)) { fail <- fail + 1L; next }

  for (m in c("naive", "bca", "bcm", "cs", if (use_onestep) "onestep")) {
    est[b, m, ] <- fit_default$coefficients[[m]]
    se [b, m, ] <- sqrt(diag(fit_default$vcov[[m]]))
  }
  # corrected-variance copies of bca / bcm (point estimate is identical,
  # only the variance differs)
  for (m in c("bca", "bcm")) {
    est[b, paste0(m, "_corr"), ] <- fit_corr$coefficients[[m]]
    se [b, paste0(m, "_corr"), ] <- sqrt(diag(fit_corr$vcov[[m]]))
  }

  if (b %% 50 == 0)
    message(sprintf("  rep %4d / %d  (elapsed: %.1fs)",
                    b, B, as.numeric(Sys.time() - t0, units = "secs")))
}
elapsed <- as.numeric(Sys.time() - t0, units = "secs")
message(sprintf("Done. %d reps, %d failed, %.1fs elapsed (%.2fs/rep).",
                B, fail, elapsed, elapsed / B))

# --- Per-method, per-parameter summary ---
psi0_mat <- matrix(psi0, nrow = B, ncol = 3, byrow = TRUE)

rows <- list()
for (m in methods) {
  bias_m   <- colMeans(est[, m, ] - psi0_mat, na.rm = TRUE)
  sd_emp   <- apply(est[, m, ], 2, sd, na.rm = TRUE)
  se_mean  <- colMeans(se[, m, ], na.rm = TRUE)
  rmse_m   <- sqrt(colMeans((est[, m, ] - psi0_mat)^2, na.rm = TRUE))

  lower <- est[, m, ] - zcrit * se[, m, ]
  upper <- est[, m, ] + zcrit * se[, m, ]
  cover <- (lower <= psi0_mat) & (psi0_mat <= upper)
  cov_rate <- colMeans(cover, na.rm = TRUE)

  for (k in seq_along(psi0)) {
    rows[[length(rows) + 1L]] <- data.frame(
      method     = m,
      parameter  = names(psi0)[k],
      truth      = unname(psi0[k]),
      bias       = bias_m[k],
      sd_emp     = sd_emp[k],
      se_mean    = se_mean[k],
      rmse       = rmse_m[k],
      cov_95     = cov_rate[k],
      stringsAsFactors = FALSE
    )
  }
}
tab <- do.call(rbind, rows)
rownames(tab) <- NULL

cat("\n=== Coverage simulation: Poisson, binary misclassification ===\n")
cat(sprintf("DGP: n = %d, true psi = (%s); p01 = %.2f, p10 = %.2f, pi = %.2f\n",
            n,
            paste(sprintf("%s = %g", names(psi0), psi0), collapse = ", "),
            p01, p10, pi_z))
cat(sprintf("Replications: %d (Wald CI level: %.0f%%)\n", B, 100 * level))
cat(sprintf("Failures: %d\n\n", fail))

print(format(tab, digits = 4), row.names = FALSE)

out_csv <- file.path("simulations", "coverage_mcglm.csv")
write.csv(tab, out_csv, row.names = FALSE)
message(sprintf("Results written to %s", out_csv))

invisible(tab)
