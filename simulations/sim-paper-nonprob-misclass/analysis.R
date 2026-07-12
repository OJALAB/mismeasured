# =====================================================================
# Summaries and figures for the mass-imputation drifting-regime study.
#
# Reads results-raw.rds (produced by run-simulation.R) and writes:
#   summary.csv          per-cell, per-method metrics
#   fig-scaled-bias.pdf  sqrt(n_NP) * bias vs n_NP -- the drifting-regime
#                        diagnostic: naive should stay flat at a nonzero
#                        level, BCA/BCM/CS/oracle should converge to 0
#   fig-coverage.pdf     95% CI coverage vs n_NP
# =====================================================================

sim_dir <- "simulations/sim-paper-nonprob-misclass"
if (!dir.exists(sim_dir)) sim_dir <- "."

res <- readRDS(file.path(sim_dir, "results-raw.rds"))

zcrit <- qnorm(0.975)
res$err     <- res$ybar_hat - res$ybar_true
res$covered <- abs(res$err) <= zcrit * res$se_hat

cells <- unique(res[, c("family", "kappa01", "kappa10",
                        "n_np", "ratio_p", "eta", "method", "estimator")])

summ <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
  key <- cells[i, ]
  sub <- merge(res, key)
  data.frame(
    key,
    n_reps      = nrow(sub),
    bias        = mean(sub$err),
    bias_scaled = sqrt(key$n_np) * mean(sub$err),
    mc_se_bias  = sd(sub$err) / sqrt(nrow(sub)),
    sd_emp      = sd(sub$ybar_hat),
    se_mean     = mean(sub$se_hat),
    rmse        = sqrt(mean(sub$err^2)),
    cov_95      = mean(sub$covered),
    share_V_P   = mean(sub$V_P / (sub$V_P + sub$V_NP + sub$V_V + sub$V_cov)),
    row.names   = NULL
  )
}))
summ <- summ[order(summ$family, summ$kappa01, summ$ratio_p,
                   summ$eta, summ$estimator, summ$n_np, summ$method), ]

write.csv(summ, file.path(sim_dir, "summary.csv"), row.names = FALSE)
cat("Wrote", file.path(sim_dir, "summary.csv"), "\n\n")
print(format(summ, digits = 3), row.names = FALSE)

# ---- figures ---------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  summ$kappa <- sprintf("kappa = (%g, %g)", summ$kappa01, summ$kappa10)
  summ$arm   <- sprintf("%s, n_P/n_NP = %g, eta %s, %s",
                        summ$family, summ$ratio_p, summ$eta,
                        summ$estimator)

  p1 <- ggplot(summ,
               aes(n_np, bias_scaled, colour = method, group = method)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
    geom_line() + geom_point() +
    geom_errorbar(aes(ymin = bias_scaled - 2 * sqrt(n_np) * mc_se_bias,
                      ymax = bias_scaled + 2 * sqrt(n_np) * mc_se_bias),
                  width = 0.05) +
    scale_x_log10() +
    facet_grid(arm ~ kappa, scales = "free_y") +
    labs(x = expression(n[NP]),
         y = expression(sqrt(n[NP]) %*% bias),
         title = "Mass imputation under drifting misclassification",
         subtitle = paste("Theory: naive flat at a nonzero level;",
                          "BCA/BCM/CS converge to 0")) +
    theme_bw()
  ggsave(file.path(sim_dir, "fig-scaled-bias.pdf"), p1,
         width = 11, height = 2.2 * length(unique(summ$arm)) + 1,
         limitsize = FALSE)

  p2 <- ggplot(summ, aes(n_np, cov_95, colour = method, group = method)) +
    geom_hline(yintercept = 0.95, linetype = 2, colour = "grey50") +
    geom_line() + geom_point() +
    scale_x_log10() +
    facet_grid(arm ~ kappa) +
    labs(x = expression(n[NP]), y = "95% CI coverage",
         title = "Coverage of the plug-in variance V_P + V_NP (+ V_V)") +
    theme_bw()
  ggsave(file.path(sim_dir, "fig-coverage.pdf"), p2,
         width = 11, height = 2.2 * length(unique(summ$arm)) + 1,
         limitsize = FALSE)
  cat("\nWrote fig-scaled-bias.pdf and fig-coverage.pdf\n")
} else {
  cat("\nggplot2 not installed; skipping figures.\n")
}
