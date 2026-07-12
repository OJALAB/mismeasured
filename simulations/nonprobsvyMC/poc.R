# =====================================================================
# Proof of concept for nonprobsvyMC (see PLAN.md, Option A / M0).
#
# Demonstrates the full nonprob_mc() pipeline against nonprobsvy:
#   1. DGP: finite population, informative non-probability sample with
#      a misclassified binary covariate, probability sample (SRSWOR)
#      with the true covariate and known design.
#   2. Reference: nonprobsvy::nonprob() MI with the PROXY covariate
#      (what a user can do today) -- biased.
#   3. Bridge: parse `y ~ mc(z, Pi) + x`, fit mismeasured::mcglm(),
#      predict on the svydesign, svymean() -> point + var_prob,
#      corrected-score var_nonprob (Theorem thm:mass-imp-cs).
#   4. Oracle: nonprob() MI with the TRUE covariate.
#
# Expected output: naive ~ biased; mcglm-cs ~ oracle; SEs comparable.
# =====================================================================

suppressPackageStartupMessages({
  library(nonprobsvy)   # >= 0.3.0
  library(mismeasured)
  library(survey)
})

set.seed(20260712)

# ---- 1. DGP ----------------------------------------------------------
N    <- 100000
n_np <- 5000    # expected non-probability sample size
n_p  <- 1000    # probability sample size (SRSWOR)

x1    <- rnorm(N)
z     <- rbinom(N, 1, 0.4)
psi0  <- c(gamma = 0.8, alpha0 = -0.5, alpha1 = 0.7)
y     <- rpois(N, exp(psi0[1] * z + psi0[2] + psi0[3] * x1))
ybar_true <- mean(y)

# fixed misclassification, strong enough to see the bias clearly
p01 <- 0.15; p10 <- 0.20
Pi  <- matrix(c(1 - p01, p01, p10, 1 - p10), 2, 2)  # Pi[j,l]=P(zhat=j-1|z=l-1)
z_hat <- z
z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01)
z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)

# informative NP selection on x only (see sim-paper-nonprob-misclass:
# selecting on z would break the z-independent-of-q assumption in S_NP)
a0  <- uniroot(function(a) mean(plogis(a + 0.5 * x1)) - n_np / N,
               c(-20, 5))$root
s_np   <- rbinom(N, 1, plogis(a0 + 0.5 * x1)) == 1
df_np  <- data.frame(y = y, z = z_hat, x1 = x1)[s_np, ]   # proxy z!

# probability sample: SRSWOR, true z observed, y NOT observed
s_p   <- sample.int(N, n_p)
df_p  <- data.frame(z = z[s_p], x1 = x1[s_p])
design <- svydesign(ids = ~1, data = df_p, weights = rep(N / n_p, n_p),
                    fpc = rep(N, n_p))

# ---- 2. reference: nonprobsvy naive MI with the proxy ---------------
fit_naive <- nonprob(
  outcome        = y ~ z + x1,
  data           = df_np,
  svydesign      = design,
  method_outcome = "glm",
  family_outcome = "poisson"
)

# ---- 3. the nonprob_mc() pipeline ------------------------------------
# Everything below uses only PUBLIC API (mismeasured >= 0.7.0):
# mc_parse_formula(), mcglm() formula interface, predict(newdata),
# sandwich::estfun/bread/meat. Zero ':::' -- the CRAN-compliance
# acceptance check for nonprobsvyMC (PLAN.md Sections 3-4).

# (a) the mc() parser nonprobsvyMC uses to sanitize the outcome formula
#     before nonprobsvy's own formula machinery ever sees it
parsed <- mc_parse_formula(y ~ mc(z, Pi) + x1, df_np)
stopifnot(parsed$K == 2L, identical(parsed$Pi, Pi))

# (b) corrected fit on the non-probability sample (formula interface)
fit_mc <- mcglm(y ~ mc(z, Pi) + x1, data = df_np, family = "poisson",
                method = c("naive", "bca", "bcm", "cs"))
psi_cs <- coef(fit_mc, method = "cs")

# (c) predict on the probability sample (TRUE z in df_p$z) and svymean()
design <- update(design,
                 y_pred = predict(fit_mc, newdata = df_p, type = "response"))
mn     <- svymean(~y_pred, design)
var_prob <- as.numeric(attr(mn, "var"))

# (d) var_nonprob: corrected-score linearization (thm:mass-imp-cs)
#     A^{-1} = bread, S = meat, G_P from the design weights
fam   <- stats::poisson()
B     <- sandwich::bread(fit_mc, method = "cs")     # (I_hat + M_hat)^{-1}
S     <- sandwich::meat(fit_mc, method = "cs")      # mean phi phi'
xi_p  <- cbind(df_p$z, 1, df_p$x1)                  # (z, intercept, x1)
w     <- weights(design)
eta_p <- predict(fit_mc, newdata = df_p, type = "link")
G_P   <- colSums(w * fam$mu.eta(eta_p) * xi_p) / sum(w)
var_nonprob <- as.numeric(t(G_P) %*% B %*% S %*% t(B) %*% G_P) /
  nobs(fit_mc)

est_mc <- as.numeric(mn)
se_mc  <- sqrt(var_prob + var_nonprob)

# ---- 4. oracle: nonprob() with the true z ----------------------------
df_np_true <- data.frame(y = y, z = z, x1 = x1)[s_np, ]
fit_oracle <- nonprob(
  outcome        = y ~ z + x1,
  data           = df_np_true,
  svydesign      = design,
  method_outcome = "glm",
  family_outcome = "poisson"
)

# ---- report ----------------------------------------------------------
row <- function(name, est, se)
  sprintf("%-22s %8.4f %9.4f %9.4f", name, est, est - ybar_true, se)
cat(sprintf("true population mean:  %.4f\n\n", ybar_true))
cat(sprintf("%-22s %8s %9s %9s\n", "estimator", "estimate", "bias", "SE"))
cat(row("nonprob glm (proxy z)",
        fit_naive$output$mean, fit_naive$output$SE), "\n")
cat(row("nonprob_mc mcglm-cs", est_mc, se_mc), "\n")
cat(row("nonprob glm (true z)",
        fit_oracle$output$mean, fit_oracle$output$SE), "\n")
cat(sprintf("\nvar decomposition (mcglm-cs): var_prob = %.3e, var_nonprob = %.3e\n",
            var_prob, var_nonprob))
cat(sprintf("psi_cs = (%s)  vs  psi0 = (%s)\n",
            paste(round(psi_cs, 3), collapse = ", "),
            paste(round(psi0, 3), collapse = ", ")))

# ---- 5. Monte Carlo: where the correction matters --------------------
# For the OVERALL mean the intercept absorbs most of the
# misclassification bias (naive ~ -0.005 vs cs ~ +0.002 in this DGP).
# The correction bites on DOMAIN means: predicting the z = 1
# subpopulation with the attenuated naive gamma is badly biased.
# Set POC_B=0 to skip (default B = 200, ~1-2 min on 6 cores).
B <- as.integer(Sys.getenv("POC_B", "200"))
if (B > 0) {
  suppressPackageStartupMessages(library(parallel))
  fam <- stats::poisson()
  one <- function(b) {
    set.seed(b)
    x1 <- rnorm(N); z <- rbinom(N, 1, 0.4)
    y  <- rpois(N, exp(psi0[1] * z + psi0[2] + psi0[3] * x1))
    zh <- z
    zh[z == 0] <- rbinom(sum(z == 0), 1, p01)
    zh[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)
    a0 <- uniroot(function(a) mean(plogis(a + 0.5 * x1)) - n_np / N,
                  c(-20, 5))$root
    s  <- rbinom(N, 1, plogis(a0 + 0.5 * x1)) == 1
    sp <- sample.int(N, n_p)
    fit <- mcglm(y[s], z_hat = zh[s], x = cbind(1, x1[s]),
                 family = "poisson", method = c("naive", "cs"), Pi = Pi)
    d1  <- z[sp] == 1
    xi_all <- cbind(z[sp], 1, x1[sp])
    xi_d1  <- xi_all[d1, , drop = FALSE]
    c(all_naive = mean(fam$linkinv(xi_all %*% fit$coefficients$naive)),
      all_cs    = mean(fam$linkinv(xi_all %*% fit$coefficients$cs)),
      all_true  = mean(y),
      d1_naive  = mean(fam$linkinv(xi_d1 %*% fit$coefficients$naive)),
      d1_cs     = mean(fam$linkinv(xi_d1 %*% fit$coefficients$cs)),
      d1_true   = mean(y[z == 1]))
  }
  r <- do.call(rbind, mclapply(seq_len(B), one,
                               mc.cores = max(1L, detectCores() - 2L)))
  cat(sprintf("\n--- Monte Carlo, B = %d ---\n", B))
  cat(sprintf("overall mean bias: naive %+ .4f | cs %+ .4f  (MC-SE %.4f)\n",
              mean(r[, 1] - r[, 3]), mean(r[, 2] - r[, 3]),
              sd(r[, 1] - r[, 3]) / sqrt(B)))
  cat(sprintf("domain z=1  bias: naive %+ .4f | cs %+ .4f  (MC-SE %.4f)\n",
              mean(r[, 4] - r[, 6]), mean(r[, 5] - r[, 6]),
              sd(r[, 4] - r[, 6]) / sqrt(B)))
}
