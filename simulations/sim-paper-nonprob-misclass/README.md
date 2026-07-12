# Mass imputation under drifting misclassification

Simulation study for the corrected-score mass-imputation estimator
(`theory/ml-glm/paper.tex`, Section "Corrected-score mass imputation
estimator", Theorem `thm:mass-imp-cs`) under the **drifting regime**
`p01 = kappa01 / sqrt(n_NP)`, `p10 = kappa10 / sqrt(n_NP)`.

## Design

- **Population**: `N = 20 * n_NP`, `x ~ N(0,1)`, `theta ~ Bern(0.4)`
  independent of `x`, `Y ~ GLM(mu(gamma*theta + alpha0 + alpha1*x))`,
  canonical link (Poisson / logistic), `psi0 = (0.8, -0.5, 0.7)` as in
  the fixed-regime study (`simulations/coverage_mcglm.R`).
- **S_NP** (non-probability): informative selection on `x` only, not
  on `Y` (transportability holds) and — deliberately — not on `theta`:
  selecting on `(x, theta)` jointly would induce dependence between
  `theta` and `x` within `S_NP` (collider stratification) and shift
  the in-sample prevalence, violating the `theta` independent of `q`
  assumption behind `c1`/`c2`. Selection on `theta` (`cfg$sel_theta`)
  is available as an assumption-violation arm. `S_NP` carries `Y` and
  the misclassified proxy `theta_hat`.
- **S_P** (probability): Poisson sampling, `pi_i` proportional to
  `exp(0.3 x)`, known weights. Carries the true `xi = (theta, 1, x)`.
  Poisson sampling makes joint inclusion probabilities factorize, so
  the paper's HT variance estimator reduces to a single sum.
- **Estimators**: mass-imputation mean with `psi_hat` from `mcglm()`
  (naive, BCA, BCM, CS) plus an oracle fit on the true `theta`. Each
  is computed in two forms: `ht` (`1/N`, the form in Theorem
  `thm:mass-imp-cs` / Kim et al. 2021) and `hajek` (`1/sum(w)`,
  ratio-stabilized). Under Poisson sampling the realized `sum(w)` is
  random, so the HT form carries an extra variance component
  proportional to the level of `mu`; Hajek removes it and typically
  has a much smaller `V_P`, which also makes the misclassification
  bias visible at smaller `n_P/n_NP` ratios.
- **Misclassification parameters**: known, or estimated from a
  validation subsample `S_V` (10% of `S_NP`). Under drift the number
  of misclassified units in `S_V` is only `O(sqrt(n_NP))`, so the
  `o_P(n^{-1/2})` condition of Prop. `prop:est_prob` sits at the
  boundary — degradation in this arm is a finding, not a bug.
- **Key diagnostic** (`fig-scaled-bias.pdf`): `sqrt(n_NP) * bias` vs
  `n_NP` with `kappa` fixed. Theory predicts a flat nonzero line for
  the naive estimator and convergence to 0 for BCA/BCM/CS.
- **Sample-size ratio**: `n_P/n_NP` in {1, 0.1}. With `n_P << n_NP`
  the `O(n_NP^{-1/2})` misclassification bias is dominated by the
  `O(n_P^{-1/2})` design noise (`V ~ V_P`), so raw-bias effects are
  visible mainly in the `ratio = 1` arm.

## Variance estimator — implemented here, not in the package

The `mismeasured` package does **not** provide a variance estimator
for the mass-imputation estimator. `functions.R::mi_variance()`
implements the paper's plug-in formulas directly:

- `V_P`: design variance under Poisson sampling — for `ht` the paper's
  HT form with *uncentered* `mu_i^2` (correct for the `1/N` estimator:
  it includes the random-sample-size component); for `hajek` the
  standard ratio linearization with residuals `mu_i - Ybar_hat` and
  `N_hat = sum(w)`. Note the paper states the linearization only for
  the `1/N` form; the Hajek variant is the standard ratio extension
  and should be flagged as such if reported in the text;
- `V_NP = G_P' A^{-1} S (A^{-1})' G_P / n_NP`, with the naive sandwich
  (`A = I_hat`, `S` = outer score) for naive/BCA/BCM — their
  asymptotic variance in the drifting regime per Theorem `thm:mis-bc`
  — and the corrected-score form (`A = I_hat + M_hat`, `S` = outer
  `phi`) for CS;
- for CS with estimated probabilities: `V_V` (validation-sample term)
  and the `Cov(phi_i, s_i)` term (because `S_V` is a subsample of
  `S_NP`), from the subsection "Corrected-score under estimated
  misclassification probabilities".

If these formulas are later added to the package, this file is the
reference implementation to check against.

## Files

| file | purpose |
|---|---|
| `functions.R` | DGP, samples, MI point estimator, inline variance estimator, one-replicate runner |
| `run-simulation.R` | grid driver (parallel), writes `results-raw.{rds,csv}` |
| `analysis.R` | writes `summary.csv`, `fig-scaled-bias.pdf`, `fig-coverage.pdf` |

## Running

From the package root:

```bash
# smoke test (~1 min): tiny grid, B = 4
SIM_QUICK=1 Rscript simulations/sim-paper-nonprob-misclass/run-simulation.R

# full run (long; use a server): 96 cells x 500 reps
SIM_B=500 SIM_CORES=16 Rscript simulations/sim-paper-nonprob-misclass/run-simulation.R

Rscript simulations/sim-paper-nonprob-misclass/analysis.R
```

The full grid is `2 families x 3 kappa x 4 n_NP x 2 ratios x 2 eta
arms`. The `n_NP = 64000` cells dominate the runtime; trim the grid in
`run-simulation.R` for a first pass (e.g. drop `binomial` or cap
`n_NP` at 16000).

Monte Carlo precision: `summary.csv` reports `mc_se_bias`; check that
it is small relative to `bias` before reading anything into the
scaled-bias plot (the analysis script draws +-2 MC-SE error bars).
