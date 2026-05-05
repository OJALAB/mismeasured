# mismeasured 0.5.2

## Bug fixes

* **Improved MC-SIMEX no longer crashes on asymmetric \eqn{\Pi}**: the
  internal helper `.mat_power_r()` raised
  \emph{`unimplemented complex function`} whenever the eigendecomposition
  of \eqn{\Pi} produced complex eigenvalues. This was triggered by
  K-level estimated-\eqn{\Pi} workflows (e.g. multicategory misclassification
  with a small validation sample). Two underlying bugs were fixed:
  - `sign(d)` was called on complex eigenvalues, which is undefined.
  - `abs(d)^power * sign(d)` is not equal to `d^power` for negative or
    complex `d`; the principal branch `d^power` is now used instead.
  Integer powers (e.g. \eqn{\Pi^2}, the only one needed for the
  default lambda = 1 in improved MC-SIMEX) are now computed by exact
  matrix multiplication, avoiding eigendecomposition entirely. In a
  reproducer (Scenario III of the paper, K = 4, n = 5000, validation
  subsample of 500), the number of successful Monte Carlo replicates
  went from 1/30 to 30/30.

# mismeasured 0.5.1

## New features

* **Asymptotic inference for `mcglm()`**: each fitted method (`naive`,
  `bca`, `bcm`, `cs`, `onestep`) now carries its own asymptotic
  variance-covariance matrix and standard errors. Variance estimators
  follow Battaglia, Christensen, Hansen and Sacher (2025):
  the drifting-regime sandwich \eqn{A^{-1} C A^{-1}} for naive/BCA/BCM
  (Theorems on two-step and BC), the Z-estimator sandwich
  \eqn{J^{-1} S J^{-T}} for CS, and the inverse Hessian for the
  one-step mixture-likelihood estimator.

* **`vcov_corrected` argument**: when set to `TRUE`, the BCA/BCM
  variance accounts for the additional uncertainty due to the estimated
  drift correction (joint score-and-drift sandwich). The default
  (`FALSE`) returns the asymptotic variance from the paper's theorems.

* **glm-style S3 methods for `mcglm` objects**:
  `summary()` (Wald table per method with z-values and p-values),
  `print()` (call / family / coefficient table / null + residual
  deviance / AIC), `vcov()`, `coef()`, `confint()`, `fitted()`,
  `predict()`, `residuals()` (response/pearson/deviance), `logLik()`
  (naive and onestep), `AIC()`, `BIC()`, `nobs()`, `family()`,
  `formula()`, `model.matrix()`, plus `se.mcglm()` for standard errors.
  All accept a `method = ...` argument to select the estimator.

# mismeasured 0.5.0

## Bug fixes

* **`RcppExports.R` invalid defaults**: `compileAttributes()` generated
  `integerVector(0)`, `numericVector(0)`, and `matrix(0, 0)` as default
  argument values in the R wrappers for `mcsimex_multi_sim_cpp` and
  `simex_sim_cpp`. These are not valid R expressions, causing runtime errors
  when the functions were called with default arguments. Fixed to `integer(0)`,
  `numeric(0)`, and `matrix(0, nrow = 0, ncol = 0)`.

## Documentation

* Updated README to use the `mcglm()` formula interface
  (`y ~ mc(z, Pi) + x1`) instead of the matrix interface.

# mismeasured 0.4.0

## New features

* **`mcglm()` for bias-corrected GLMs**: added `mcglm()` implementing
  analytical bias correction for GLMs with misclassified covariates following
  Battaglia, Christensen, Hansen and Sacher (2025). Supports five estimation
  methods (naive, BCA, BCM, corrected score, one-step), binary and multicategory
  misclassified covariates, Poisson/Binomial/Gaussian families, and multinomial
  response models.

* **Formula interface for `mcglm()`**: `y ~ mc(z, Pi) + x1 + x2` syntax with
  automatic extraction of the misclassification matrix and optional auto-estimation
  of true prevalence `pi_z` from observed data and Pi.

# mismeasured 0.3.0

## New features

* **K-level improved MC-SIMEX**: `method = "improved"` now supports a single
  K-level misclassified covariate with a fixed, known misclassification matrix.
  The implementation estimates latent category probabilities from observed
  category frequencies, builds the K-level dummy-coefficient correction matrix,
  and applies the corresponding intercept adjustment. Binary fits keep the
  existing `c.lambda` output; K-level fits expose `pi.vec` and
  `correction.matrix`.

# mismeasured 0.2.1

## Bug fixes

* **MC design-matrix construction preserves formula semantics**: the MC-SIMEX
  path previously used `all.vars()` + `reformulate()` to build the non-mc
  covariate matrix, which silently dropped `I()` transformations, interactions,
  offsets, and explicit intercept choices. For example, `y ~ mc(z, Pi) + I(x^2)`
  was fitted as `y ~ z + x`. Now uses `model.matrix()` on the clean formula and
  removes mc-variable columns via the `assign` attribute, preserving all formula
  semantics. (#6)

* **`predict.simex()` factor level mismatch**: predicting with `newdata`
  containing only a subset of mc factor levels (e.g., only level `"1"`) could
  misassign the baseline, producing wrong predictions. Now enforces training
  factor levels stored on the fitted object. The multi-mc predict path also
  switched from `model.matrix()` (wrong column order) to the same manual dummy
  construction used during fitting. (#6)

* **`refit()` design-matrix mismatch for multi-MC models**: `refit()` used
  `naive.model$x` for multi-mc and response+mc fits, but that matrix has
  different column ordering than the package's manual `[dummies | x_mat]`
  layout. Fitted values and variance estimates could be computed against the
  wrong columns. Now uses stored `xi.hat` directly. (#6)

## Improvements

* **MC matrix validation**: misclassification matrices now reject entries
  outside `[0, 1]` (previously only checked column sums). (#6)

## Packaging

* Removed tracked object files (`src/*.o`) from git.
* Added `theory/`, `src/symbols.rds`, `.DS_Store` to `.Rbuildignore`/`.gitignore`.
* Lowered `CXX_STD` from `CXX17` to `CXX11` for broader platform portability.

# mismeasured 0.2.0

## Documentation

* Added package vignette "Correcting measurement error and misclassification
  with mismeasured" with 10 worked examples covering all error types, families,
  and features.

## New features

* **Multiple `mc()` terms**: formulas like `y ~ mc(z1, Pi1) + mc(z2, Pi2) + x`
  are now supported with the standard MC-SIMEX method. Each misclassified
  covariate is resampled independently. A new C++ function
  `mcsimex_multi_sim_cpp` handles the multi-variable simulation.

* **Response misclassification**: `mc()` can now be used on the left-hand side
  of the formula, e.g. `mc(y, Pi_y) ~ mc(z, Pi_z) + x`. The response is
  resampled alongside covariates at each simulation replicate. Supports
  response-only mc, or combined response + covariate mc. Runs in C++.

* **Berkson measurement error**: `me(x, sd, type = "berkson")` implements the
  Berkson SIMEX correction, where the simulation step subtracts noise instead
  of adding it.

* **Non-zero mean measurement error**: `me(x, sd, mean = mu)` supports
  systematic (non-zero mean) measurement error in the SIMEX simulation.

* **Frequency weight tests**: added test suite for the `weights` argument
  covering expanded-data equivalence and family-specific cases.

## Bug fixes
 
* **ME-SIMEX jackknife variance**: the variance estimator for `me()` terms was
  using the naive model's vcov at all lambda levels instead of the mean of
  model-based vcov across B simulation refits. This caused severe SE
  underestimation (~63% CI coverage). Fixed by computing `(X_sim' W X_sim)^{-1}`
  in C++ at each replicate and accumulating the average. SE ratio vs CRAN
  `simex` is now ~1.0.

* **`glm()` weights NSE bug**: passing `weights` to `glm()` inside
  `.simex_continuous()` and `.simex_discrete()` failed because `glm()` uses
  non-standard evaluation and found `stats::weights` (a function) instead of
  the local variable. Fixed by storing weights as `.wt` in the data frame.

## Internal

* `mcsimex_multi_sim_cpp` extended with optional response mc parameters
  (`y_z_hat_r`, `Pi_y_r`, `K_y`) for C++-level response resampling.
* `simex_sim_cpp` extended with `error_type_r` and `me_mean_r` parameters
  for Berkson and non-zero mean error.
* `simex_sim_cpp` now returns a list with `theta` (coefficient matrix) and
  `vcov_model` (per-lambda average model vcov) for correct jackknife variance.
* `refit()` refactored via `.refit_build_xi_hat()` to handle single-mc,
  multi-mc, and response-mc configurations.
