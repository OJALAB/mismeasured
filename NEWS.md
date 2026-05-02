# mismeasured 0.2.0

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
