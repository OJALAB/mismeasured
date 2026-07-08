# Jackknife variance estimation for SIMEX (continuous measurement error)

Computes Var_model (mean of model-based vcov across B refits) minus
Var_sim (empirical covariance of B coefficient estimates) at each
lambda, then extrapolates to lambda = -1.

## Usage

``` r
jackknife_variance_simex(
  theta_list,
  naive_fit,
  lambda,
  jk_method,
  vcov_model_list = NULL
)
```

## Details

This follows the Stefanski & Cook (1995) measurement-error jackknife,
where Var_model at each lambda is the average of vcov(glm_b) across the
B simulation replicates. The per-lambda model vcov averages are computed
in C++ and passed in via `vcov_model_list`.
