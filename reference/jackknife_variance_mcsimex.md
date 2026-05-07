# Jackknife variance estimation for MC-SIMEX

Computes Var_jk = Var_model - Var_sim at each lambda, then extrapolates
each variance element to lambda = -1.

## Usage

``` r
jackknife_variance_mcsimex(
  theta_list,
  psi_naive,
  naive_fit,
  lambda,
  jk_method,
  fam,
  xi_hat,
  y,
  wt
)
```
