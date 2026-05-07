# Dispatch wrapper: analytical (default) or numerical Jacobian for K-class drift

Dispatch wrapper: analytical (default) or numerical Jacobian for K-class
drift

## Usage

``` r
.mcglm_compute_Mhat_multi(
  psi,
  x,
  K,
  mu_fun,
  Pi,
  pi_z,
  wt = NULL,
  jacobian = c("analytical", "numerical"),
  mu_dot_fun = NULL
)
```
