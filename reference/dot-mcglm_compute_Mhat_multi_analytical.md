# Analytical Jacobian M_hat for multicategory drift

Closed-form block expression for the K-class Jacobian, with B = Pi
column-scaled by pi_z and w = rowSums(B) = Pr(Z_hat = j). Replaces the
forward-difference numerical Jacobian.

## Usage

``` r
.mcglm_compute_Mhat_multi_analytical(
  psi,
  x,
  K,
  mu_dot_fun,
  Pi,
  pi_z,
  wt = NULL
)
```
