# Sandwich variance for BCA/BCM estimators (multicategory)

Sandwich variance for BCA/BCM estimators (multicategory)

## Usage

``` r
.mcglm_vcov_bc_multi(
  psi_bc,
  y,
  xi_hat,
  z_hat,
  x,
  K,
  family,
  Pi = NULL,
  pi_z = NULL,
  psi_naive = NULL,
  type = c("bca", "bcm"),
  corrected = FALSE,
  wt = NULL,
  jacobian = c("analytical", "numerical")
)
```
