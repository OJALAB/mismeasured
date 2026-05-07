# BCM for multicategory misclassification

BCM for multicategory misclassification

## Usage

``` r
.mcglm_fit_bcm_multi(
  psi_naive,
  y,
  xi_hat,
  z_hat,
  x,
  K,
  family,
  Pi,
  pi_z,
  iterate = FALSE,
  max_iter = 50,
  tol = 1e-08,
  wt = NULL,
  jacobian = c("analytical", "numerical")
)
```
