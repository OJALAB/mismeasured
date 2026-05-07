# Multiplicative bias correction (BCM) for binary misclassification

Multiplicative bias correction (BCM) for binary misclassification

## Usage

``` r
.mcglm_fit_bcm_bin(
  psi_naive,
  y,
  xi_hat,
  x,
  family,
  c1,
  c2,
  iterate = FALSE,
  max_iter = 50,
  tol = 1e-08,
  wt = NULL
)
```
