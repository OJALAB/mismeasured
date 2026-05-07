# Build xi_hat design matrix for mcglm

For binary Z: xi_hat = cbind(z_hat, x). For multicategory Z: xi_hat =
cbind(d_hat, x) where d_hat is the dummy-encoding of z_hat with baseline
= 0.

## Usage

``` r
.mcglm_build_xi_hat(z_hat, x, K)
```

## Arguments

- z_hat:

  Integer proxy covariate vector.

- x:

  Covariate matrix (n x r).

- K:

  Number of Z categories.

## Value

n x p matrix.
