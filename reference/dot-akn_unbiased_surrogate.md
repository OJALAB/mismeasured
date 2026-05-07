# Build AKN98 unbiased surrogate

Build AKN98 unbiased surrogate

## Usage

``` r
.akn_unbiased_surrogate(z_hat, Q_inv, p0, K)
```

## Value

n x (K-1) matrix x_akn with E\[x_akn\[i, \] \| Z_i = l\] = e_l for l =
1, ..., K-1 and = 0 for l = 0.
