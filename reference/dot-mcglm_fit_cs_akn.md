# Fit AKN98 corrected-score estimator

Works for any K \>= 2 and for Gaussian, Poisson, Binomial families.
Multinomial is unsupported (AKN98 does not cover multinomial responses).

## Usage

``` r
.mcglm_fit_cs_akn(psi_init, y, x_akn, x, K, family, wt = NULL)
```
