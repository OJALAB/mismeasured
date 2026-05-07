# AKN98 corrected information

Returned matrix is the \*mean\* per observation; the Jacobian of the
mean score colMeans(phi) wrt psi is -.akn_information() (negative
because the score is residual-based).

## Usage

``` r
.akn_information(psi, y, x_akn, x, K, fam, wt = NULL)
```
