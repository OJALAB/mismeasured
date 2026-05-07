# Compute matrix power Pi^power

Uses exact matrix multiplication for non-negative integer powers and the
principal branch of the eigendecomposition otherwise. Pi need not be
symmetric; complex eigenvalues are handled via complex arithmetic (the
result is real for any column-stochastic Pi raised to a real power,
since complex eigenvalues come in conjugate pairs).

## Usage

``` r
.mat_power_r(Pi, power)
```
