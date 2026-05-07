# Estimate pi_z (true prevalence) from observed proxy and Pi

Uses Bayesian inversion: if pi_obs = Pi pi_z = solve(Pi)

## Usage

``` r
.mcglm_estimate_pi_z(z_hat, Pi)
```

## Arguments

- z_hat:

  Integer vector of observed proxy values (0-based).

- Pi:

  K x K misclassification matrix.

## Value

Numeric vector of length K (estimated true prevalences).
