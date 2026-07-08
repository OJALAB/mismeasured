# Resample misclassified categories at inflation level lambda

Draws z\* \| z_hat from the columns of Pi^lambda, replicating the
simulation step of MC-SIMEX so that a refit on z\* estimates
Var(theta(lambda)) rather than Var(theta(0)).

## Usage

``` r
.resample_z_lambda(z_hat, Pi, lambda, K)
```
