# Specify a variable measured with error

Use inside a
[`simex`](https://ojalab.github.io/mismeasured/reference/simex.md)
formula to mark a covariate (or response) as measured with additive
Gaussian error.

## Usage

``` r
me(variable, sd, type = "classical", mean = 0)
```

## Arguments

- variable:

  bare name of the variable in the data.

- sd:

  measurement error standard deviation. Can be a scalar (homoscedastic),
  a numeric vector of length n, or a bare column name in the data
  (heteroscedastic).

- type:

  type of measurement error: `"classical"` (default) for the standard
  additive error model \\W = X + U\\, or `"berkson"` for the Berkson
  error model \\X = W + U\\ where the true value is a noisy version of
  the observed.

- mean:

  mean of the measurement error distribution (default 0). Use for
  systematic (non-zero mean) measurement error.

## Value

An object of class `"me_term"` (used internally by
[`simex`](https://ojalab.github.io/mismeasured/reference/simex.md); not
intended to be called directly).

## See also

[`mc`](https://ojalab.github.io/mismeasured/reference/mc.md),
[`simex`](https://ojalab.github.io/mismeasured/reference/simex.md)

## Examples

``` r
if (FALSE) { # \dontrun{
simex(y ~ me(x, 0.5) + w, data = df)
simex(y ~ me(x, sd_x) + w, data = df)            # heteroscedastic
simex(y ~ me(x, 0.5, type = "berkson") + w, data = df)  # Berkson error
simex(y ~ me(x, 0.5, mean = 0.1) + w, data = df)       # non-zero mean
} # }
```
