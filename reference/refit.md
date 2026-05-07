# Refit a SIMEX model with a different extrapolation function

Re-extrapolates the stored simulation results using a different fitting
method, without re-running the simulation step.

## Usage

``` r
refit(object, ...)
```

## Arguments

- object:

  a `simex` object.

- ...:

  additional arguments, including `extrapolation` (`"quadratic"`,
  `"linear"`, or `"loglinear"`) and `jackknife` (logical or character).

## Value

An updated `simex` object with new extrapolation results.
