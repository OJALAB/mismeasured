# Print method for `mcglm` fits

Mirrors the layout of `print.glm`: shows the call, the family, a
coefficient table (one column per estimation method), the residual
degrees of freedom and (when available) the deviance and AIC of the
naive GLM fit.

## Usage

``` r
# S3 method for class 'mcglm'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `"mcglm"`.

- digits:

  Number of digits to display.

- ...:

  Unused.
