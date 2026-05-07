# Wald confidence intervals for an `mcglm` fit

Wald confidence intervals for an `mcglm` fit

## Usage

``` r
# S3 method for class 'mcglm'
confint(object, parm, level = 0.95, method = NULL, ...)
```

## Arguments

- object:

  An `mcglm` object.

- parm:

  Optional names or indices of parameters to return.

- level:

  Coverage level.

- method:

  Optional estimation method (defaults to the last one fit, e.g.
  `"onestep"` if present).

- ...:

  Unused.
