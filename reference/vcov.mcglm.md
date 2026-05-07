# Variance-covariance matrix of an `mcglm` fit

Variance-covariance matrix of an `mcglm` fit

## Usage

``` r
# S3 method for class 'mcglm'
vcov(object, method = NULL, ...)
```

## Arguments

- object:

  An `mcglm` object.

- method:

  Optional character: which estimator to use. Defaults to the last
  method fit (e.g. `"onestep"` when present, else the most refined
  correction available).

- ...:

  Unused.

## Value

A \\p \times p\\ sandwich variance-covariance matrix.
