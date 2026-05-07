# Generate misclassified data

Takes a data frame of factor variables and produces misclassified
versions using the misclassification matrix raised to power `k` via
eigendecomposition.

## Usage

``` r
misclass(data.org, mc.matrix, k = 1)
```

## Arguments

- data.org:

  data frame containing factor variables.

- mc.matrix:

  a named list of misclassification matrices. Names must correspond to
  variable names in `data.org`. Column names must match factor levels.

- k:

  exponent for the misclassification matrix (default: 1).

## Value

A data frame with the misclassified variables.

## Examples

``` r
x1 <- factor(rbinom(100, 1, 0.5))
p1 <- matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2)
colnames(p1) <- levels(x1)
x <- data.frame(x1 = x1)
x.mc <- misclass(x, list(x1 = p1), k = 1)
```
