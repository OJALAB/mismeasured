# Construct a block diagonal matrix

Construct a block diagonal matrix

## Usage

``` r
diag.block(d, n)
```

## Arguments

- d:

  a list of matrices/vectors, or a single matrix/vector.

- n:

  number of repetitions (used when `d` is not a list).

## Value

A block diagonal matrix.

## Examples

``` r
a <- matrix(1, 2, 2)
b <- matrix(2, 2, 3)
diag.block(list(a, b))
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1    1    0    0    0
#> [2,]    1    1    0    0    0
#> [3,]    0    0    2    2    2
#> [4,]    0    0    2    2    2
```
