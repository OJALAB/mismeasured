
# mismeasured

High-performance SIMEX and MC-SIMEX correction for measurement error and
misclassification in generalized linear models.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("OJALAB/mismeasured")
```

## Overview

The `mismeasured` package implements the Simulation-Extrapolation
(SIMEX) and Misclassification SIMEX (MC-SIMEX) algorithms for correcting
bias due to measurement error or misclassification in regression models.
It provides:

- A **unified `simex()` function** with brms-style formula terms: `me()`
  for continuous measurement error and `mc()` for discrete
  misclassification
- **C++ simulation engine** via Rcpp/RcppEigen (100-300x faster than
  pure-R implementations)
- **Standard and improved MC-SIMEX** (Sevilimedu & Yu, 2026) with exact
  fixed-matrix correction
- Standard S3 methods: `summary()`, `plot()`, `predict()`, `confint()`,
  `refit()`, etc.

## Quick start

### Continuous measurement error

When a covariate is measured with additive Gaussian error, wrap it with
`me(variable, sd)`:

``` r
library(mismeasured)

set.seed(42)
n <- 2000
x_true <- rnorm(n)
y <- 1 + 2 * x_true + rnorm(n, sd = 0.5)
x_obs <- x_true + rnorm(n, sd = 0.5)  # observed with error
df <- data.frame(y = y, x = x_obs)

fit <- simex(y ~ me(x, 0.5), data = df, B = 200)
summary(fit)
#> 
#> Call:
#> simex(formula = y ~ me(x, 0.5), data = df, B = 200)
#> 
#> Family: gaussian 
#> SIMEX variable(s): x 
#> Extrapolation: quadratic 
#> Lambda grid: 0, 0.5, 1, 1.5, 2 
#> B = 200 , n = 2000 
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -3.38007 -0.78027 -0.01912  0.79010  3.25558 
#> 
#> Naive coefficients:
#> (Intercept)           x 
#>      0.9843      1.5779 
#> 
#> SIMEX corrected coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  0.98947    0.03033   32.62   <2e-16 ***
#> x            1.90933    0.02703   70.63   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Misclassification

When a discrete covariate is subject to misclassification, wrap it with
`mc(variable, matrix)`:

``` r
z_true <- rbinom(n, 1, 0.4)
y2 <- rpois(n, exp(0.5 + 0.8 * z_true + 0.3 * x_true))

# Misclassify z
z_star <- z_true
z_star[z_true == 0] <- rbinom(sum(z_true == 0), 1, 0.10)
z_star[z_true == 1] <- 1 - rbinom(sum(z_true == 1), 1, 0.15)

Pi <- matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)
df2 <- data.frame(y = y2, z = factor(z_star), x = x_true)

# Improved MC-SIMEX (default) -- only needs B=1
fit_mc <- simex(y ~ mc(z, Pi) + x, family = poisson(), data = df2)
summary(fit_mc)
#> 
#> Call:
#> simex(formula = y ~ mc(z, Pi) + x, family = poisson(), data = df2)
#> 
#> Family: poisson 
#> MC-SIMEX variable: z 
#> Method: improved 
#> Extrapolation: exact (improved) 
#> Lambda grid: 0, 1 
#> B = 1 , n = 2000 
#> Estimated P(X=1): 0.4187 
#> Correction factor(s): 1.7676 
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -8.3727 -1.2950 -0.3292  0.8381  8.2882 
#> 
#> Naive coefficients:
#>           1 (Intercept)           x 
#>      0.6257      0.5556      0.3022 
#> 
#> MC-SIMEX corrected coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> 1            0.88919    0.05797  15.340   <2e-16 ***
#> (Intercept)  0.46438    0.04795   9.686   <2e-16 ***
#> x            0.30831    0.02953  10.440   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Use `method = "standard"` for the classic extrapolation-based MC-SIMEX:

``` r
fit_std <- simex(y ~ mc(z, Pi) + x, family = poisson(), data = df2,
                 method = "standard", B = 200)
summary(fit_std)
#> 
#> Call:
#> simex(formula = y ~ mc(z, Pi) + x, family = poisson(), data = df2, 
#>     method = "standard", B = 200)
#> 
#> Family: poisson 
#> MC-SIMEX variable: z 
#> Method: standard 
#> Extrapolation: quadratic 
#> Lambda grid: 0, 0.5, 1, 1.5, 2 
#> B = 200 , n = 2000 
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -6.9372 -1.1856 -0.1909  0.9679  8.4095 
#> 
#> Naive coefficients:
#>           1 (Intercept)           x 
#>      0.6257      0.5556      0.3022 
#> 
#> MC-SIMEX corrected coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> 1            0.82262    0.04271   19.26   <2e-16 ***
#> (Intercept)  0.44074    0.03189   13.82   <2e-16 ***
#> x            0.29552    0.01787   16.53   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## Formula syntax

| Term | Meaning | Example |
|----|----|----|
| `me(x, 0.5)` | `x` measured with Gaussian error, sd = 0.5 | `y ~ me(x, 0.5) + w` |
| `me(x, sd_x)` | Heteroscedastic error, sd from column `sd_x` | `y ~ me(x, sd_x) + w` |
| `mc(z, Pi)` | `z` is a misclassified factor, Pi is the K x K misclassification matrix | `y ~ mc(z, Pi) + x` |

## References

- Cook, J.R. and Stefanski, L.A. (1994). Simulation-extrapolation
  estimation in parametric measurement error models. *JASA*, 89,
  1314–1328.
- Kuchenhoff, H., Mwalili, S.M. and Lesaffre, E. (2006). A general
  method for dealing with misclassification in regression: The
  misclassification SIMEX. *Biometrics*, 62(1), 85–96.
- Sevilimedu, V. and Yu, L. (2026). An improved misclassification
  simulation extrapolation (MC-SIMEX) algorithm. *Statistics in
  Medicine*, 45, e70418.
