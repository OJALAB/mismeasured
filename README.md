
# mismeasured

<!-- badges: start -->

[![R-CMD-check](https://github.com/OJALAB/mismeasured/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/OJALAB/mismeasured/actions/workflows/R-CMD-check.yaml)
[![Codecov](https://codecov.io/gh/OJALAB/mismeasured/graph/badge.svg)](https://codecov.io/gh/OJALAB/mismeasured)
<!-- badges: end -->

Bias correction for generalized linear models with measurement error and
misclassification.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("OJALAB/mismeasured")
```

## Overview

The `mismeasured` package provides two complementary approaches for
correcting bias in GLMs when covariates are measured with error or
subject to misclassification:

**`simex()`** — Simulation-Extrapolation (SIMEX / MC-SIMEX) with a
formula interface:

- `me()` terms for continuous measurement error, `mc()` terms for
  discrete misclassification
- C++ simulation engine via Rcpp/RcppEigen (100–300x faster than pure-R
  implementations)
- Standard and improved MC-SIMEX (Sevilimedu & Yu, 2026) with exact
  fixed-matrix correction

**`mcglm()`** — Analytical bias correction for GLMs with misclassified
covariates (Battaglia, Christensen, Hansen & Sacher, 2025):

- **Naive** — uncorrected GLM on the proxy covariate
- **BCA** — additive bias correction
- **BCM** — multiplicative bias correction (iterated BCM converges to
  CS)
- **CS** — corrected-score estimator
- **One-step** — joint mixture-likelihood via automatic differentiation
  (RTMB)
- Supports binary and multicategory misclassified covariates,
  Poisson/Binomial/Gaussian families, and multinomial response models

## Quick start

### Continuous measurement error (SIMEX)

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
#> (Intercept)  0.98947    0.02423   40.83   <2e-16 ***
#> x            1.90933    0.02367   80.66   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Misclassification (MC-SIMEX)

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
#> 1            0.88919    0.05082   17.50   <2e-16 ***
#> (Intercept)  0.46438    0.02977   15.60   <2e-16 ***
#> x            0.30831    0.01416   21.77   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Bias-corrected GLM (mcglm)

For analytical bias corrections when misclassification rates are known:

``` r
set.seed(42)
n <- 5000
x <- cbind(1, rnorm(n))
z <- rbinom(n, 1, 0.4)
eta <- 0.8 * z - 0.5 * x[, 1] + 0.7 * x[, 2]
y <- rpois(n, exp(eta))

# Introduce misclassification
p01 <- 0.10; p10 <- 0.15
z_hat <- z
z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01)
z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)

fit <- mcglm(y, z_hat, x, family = "poisson",
             method = c("naive", "bca", "bcm", "cs"),
             p01 = p01, p10 = p10, pi_z = 0.4)
fit
#> Bias-corrected GLM with misclassified covariate
#>   n = 5000, p = 3, K = 2
#> 
#>          NAIVE     BCA     BCM      CS
#> gamma   0.5916  0.7392  0.7884  0.7909
#> alpha0 -0.4003 -0.4806 -0.5074 -0.5136
#> alpha1  0.7076  0.7069  0.7067  0.7067
```

## Formula syntax (simex)

| Term | Meaning | Example |
|----|----|----|
| `me(x, 0.5)` | `x` measured with Gaussian error, sd = 0.5 | `y ~ me(x, 0.5) + w` |
| `me(x, sd_x)` | Heteroscedastic error, sd from column `sd_x` | `y ~ me(x, sd_x) + w` |
| `mc(z, Pi)` | `z` is a misclassified factor, Pi is the K x K misclassification matrix | `y ~ mc(z, Pi) + x` |

## References

- Battaglia, L., Christensen, T., Hansen, S. and Sacher, S. (2025).
  Inference for regression with variables generated by AI or machine
  learning. *arXiv preprint arXiv:2402.15585*.
- Cook, J.R. and Stefanski, L.A. (1994). Simulation-extrapolation
  estimation in parametric measurement error models. *JASA*, 89,
  1314–1328.
- Kuechenhoff, H., Mwalili, S.M. and Lesaffre, E. (2006). A general
  method for dealing with misclassification in regression: The
  misclassification SIMEX. *Biometrics*, 62(1), 85–96.
- Sevilimedu, V. and Yu, L. (2026). An improved misclassification
  simulation extrapolation (MC-SIMEX) algorithm. *Statistics in
  Medicine*, 45, e70418.
