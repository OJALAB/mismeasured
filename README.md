
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
- Asymptotic inference (sandwich SE, Wald CIs) and the usual glm-style
  S3 methods (`summary`, `vcov`, `confint`, `fitted`, `predict`,
  `residuals`, `logLik`, `AIC`, …) are provided for every method

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
z_true <- rbinom(n, 1, 0.4)
x1 <- rnorm(n)
y <- rpois(n, exp(0.8 * z_true - 0.5 + 0.7 * x1))

# Introduce misclassification
p01 <- 0.10; p10 <- 0.15
z_hat <- z_true
z_hat[z_true == 0] <- rbinom(sum(z_true == 0), 1, p01)
z_hat[z_true == 1] <- 1 - rbinom(sum(z_true == 1), 1, p10)

Pi <- matrix(c(1 - p01, p01, p10, 1 - p10), 2, 2)
df3 <- data.frame(y = y, z = factor(z_hat), x1 = x1)

fit <- mcglm(y ~ mc(z, Pi) + x1, data = df3, family = "poisson",
             method = c("naive", "bca", "bcm", "cs"), pi_z = 0.4)
fit
#> 
#> Call:
#> mcglm(formula = y ~ mc(z, Pi) + x1, data = df3, family = "poisson", 
#>     method = c("naive", "bca", "bcm", "cs"), pi_z = 0.4)
#> 
#> Family: poisson  |  n = 5000, K = 2, p = 3
#> Methods: naive, bca, bcm, cs
#> 
#> Coefficients:
#>         NAIVE    BCA      BCM      CS     
#> gamma    0.6272   0.7843   0.8367   0.8400
#> alpha0  -0.4113  -0.4988  -0.5280  -0.5353
#> alpha1   0.7134   0.7136   0.7137   0.7137
#> 
#> Degrees of Freedom: 5000 Total (i.e. Null);  4997 Residual
#> Null Deviance:     9252 
#> Residual Deviance: 5660  | AIC (naive): 12590
```

#### Inference and glm-style methods

`mcglm()` returns asymptotic standard errors for every fitted method
(sandwich estimators from Battaglia et al., 2025). All the usual GLM S3
methods are available; pass `method =` to select an estimator.

``` r
# Wald table per method (estimate, SE, z, p)
summary(fit)
#> 
#> Call:
#> mcglm(formula = y ~ mc(z, Pi) + x1, data = df3, family = "poisson", 
#>     method = c("naive", "bca", "bcm", "cs"), pi_z = 0.4)
#> 
#> Family: poisson  |  n = 5000, K = 2, p = 3
#> Methods: naive, bca, bcm, cs
#> 
#> --- NAIVE ---
#>        Estimate Std. Error z value Pr(>|z|)    
#> gamma   0.62722    0.02866   21.89   <2e-16 ***
#> alpha0 -0.41131    0.02356  -17.46   <2e-16 ***
#> alpha1  0.71335    0.01470   48.52   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> --- BCA ---
#>        Estimate Std. Error z value Pr(>|z|)    
#> gamma   0.78434    0.02951   26.58   <2e-16 ***
#> alpha0 -0.49881    0.02539  -19.64   <2e-16 ***
#> alpha1  0.71364    0.01461   48.85   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> --- BCM ---
#>        Estimate Std. Error z value Pr(>|z|)    
#> gamma   0.83669    0.02994   27.95   <2e-16 ***
#> alpha0 -0.52797    0.02608  -20.24   <2e-16 ***
#> alpha1  0.71374    0.01461   48.87   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> --- CS ---
#>        Estimate Std. Error z value Pr(>|z|)    
#> gamma   0.83996    0.03956   21.23   <2e-16 ***
#> alpha0 -0.53531    0.02980  -17.96   <2e-16 ***
#> alpha1  0.71374    0.01468   48.61   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual deviance (naive): 5660  on 4997 degrees of freedom
#> AIC (naive): 12590
#> 
#> Bias correction (difference from naive):
#>         bca         bcm         cs        
#> gamma    0.1571192   0.2094717   0.2127385
#> alpha0  -0.0874993  -0.1166543  -0.1240022
#> alpha1   0.0002896   0.0003861   0.0003856

# Variance-covariance matrix and confidence intervals for the CS estimator
vcov(fit, method = "cs")
#>                gamma        alpha0        alpha1
#> gamma   1.564896e-03 -0.0009691680  7.824985e-06
#> alpha0 -9.691680e-04  0.0008882385 -1.399587e-04
#> alpha1  7.824985e-06 -0.0001399587  2.155703e-04
confint(fit, method = "cs", level = 0.95)
#>             2.5 %     97.5 %
#> gamma   0.7624216  0.9174891
#> alpha0 -0.5937283 -0.4769014
#> alpha1  0.6849617  0.7425153

# Standard glm helpers, dispatched per-method
coef(fit, method = "bca")
#>      gamma     alpha0     alpha1 
#>  0.7843361 -0.4988120  0.7136425
head(fitted(fit, method = "cs"))
#>         1         2         3         4         5         6 
#> 2.1069839 0.5835903 0.5485658 1.8041813 0.8913723 0.5784742
head(residuals(fit, method = "cs", type = "pearson"))
#>           1           2           3           4           5           6 
#> -0.07370344  1.85410728 -0.74065227 -0.59870637 -0.94412515 -0.76057489
AIC(fit)         # naive log-likelihood when no onestep was fit
#> [1] 12592.08
nobs(fit); family(fit)$family
#> [1] 5000
#> [1] "poisson"
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
