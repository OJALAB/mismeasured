
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
- Standard and improved MC-SIMEX (Sevilimedu & Yu, 2026) with
  closed-form fixed-matrix correction (exact for identity-link linear
  models)

**`mcglm()`** — Analytical bias correction for GLMs with misclassified
covariates (Yi, Yan, Liao & Spiegelman, 2019; Battaglia, Christensen,
Hansen & Sacher, 2025):

- **SUB** (default) — subtraction correction; **EC** — expectation
  correction; **IL** — induced likelihood
- **CS-AKN** — corrected score of Akazawa, Kinukawa & Nakamura (1998);
  needs only the misclassification matrix and allows the true category
  to depend on the covariates
- **CS** — drift-corrected score; **BCA** / **BCM** — additive and
  multiplicative bias corrections
- **One-step** — joint mixture-likelihood via automatic differentiation
  (RTMB)
- **Validation samples**: misclassification probabilities estimated from
  internal or external audits (`validation_sample()`, `estimate_mc()`),
  with design weights, strata, a covariate-dependent prevalence,
  shrinkage for sparse audits, and their uncertainty propagated into the
  standard errors
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
#> 1            0.82262    0.04275   19.24   <2e-16 ***
#> (Intercept)  0.44074    0.03191   13.81   <2e-16 ***
#> x            0.29552    0.01788   16.53   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Bias-corrected GLM (mcglm)

For analytical bias corrections when misclassification rates are known
(see `vignette("validation", "mismeasured")` when they are estimated
from an audit):

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
             method = c("naive", "sub", "il", "cs_akn"), pi_z = 0.4)
fit
#> 
#> Call:
#> mcglm(formula = y ~ mc(z, Pi) + x1, data = df3, family = "poisson", 
#>     method = c("naive", "sub", "il", "cs_akn"), pi_z = 0.4)
#> 
#> Family: poisson  |  n = 5000, K = 2, p = 3
#> Methods: naive, sub, il, cs_akn
#> 
#> Coefficients:
#>              NAIVE    CS_AKN   SUB      IL     
#> gamma         0.6272   0.8401   0.8400   0.8372
#> (Intercept)  -0.4113  -0.5358  -0.5350  -0.5341
#> x1            0.7134   0.7143   0.7134   0.7126
#> 
#> Degrees of Freedom: 5000 Total (i.e. Null);  4997 Residual
#> Null Deviance:     9252 
#> Residual Deviance: 5660  | AIC (naive): 12590
```

#### Inference and glm-style methods

`mcglm()` returns asymptotic standard errors for every fitted method.
All the usual GLM S3 methods are available; pass `method =` to select an
estimator.

``` r
# Wald table per method (estimate, SE, z, p)
summary(fit)
#> 
#> Call:
#> mcglm(formula = y ~ mc(z, Pi) + x1, data = df3, family = "poisson", 
#>     method = c("naive", "sub", "il", "cs_akn"), pi_z = 0.4)
#> 
#> Family: poisson  |  n = 5000, K = 2, p = 3
#> Methods: naive, sub, il, cs_akn
#> z categories (Pi assumed in this order): 0 (baseline), 1
#> 
#> --- NAIVE ---
#>             Estimate Std. Error z value Pr(>|z|)    
#> gamma        0.62722    0.02866   21.89   <2e-16 ***
#> (Intercept) -0.41131    0.02356  -17.46   <2e-16 ***
#> x1           0.71335    0.01470   48.52   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> --- SUB ---
#>             Estimate Std. Error z value Pr(>|z|)    
#> gamma        0.83997    0.03955   21.24   <2e-16 ***
#> (Intercept) -0.53497    0.02977  -17.97   <2e-16 ***
#> x1           0.71335    0.01470   48.52   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> --- IL ---
#>             Estimate Std. Error z value Pr(>|z|)    
#> gamma        0.83720    0.03434   24.38   <2e-16 ***
#> (Intercept) -0.53413    0.02720  -19.64   <2e-16 ***
#> x1           0.71256    0.01489   47.86   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> --- CS_AKN ---
#>             Estimate Std. Error z value Pr(>|z|)    
#> gamma        0.84005    0.03960   21.21   <2e-16 ***
#> (Intercept) -0.53584    0.03012  -17.79   <2e-16 ***
#> x1           0.71425    0.01520   47.00   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual deviance (naive): 5660  on 4997 degrees of freedom
#> AIC (naive): 12590
#> 
#> Bias correction (difference from naive):
#>              cs_akn      sub         il        
#> gamma         2.128e-01   2.128e-01   2.100e-01
#> (Intercept)  -1.245e-01  -1.237e-01  -1.228e-01
#> x1            8.984e-04   2.226e-10  -7.883e-04

# Variance-covariance matrix and confidence intervals for SUB
vcov(fit, method = "sub")
#>                     gamma   (Intercept)            x1
#> gamma        1.564521e-03 -0.0009690157  7.023513e-06
#> (Intercept) -9.690157e-04  0.0008861499 -1.404931e-04
#> x1           7.023513e-06 -0.0001404931  2.161381e-04
confint(fit, method = "sub", level = 0.95)
#>                  2.5 %     97.5 %
#> gamma        0.7624468  0.9174957
#> (Intercept) -0.5933166 -0.4766271
#> x1           0.6845382  0.7421676

# Standard glm helpers, dispatched per-method
coef(fit, method = "il")
#>       gamma (Intercept)          x1 
#>   0.8372042  -0.5341347   0.7125646
head(fitted(fit, method = "sub"))
#>         1         2         3         4         5         6 
#> 2.1072387 0.5837916 0.5487733 1.8045507 0.8914757 0.5786764
head(residuals(fit, method = "sub", type = "pearson"))
#>           1           2           3           4           5           6 
#> -0.07387449  1.85352421 -0.74079236 -0.59892008 -0.94417991 -0.76070784
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
