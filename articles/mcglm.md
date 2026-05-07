# Bias-corrected GLMs with misclassified covariates: mcglm() and mclm()

## Scope

This vignette describes the analytical bias-correction estimators in
**mismeasured**:
[`mcglm()`](https://ojalab.github.io/mismeasured/reference/mcglm.md) for
any supported family and
[`mclm()`](https://ojalab.github.io/mismeasured/reference/mclm.md), the
Gaussian-fixed wrapper of
[`mcglm()`](https://ojalab.github.io/mismeasured/reference/mcglm.md).
Five estimators are implemented for GLMs in which a discrete covariate
is observed only through a misclassified proxy: naive, additive bias
correction (BCA), multiplicative bias correction (BCM), corrected score
(CS, two formulations), and a one-step joint mixture-likelihood
estimator.

For simulation-extrapolation methods, see the companion vignette
[`vignette("simex", "mismeasured")`](https://ojalab.github.io/mismeasured/articles/simex.md).

## Setup and notation

Let $`Y_i`$ be a scalar response, $`x_i \in \mathbb{R}^r`$ a vector of
cleanly observed covariates (typically including an intercept), and
$`Z_i \in \{0, 1, \dots, K-1\}`$ a latent categorical regressor observed
only through a proxy $`\hat Z_i`$. For binary problems ($`K = 2`$) take
$`\xi_i = (Z_i, x_i^\top)^\top`$ and
$`\hat\xi_i = (\hat Z_i, x_i^\top)^\top`$; for $`K \ge 3`$, $`Z_i`$ is
dummy-encoded with baseline $`0`$ so
$`\xi_i = (d_i^\top, x_i^\top)^\top`$. The parameter is
$`\psi = (\gamma^\top, \alpha^\top)^\top`$ with
$`\gamma \in \mathbb{R}^{K-1}`$ on the (dummy-coded) latent regressor
and $`\alpha \in \mathbb{R}^r`$ on $`x_i`$.

The conditional mean is
$`\mathbb{E}[Y_i \mid Z_i, x_i] = \mu(\psi^\top \xi_i)`$ for
$`\mu = g^{-1}`$ the inverse canonical link.

**Misclassification mechanism.** For $`K = 2`$ this is summarized by
$`p_{01} = \Pr(\hat Z = 1 \mid Z = 0)`$,
$`p_{10} = \Pr(\hat Z = 0 \mid Z = 1)`$, and prevalence
$`\pi = \Pr(Z = 1)`$. For $`K \ge 3`$, by the $`K \times K`$
column-stochastic matrix $`\Pi_{j\ell} = \Pr(\hat Z = j \mid Z = \ell)`$
and the prevalence vector $`(\pi_0, \dots, \pi_{K-1})`$. We assume
nondifferential misclassification: \$\hat Z_i
\\\rule\[1pt\]{8pt}{0.4pt}\\\\\\\rule\[-3pt\]{0.4pt}{8pt}\rule\[1pt\]{8pt}{0.4pt}\\\\\\\rule\[-3pt\]{0.4pt}{8pt}\\
(Y_i, x_i) \mid Z_i\$.

## The naive proxy-score estimator and its bias

The naive estimator solves the proxy-score equation
``` math
\hat U_n(\psi) = \frac{1}{n} \sum_{i=1}^n \hat\xi_i \{Y_i - \mu(\psi^\top \hat\xi_i)\} = 0.
```
This is just GLM estimation with $`\hat Z_i`$ in place of $`Z_i`$. It is
typically inconsistent for $`\psi_0`$. For binary $`Z`$, the package
theory (Warmuz & Beresewicz, 2026, Proposition on score drift) shows
``` math
\mathbb{E}[\hat U_n(\psi_0)] = \mathbb{E}[m_i(\psi_0)],
```
where the per-observation drift is
``` math
m_i(\psi) = \begin{pmatrix} -c_1\, \delta_i(\psi) \\ -c_2\, \delta_i(\psi)\, x_i \end{pmatrix},
\qquad
\delta_i(\psi) = \mu(\gamma + \alpha^\top x_i) - \mu(\alpha^\top x_i),
```
with constants $`c_1 = p_{01}(1-\pi)`$ and
$`c_2 = p_{01}(1-\pi) - p_{10}\pi`$. The multicategory version replaces
$`(\delta_i, c_1, c_2)`$ with explicit double sums over $`(\Pi, \pi_z)`$
(Appendix B of the package theory note).

## The five estimators

| Method | Idea | Required nuisance | Asymptotic regime |
|----|----|----|----|
| `naive` | proxy GLM, no correction | – | inconsistent under fixed $`\Pi`$ |
| `bca` | additive: $`\hat\psi_{\mathrm{bca}} = \hat\psi - \hat I(\hat\psi)^{-1}\hat m(\hat\psi)`$ | $`(c_1, c_2)`$ or $`(\Pi, \pi)`$ | drifting $`\Pi = O(n^{-1/2})`$ |
| `bcm` | multiplicative: $`\hat\psi_{\mathrm{bcm}} = \hat\psi - (\hat I + \hat M)^{-1} \hat m`$ | same | drifting $`\Pi = O(n^{-1/2})`$ |
| `cs` | corrected score: solve $`\sum_i [\hat\xi_i\{Y_i - \mu(\psi^\top\hat\xi_i)\} - m_i(\psi)] = 0`$ | same | **fixed $`\Pi`$** |
| `cs_akn` | AKN98 corrected score via unbiased surrogate $`x = Q^{-1}(u - p_0)`$ | $`\Pi`$ only | **fixed $`\Pi`$** |
| `onestep` | maximize $`\prod_i \sum_\ell \pi_\ell \Pi_{\hat z_i+1, \ell+1} f(Y_i \mid Z = \ell, \psi)`$ via RTMB | (see below) | fixed $`\Pi`$ |

`bca` and `bcm` remove the leading $`n^{-1/2}`$ bias under a
drifting-sequence regime where misclassification probabilities shrink
with $`n`$. `cs` and `cs_akn` are unbiased estimating equations under
fixed $`\Pi`$ – they remain consistent even when misclassification is
moderate or large. `onestep` is a likelihood-based alternative
implemented through automatic differentiation (`RTMB`); it estimates
mixture weights jointly with $`\psi`$ unless `fix_omega = TRUE`.

### Two CS formulations

The package implements *two* corrected-score estimators that solve the
same population equation but differ in finite-sample form:

- **`cs` (drift)** keeps the proxy as-is and adds the explicit drift
  $`m_i(\psi)`$ to the naive score. Needs $`(\Pi, \pi_z)`$. (Warmuz &
  Beresewicz, 2026, Theorem on the corrected-score estimator; Battaglia
  et al., 2025 for the linear-model origin.)
- **`cs_akn` (unbiased surrogate)** transforms the regressor via
  $`x = Q^{-1}(u - p_0)`$ where
  $`Q = [p_1 - p_0, \ldots, p_{K-1} - p_0]`$, then substitutes into the
  log-likelihood. Needs only $`\Pi`$ – $`\pi_z`$ cancels out. (Akazawa,
  Kinukawa & Nakamura, 1998.)

The two are population-equivalent at $`\psi_0`$. Empirically they agree
to $`\sim 10^{-4}`$ for $`K = 2`$ and within sampling noise for
$`K \ge 3`$; `cs_akn`’s SD is slightly larger when $`\Pi`$ has
off-diagonals away from zero (the $`Q^{-1}`$ amplifies
single-observation contributions).

### Inference

All five estimators report an asymptotic variance through
[`vcov()`](https://rdrr.io/r/stats/vcov.html):

- naive / BCA / BCM: $`A^{-1} C A^{-1}`$ where
  $`A = \mathbb{E}[\dot\mu(\psi_0^\top\xi_i)\xi_i\xi_i^\top]`$ and
  $`C = \mathbb{E}[\varepsilon_i^2 \xi_i \xi_i^\top]`$.
- CS / cs_akn: $`J^{-1} S J^{-\top}`$ with $`J = -(\hat I + \hat M)`$
  (cs) or $`J = -I_{\mathrm{akn}}`$ (cs_akn) and
  $`S = \mathbb{E}[\phi_i \phi_i^\top]`$.
- onestep: inverse Hessian of the integrated log-likelihood.

Set `vcov_corrected = TRUE` in
[`mcglm()`](https://ojalab.github.io/mismeasured/reference/mcglm.md) to
get a more conservative BCA/BCM sandwich that accounts for the
additional uncertainty from estimating $`\hat m(\hat\psi)`$.

## Worked examples

### Binary misclassification, Poisson response

``` r

n  <- 5000
Pi <- matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2)   # p01 = p10 = 0.10
z  <- rbinom(n, 1, 0.4)
z_hat <- z
z_hat[z == 0] <- rbinom(sum(z == 0), 1, 0.10)
z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.10)
x1 <- rnorm(n)
y  <- rpois(n, exp(0.8 * z - 0.5 + 0.7 * x1))

df <- data.frame(y = y, z = factor(z_hat), x1 = x1)

fit <- mcglm(y ~ mc(z, Pi) + x1, data = df, family = "poisson",
             method = c("naive", "bca", "bcm", "cs", "cs_akn"))
fit
#> 
#> Call:
#> mcglm(formula = y ~ mc(z, Pi) + x1, data = df, family = "poisson", 
#>     method = c("naive", "bca", "bcm", "cs", "cs_akn"))
#> 
#> Family: poisson  |  n = 5000, K = 2, p = 3
#> Methods: naive, bca, bcm, cs, cs_akn
#> 
#> Coefficients:
#>         NAIVE    BCA      BCM      CS       CS_AKN 
#> gamma    0.6803   0.8182   0.8548   0.8554   0.8559
#> alpha0  -0.4362  -0.5043  -0.5223  -0.5265  -0.5253
#> alpha1   0.6906   0.6885   0.6879   0.6880   0.6849
#> 
#> Degrees of Freedom: 5000 Total (i.e. Null);  4997 Residual
#> Null Deviance:     9333 
#> Residual Deviance: 5605  | AIC (naive): 12610
```

The naive $`\hat\gamma`$ is biased toward zero; all four corrections
recover $`\gamma_0 = 0.8`$.

``` r

summary(fit)
#> 
#> Call:
#> mcglm(formula = y ~ mc(z, Pi) + x1, data = df, family = "poisson", 
#>     method = c("naive", "bca", "bcm", "cs", "cs_akn"))
#> 
#> Family: poisson  |  n = 5000, K = 2, p = 3
#> Methods: naive, bca, bcm, cs, cs_akn
#> 
#> --- NAIVE ---
#>        Estimate Std. Error z value Pr(>|z|)    
#> gamma   0.68035    0.02846   23.90   <2e-16 ***
#> alpha0 -0.43620    0.02370  -18.41   <2e-16 ***
#> alpha1  0.69055    0.01553   44.46   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> --- BCA ---
#>        Estimate Std. Error z value Pr(>|z|)    
#> gamma   0.81819    0.02908   28.14   <2e-16 ***
#> alpha0 -0.50429    0.02517  -20.04   <2e-16 ***
#> alpha1  0.68848    0.01568   43.90   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> --- BCM ---
#>        Estimate Std. Error z value Pr(>|z|)    
#> gamma   0.85477    0.02932   29.15   <2e-16 ***
#> alpha0 -0.52234    0.02560  -20.41   <2e-16 ***
#> alpha1  0.68790    0.01574   43.70   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> --- CS ---
#>        Estimate Std. Error z value Pr(>|z|)    
#> gamma   0.85535    0.03635   23.53   <2e-16 ***
#> alpha0 -0.52649    0.02790  -18.87   <2e-16 ***
#> alpha1  0.68795    0.01594   43.16   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> --- CS_AKN ---
#>        Estimate Std. Error z value Pr(>|z|)    
#> gamma   0.85593    0.03639   23.52   <2e-16 ***
#> alpha0 -0.52533    0.02806  -18.72   <2e-16 ***
#> alpha1  0.68490    0.01669   41.04   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual deviance (naive): 5605  on 4997 degrees of freedom
#> AIC (naive): 12610
#> 
#> Bias correction (difference from naive):
#>         bca        bcm        cs         cs_akn   
#> gamma    0.137840   0.174425   0.175005   0.175581
#> alpha0  -0.068093  -0.086142  -0.090290  -0.089125
#> alpha1  -0.002073  -0.002656  -0.002603  -0.005655
```

### Binary misclassification, logistic response

``` r

n  <- 5000
Pi <- matrix(c(0.85, 0.15, 0.10, 0.90), 2, 2)
z  <- rbinom(n, 1, 0.4); zh <- z
zh[z == 0] <- rbinom(sum(z == 0), 1, 0.15)
zh[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.10)
x1 <- rnorm(n)
y  <- rbinom(n, 1, plogis(0.9 * z - 0.4 + 0.6 * x1))

dfb <- data.frame(y = y, z = factor(zh), x1 = x1)
fitb <- mcglm(y ~ mc(z, Pi) + x1, data = dfb, family = "binomial",
              method = c("naive", "cs", "cs_akn"))
coef(fitb, method = "naive")
#>      gamma     alpha0     alpha1 
#>  0.6661615 -0.3362542  0.5299033
coef(fitb, method = "cs")
#>      gamma     alpha0     alpha1 
#>  0.9279532 -0.4041620  0.5409704
coef(fitb, method = "cs_akn")
#>      gamma     alpha0     alpha1 
#>  0.9276973 -0.4041227  0.5399984
```

### Multicategory misclassification ($`K = 3`$)

For $`K \ge 3`$ pass the misclassification matrix and (for
`cs`/`bca`/`bcm`) the prevalence vector $`\pi`$. `cs_akn` needs only
$`\Pi`$.

``` r

n   <- 5000
Pi3 <- matrix(0.1, 3, 3); diag(Pi3) <- 0.8
pi3 <- c(0.5, 0.3, 0.2)
Z  <- sample(0:2, n, replace = TRUE, prob = pi3)
Zh <- vapply(Z, function(z) sample(0:2, 1, prob = Pi3[, z + 1]), 0L)
x1 <- rnorm(n)
ym <- rpois(n, exp(c(0, 1.0, -0.9)[Z + 1] + 0.8 - 0.7 * x1))

fit3 <- mcglm(ym, z_hat = Zh, x = cbind(1, x1), family = "poisson",
              method = c("naive", "cs", "cs_akn"),
              Pi = Pi3, pi_z = pi3)
coef(fit3, method = "cs")
#>     gamma1     gamma2     alpha0     alpha1 
#>  0.9974917 -1.2197661  0.8089698 -0.6927194
coef(fit3, method = "cs_akn")
#>     gamma1     gamma2     alpha0     alpha1 
#>  0.9985326 -1.2428888  0.8066988 -0.6920223
```

### Confidence intervals

``` r

confint(fit, method = "cs")
#>             2.5 %     97.5 %
#> gamma   0.7841097  0.9265936
#> alpha0 -0.5811827 -0.4717990
#> alpha1  0.6567124  0.7191898
confint(fit, method = "cs_akn")
#>             2.5 %     97.5 %
#> gamma   0.7845970  0.9272574
#> alpha0 -0.5803202 -0.4703317
#> alpha1  0.6521881  0.7176108
```

### One-step (joint mixture likelihood)

The one-step estimator maximizes the integrated mixture likelihood
across the latent categories. It uses `RTMB` for automatic
differentiation and is robust under fixed $`\Pi`$ for any $`K`$. Mixture
weights can be estimated jointly (default) or fixed at user-supplied
$`\pi_z`$ via `fix_omega = TRUE`.

``` r

fit_os <- mcglm(y ~ mc(z, Pi) + x1, data = df, family = "poisson",
                method = c("naive", "onestep"))
coef(fit_os, method = "onestep")
#>      gamma     alpha0     alpha1 
#>  0.8492809 -0.5019135  0.6972155
```

## `mclm()`: Gaussian-fixed convenience wrapper

[`mclm()`](https://ojalab.github.io/mismeasured/reference/mclm.md) is a
thin wrapper around
[`mcglm()`](https://ojalab.github.io/mismeasured/reference/mcglm.md)
with `family = "gaussian"`, in the spirit of base R’s
[`lm()`](https://rdrr.io/r/stats/lm.html) versus
[`glm()`](https://rdrr.io/r/stats/glm.html). All other arguments forward
unchanged and the returned object remains of class `"mcglm"`, so all S3
methods (`summary`, `vcov`, `confint`, …) work identically.

``` r

n  <- 4000
Pi <- matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2)
z  <- rbinom(n, 1, 0.4); zh <- z
zh[z == 0] <- rbinom(sum(z == 0), 1, 0.10)
zh[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.10)
x1 <- rnorm(n)
yg <- 0.8 * z - 0.5 + 0.7 * x1 + rnorm(n, 0, 0.5)

dfg <- data.frame(y = yg, z = factor(zh), x1 = x1)
fitg <- mclm(y ~ mc(z, Pi) + x1, data = dfg,
             method = c("naive", "cs", "cs_akn"))
coef(fitg, method = "cs")
#>      gamma     alpha0     alpha1 
#>  0.7913075 -0.4782274  0.6984542
coef(fitg, method = "cs_akn")
#>      gamma     alpha0     alpha1 
#>  0.7913625 -0.4782466  0.6998920
```

For Gaussian, `cs_akn` admits a closed-form linear-system solution (no
iterative solver), so the two CS formulations agree to numerical
precision.

## When to use which method

- **Reliable $`(\Pi, \pi_z)`$, well-curated proxy.** Any of `cs`,
  `cs_akn`, or `onestep`. Use `cs` for the most established formulation;
  `onestep` when you also want a likelihood-based confidence set or
  model fit comparison.
- **Reliable $`\Pi`$ but $`\pi_z`$ unknown or unstable.** Prefer
  `cs_akn` – the surrogate construction makes $`\pi_z`$ irrelevant.
- **Misclassification probabilities small and shrinking with $`n`$.**
  `bca` or `bcm` are valid (drifting-sequence regime). `bcm` removes
  second-order bias and is closer to the iterated CS update.
- **Multicategory $`K \ge 3`$.** All five estimators apply; `cs_akn`
  carries a small variance penalty (≈ 10–20% larger SD) versus drift
  `cs`.
- **Multinomial response.** Only `naive` and `onestep` (the others are
  univariate-response specific).

## Estimated misclassification parameters

The package theory (Warmuz & Beresewicz, 2026, Proposition on plug-in
robustness) shows that the corrected-score estimator remains consistent
when $`(\Pi, \pi)`$ are estimated from auxiliary data, provided the
plug-in score satisfies
$`\sup_\psi \|\Phi_n^* - \Phi_n^\circ\| = o_P(1)`$, and is
asymptotically unbiased under the tighter rate $`o_P(n^{-1/2})`$. In
practice this means a validation subsample of moderate size (e.g. 5–10%
of the main sample) is sufficient.

## Result object and S3 methods

[`mcglm()`](https://ojalab.github.io/mismeasured/reference/mcglm.md) and
[`mclm()`](https://ojalab.github.io/mismeasured/reference/mclm.md) both
return objects of class `"mcglm"`. The following S3 methods are
available:

- [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) – coefficient
  tables per method.
- [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html), `se()`,
  [`confint()`](https://rdrr.io/r/stats/confint.html) – accept a
  `method =` argument to select which estimator’s quantity to extract.
- [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) – on the
  original (proxy) data.
- [`logLik()`](https://rdrr.io/r/stats/logLik.html),
  [`AIC()`](https://rdrr.io/r/stats/AIC.html),
  [`BIC()`](https://rdrr.io/r/stats/AIC.html),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html) – only meaningful for
  `onestep`.

## References

Akazawa, K., Kinukawa, N. and Nakamura, T. (1998). A note on the
corrected score function adjusting for misclassification. *J. Japan
Statist. Soc.* 28(1), 115–123. (Foundational existence theorem and
explicit corrected scores for misclassified discrete covariates; basis
of `cs_akn`.)

Battaglia, L., Christensen, T., Hansen, S. and Sacher, S. (2025).
Inference for regression with variables generated by AI or machine
learning. *arXiv preprint* arXiv:2402.15585. (Drifting-asymptotics
framework for linear models; origin of BCA/BCM and one-step.)

Nakamura, T. (1990). Corrected score function for errors-in-variables
models: methodology and application to generalized linear models.
*Biometrika* 77(1), 127–137. (Originating reference for the
corrected-score concept under continuous additive error.)

Warmuz, M. and Beresewicz, M. (2026). Bias corrections for GLMs with
misclassified covariates. Working paper. (GLM extension of the
drifting-asymptotics regime; multicategory drift formulas; plug-in
robustness conditions used in the package.)
