# Bias-Corrected GLM Estimation with Misclassified Covariates

Fits a GLM in which one covariate is observed through a misclassified
proxy and applies the bias-correction estimators of Battaglia,
Christensen, Hansen and Sacher (2025). Five estimators are available:
*naive* (uncorrected proxy-score), *BCA* (additive bias correction),
*BCM* (multiplicative bias correction), *CS* (corrected score), and
*one-step* (joint mixture likelihood via RTMB). Both binary and
multicategory latent regressors are supported.

## Usage

``` r
mcglm(
  formula,
  data = NULL,
  family = "poisson",
  method = c("naive", "bca", "bcm", "cs"),
  p01 = NULL,
  p10 = NULL,
  pi_z = NULL,
  Pi = NULL,
  K = NULL,
  c1 = NULL,
  c2 = NULL,
  iterate = FALSE,
  jacobian = c("analytical", "numerical"),
  fix_omega = FALSE,
  vcov_corrected = FALSE,
  weights = NULL,
  J = NULL,
  homoskedastic = TRUE,
  optim_control = list(),
  z_hat = NULL,
  x = NULL
)
```

## Arguments

- formula:

  Either a formula of the form `y ~ mc(z, Pi) + x1 + x2`, where
  [`mc()`](https://ojalab.github.io/mismeasured/reference/mc.md) marks
  the misclassified covariate \\\hat Z_i\\ together with its
  misclassification matrix \\\Pi\\; or, when using the matrix interface,
  a numeric response vector of length \\n\\.

- data:

  A data frame (required for the formula interface). The variable named
  in [`mc()`](https://ojalab.github.io/mismeasured/reference/mc.md)
  should be a [`factor`](https://rdrr.io/r/base/factor.html) or
  integer-valued vector with levels coded as \\0, 1, \dots, K-1\\.

- family:

  A [`family`](https://rdrr.io/r/stats/family.html) object or one of the
  strings `"poisson"`, `"binomial"`, `"gaussian"` (any \\K\\, all five
  methods), or `"multinomial"` (`naive` and `onestep` only).

- method:

  Character vector of estimators to fit. Any subset of
  `c("naive", "bca", "bcm", "cs", "cs_akn", "onestep")`; the default is
  the four analytical estimators `c("naive", "bca", "bcm", "cs")`. The
  `"cs_akn"` entry selects the Akazawa–Kinukawa–Nakamura (1998)
  corrected-score construction (an alternative formulation of `"cs"`
  based on the unbiased-surrogate transform \\x = Q^{-1}(u - p_0)\\); it
  needs only \\\Pi\\ (not \\\pi_z\\) and is unsupported for
  `family = "multinomial"`.

- p01:

  False-positive rate \\p\_{01} = \Pr(\hat Z = 1 \mid Z = 0)\\ (`K = 2`
  only). Auto-extracted from `Pi` when supplied.

- p10:

  False-negative rate \\p\_{10} = \Pr(\hat Z = 0 \mid Z = 1)\\ (`K = 2`
  only). Auto-extracted from `Pi` when supplied.

- pi_z:

  Latent prevalence: a scalar \\\pi = \Pr(Z = 1)\\ (`K = 2`) or a
  \\K\\-vector \\(\pi_0, \dots, \pi\_{K-1})\\ of class probabilities. If
  `NULL`, `pi_z` is estimated from `z_hat` and `Pi` via \\\hat\pi =
  \Pi^{-1} \hat\pi\_{\text{obs}}\\ (clamped to \\\[0.01, 0.99\]\\).

- Pi:

  The \\K \times K\\ misclassification matrix \\\Pi\_{j\ell} = \Pr(\hat
  Z = j - 1 \mid Z = \ell - 1)\\; columns must sum to 1. Extracted
  automatically from the
  [`mc()`](https://ojalab.github.io/mismeasured/reference/mc.md) term
  when using the formula interface.

- K:

  Number of latent categories. Auto-detected from `Pi` or, failing that,
  from the levels of `z_hat`.

- c1:

  Pre-computed misclassification constant \\c_1 = p\_{01}(1 - \pi)\\
  (binary case, advanced use). Computed from `p01` and `pi_z` when
  `NULL`.

- c2:

  Pre-computed misclassification constant \\c_2 = p\_{01}(1 - \pi) -
  p\_{10} \pi\\ (binary case, advanced use). Computed from `p01`, `p10`
  and `pi_z` when `NULL`.

- iterate:

  Logical. If `TRUE`, iterate the BCA/BCM updates to convergence (the
  iterated BCM solves the corrected score equation; iterated BCA
  performs Newton steps on the additive correction).

- jacobian:

  Character. Either `"analytical"` (default) or `"numerical"`. Controls
  how the drift Jacobian \\M = \partial m / \partial \psi\\ is evaluated
  for the multicategory (\\K \ge 3\\) BCM and CS estimators. The
  analytical form is exact and faster; the numerical form (forward
  differences) is provided as a fallback for cross-checking. The binary
  case (\\K = 2\\) always uses the analytical Jacobian.

- fix_omega:

  Logical. If `TRUE`, fix the mixture weights in the one-step estimator
  at the values implied by the supplied misclassification parameters; if
  `FALSE` (default), they are estimated jointly with \\\psi\\ via a
  softmax parameterization.

- vcov_corrected:

  Logical. If `TRUE`, the BCA/BCM variance accounts for the additional
  uncertainty due to plug-in estimation of \\\hat m(\hat\psi)\\ (joint
  score-and-drift sandwich); if `FALSE` (default), the drifting-regime
  asymptotic variance \\A^{-1} C A^{-1}\\ is used. Both share the same
  point estimate.

- weights:

  Optional positive frequency weights of length \\n\\.

- J:

  Number of response categories (multinomial family only;
  auto-detected).

- homoskedastic:

  Logical. For one-step Gaussian fits, assume a single residual variance
  across mixture components; if `FALSE`, binary fits estimate one
  variance per latent class, while \\K \> 2\\ fits estimate one baseline
  variance and one pooled nonbaseline variance.

- optim_control:

  List of control parameters passed to
  [`nlminb`](https://rdrr.io/r/stats/nlminb.html) for the one-step
  estimator.

- z_hat:

  Integer vector of observed proxy values \\\hat Z_i \in \\0, \dots,
  K-1\\\\ (matrix interface only).

- x:

  Covariate matrix \\\[x_1, \dots, x_n\]^\top\\ of dimension \\n \times
  r\\; should include an intercept column when one is desired (matrix
  interface only).

## Value

An object of class `"mcglm"` with components:

- coefficients:

  Named list of coefficient vectors \\\hat\psi\\, one per method.

- vcov:

  Named list of asymptotic variance-covariance matrices, one per method
  (sandwich estimators from the paper's theorems).

- se:

  Named list of standard errors per method.

- naive_fit:

  The [`glm`](https://rdrr.io/r/stats/glm.html) object from the naive
  fit (`NULL` for multinomial).

- method, family, K, n, p:

  Metadata.

- weights:

  Frequency weights used (`NULL` if unweighted).

- loglik_onestep, npar_onestep, vcov_onestep:

  One-step log-likelihood, full optimized parameter count, and variance
  (when `onestep` is fit).

- call, formula:

  The matched call and the formula (`NULL` if the matrix interface was
  used).

See
[`summary.mcglm`](https://ojalab.github.io/mismeasured/reference/summary.mcglm.md),
[`vcov.mcglm`](https://ojalab.github.io/mismeasured/reference/vcov.mcglm.md),
[`confint.mcglm`](https://ojalab.github.io/mismeasured/reference/confint.mcglm.md),
[`predict.mcglm`](https://ojalab.github.io/mismeasured/reference/predict.mcglm.md),
and
[`logLik.mcglm`](https://ojalab.github.io/mismeasured/reference/logLik.mcglm.md)
for downstream methods.

## Model and notation

Let \\Y_i\\ denote a scalar response, \\x_i \in \mathbb{R}^r\\ a vector
of cleanly observed covariates (typically including an intercept), and
\\Z_i \in \\0, 1, \dots, K-1\\\\ a latent categorical regressor that is
observed only through a proxy \\\hat Z_i\\. For the binary case (\\K =
2\\) we let \\\xi_i = (Z_i, x_i^\top)^\top\\ and \\\hat\xi_i = (\hat
Z_i, x_i^\top)^\top\\; for \\K \> 2\\, \\Z_i\\ is dummy-encoded with
baseline category 0 and \\\xi_i = (d_i^\top, x_i^\top)^\top\\. The
parameter is \\\psi = (\gamma^\top, \alpha^\top)^\top\\, where \\\gamma
\in \mathbb{R}^{K-1}\\ are the coefficients on the (dummy-coded) latent
regressor and \\\alpha \in \mathbb{R}^r\\ the coefficients on \\x_i\\.
The conditional mean is \\E\[Y_i \| Z_i, x_i\] = \mu(\psi^\top \xi_i)\\
for \\\mu = g^{-1}\\ the inverse canonical link.

*Misclassification mechanism.* For \\K = 2\\, the mechanism is
summarized by the false-positive rate \\p\_{01} = \Pr(\hat Z = 1 \mid Z
= 0)\\, the false-negative rate \\p\_{10} = \Pr(\hat Z = 0 \mid Z = 1)\\
and the prevalence \\\pi = \Pr(Z = 1)\\. For \\K \> 2\\, by the \\K
\times K\\ matrix \\\Pi\_{j\ell} = \Pr(\hat Z = j - 1 \mid Z = \ell -
1)\\ (column-stochastic) and the prevalence vector \\(\pi_0, \dots,
\pi\_{K-1})\\. Misclassification is assumed nondifferential: \\\hat Z
\perp\\\\\\\perp (Y, x) \mid Z\\.

*Estimators.* Let \\\hat\psi\\ solve the proxy score \\\hat U_n(\psi) =
n^{-1} \sum_i \hat\xi_i \\Y_i - \mu(\psi^\top \hat\xi_i)\\ = 0\\ (this
is the naive GLM estimator). Define the per-observation drift
\\m_i(\psi) = (-c_1 \delta_i(\psi),\\ -c_2 \delta_i(\psi)\\ x_i)^\top\\
with \\\delta_i(\psi) = \mu(\gamma + \alpha^\top x_i) - \mu(\alpha^\top
x_i)\\, \\c_1 = p\_{01}(1 - \pi)\\ and \\c_2 = p\_{01}(1 - \pi) -
p\_{10} \pi\\, and let \\\hat I(\psi) = n^{-1} \sum_i \dot\mu(\psi^\top
\hat\xi_i) \hat\xi_i \hat\xi_i^\top\\ and \\\hat M(\psi) = \partial \hat
m / \partial \psi^\top\\. Then:

- BCA::

  \\\hat\psi\_{\mathrm{bca}} = \hat\psi - \hat I(\hat\psi)^{-1} \hat
  m(\hat\psi)\\

- BCM::

  \\\hat\psi\_{\mathrm{bcm}} = \hat\psi - (\hat I(\hat\psi) + \hat
  M(\hat\psi))^{-1} \hat m(\hat\psi)\\

- CS::

  solves \\n^{-1} \sum_i \phi_i(\psi) = 0\\ with \\\phi_i(\psi) =
  \hat\xi_i \\Y_i - \mu(\psi^\top \hat\xi_i)\\ - m_i(\psi)\\

- one-step::

  maximizes the integrated mixture likelihood \\\prod_i \sum\_\ell
  \pi\_\ell\\ \Pi\_{\hat z_i+1, \ell+1}\\ f(Y_i \mid Z = \ell, \psi)\\
  by automatic differentiation (RTMB).

Multicategory analogues replace \\(p\_{01}, p\_{10}, \pi)\\ with \\(\Pi,
\pi_z)\\ throughout; see Appendix B of Battaglia et al. (2025) for
explicit formulas.

*Inference.* All five estimators report an asymptotic variance through
[`vcov.mcglm`](https://ojalab.github.io/mismeasured/reference/vcov.mcglm.md):
\\A^{-1} C A^{-1}\\ for naive / BCA / BCM (Theorems on the two-step
expansion and on bias-corrected estimators under drifting
misclassification), \\J^{-1} S J^{-\top}\\ for CS (Z-estimator sandwich,
with \\J = -(\hat I + \hat M)\\, \\S = E\[\phi_i \phi_i^\top\]\\), and
the inverse Hessian of the integrated log-likelihood for one-step. See
`vcov_corrected` below for an alternative, more conservative BCA/BCM
sandwich.

## Required inputs by method

Only the chosen `method`s' inputs are validated; you do not have to
supply parameters for methods you are not fitting.

|  |  |  |
|----|----|----|
| **Method** | **Binary (K = 2)** | **Multicategory (K \> 2)** |
| `"naive"` | *none* | *none* |
| `"bca"`, `"bcm"`, `"cs"` | `c1` and `c2`, or `p01`/`p10`/`pi_z`, or `Pi` | `Pi` and `pi_z` |
| `"onestep"` (`fix_omega = FALSE`) | *none* (mixture weights estimated) | *none* |
| `"onestep"` (`fix_omega = TRUE`) | `p01`, `p10`, `pi_z` | `Pi`, `pi_z` |

*Auto-derivations performed by* `mcglm`*:*

- If `Pi` is supplied with `nrow(Pi) == 2`, `p01 = Pi[2, 1]` and
  `p10 = Pi[1, 2]` are extracted.

- If `Pi` and `z_hat` are available but `pi_z` is not, `pi_z` is
  estimated by Bayesian inversion of the observed proxy frequencies
  (`.mcglm_estimate_pi_z`).

- If `p01`, `p10`, `pi_z` are available but `c1`, `c2` are not, the
  latter are computed from the formulas above.

Supplying `c1` and `c2` directly is an "advanced" path that bypasses
these checks. The formula interface (`y ~ mc(z, Pi) + ...`) extracts
`Pi` automatically and is the recommended entry point.

## References

Nakamura, T. (1990). Corrected score function for errors-in-variables
models: methodology and application to generalized linear models.
*Biometrika*, 77(1), 127–137.
[doi:10.1093/biomet/77.1.127](https://doi.org/10.1093/biomet/77.1.127) .
Originating reference for the corrected-score concept.

Akazawa, K., Kinukawa, N. and Nakamura, T. (1998). A note on the
corrected score function adjusting for misclassification. *Journal of
the Japan Statistical Society*, 28(1), 115–123. Establishes existence of
a corrected score for misclassified discrete covariates and gives
explicit forms for normal, Poisson, and binary logistic GLMs;
foundational reference for the `"cs"` method.

Battaglia, L., Christensen, T., Hansen, S. and Sacher, S. (2025).
Inference for regression with variables generated by AI or machine
learning. *arXiv preprint arXiv:2402.15585*. Drifting-asymptotics
framework for linear models with AI/ML-generated covariates; originating
reference for the BCA/BCM and one-step estimators.

Warmuz, M. and Beresewicz, M. (2026). Bias corrections for GLMs with
misclassified covariates. Working paper. Extends the drifting-regime
BCA/BCM to GLMs, derives the corrected-score estimator under fixed
misclassification probabilities for binary and multicategory latent
regressors with explicit drift formulas, and gives the plug-in
robustness conditions used here when \\(\Pi, \pi)\\ are estimated from
auxiliary data.

## See also

[`simex`](https://ojalab.github.io/mismeasured/reference/simex.md) for
SIMEX and MC-SIMEX corrections;
[`mc`](https://ojalab.github.io/mismeasured/reference/mc.md) for the
formula-interface marker;
[`summary.mcglm`](https://ojalab.github.io/mismeasured/reference/summary.mcglm.md),
[`vcov.mcglm`](https://ojalab.github.io/mismeasured/reference/vcov.mcglm.md),
[`confint.mcglm`](https://ojalab.github.io/mismeasured/reference/confint.mcglm.md).

## Examples

``` r
# --- Formula interface (binary misclassification, Poisson) ---
set.seed(42)
n  <- 2000
Pi <- matrix(c(0.90, 0.10,    # Pi[2,1] = p01 = 0.10
               0.15, 0.85),   # Pi[1,2] = p10 = 0.15
             nrow = 2, byrow = FALSE)
z      <- rbinom(n, 1, 0.4)               # latent Z, prevalence pi = 0.4
z_hat  <- z
z_hat[z == 0] <- rbinom(sum(z == 0), 1, 0.10)        # flip 0 -> 1
z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)    # flip 1 -> 0
x1 <- rnorm(n)
y  <- rpois(n, exp(0.8 * z + -0.5 + 0.7 * x1))
df <- data.frame(y = y, z = factor(z_hat), x1 = x1)

fit <- mcglm(y ~ mc(z, Pi) + x1, data = df, family = "poisson",
             method = c("naive", "bca", "bcm", "cs"))
fit
#> 
#> Call:
#> mcglm(formula = y ~ mc(z, Pi) + x1, data = df, family = "poisson", 
#>     method = c("naive", "bca", "bcm", "cs"))
#> 
#> Family: poisson  |  n = 2000, K = 2, p = 3
#> Methods: naive, bca, bcm, cs
#> 
#> Coefficients:
#>              NAIVE    BCA      BCM      CS     
#> gamma         0.5653   0.7056   0.7538   0.7541
#> (Intercept)  -0.3709  -0.4394  -0.4629  -0.4675
#> x1            0.7121   0.7110   0.7107   0.7107
#> 
#> Degrees of Freedom: 2000 Total (i.e. Null);  1997 Residual
#> Null Deviance:     3663 
#> Residual Deviance: 2265  | AIC (naive): 5026 
summary(fit)
#> 
#> Call:
#> mcglm(formula = y ~ mc(z, Pi) + x1, data = df, family = "poisson", 
#>     method = c("naive", "bca", "bcm", "cs"))
#> 
#> Family: poisson  |  n = 2000, K = 2, p = 3
#> Methods: naive, bca, bcm, cs
#> z categories (Pi assumed in this order): 0 (baseline), 1
#> 
#> --- NAIVE ---
#>             Estimate Std. Error z value Pr(>|z|)    
#> gamma        0.56529    0.04385   12.89   <2e-16 ***
#> (Intercept) -0.37092    0.03535  -10.49   <2e-16 ***
#> x1           0.71207    0.02039   34.92   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> --- BCA ---
#>             Estimate Std. Error z value Pr(>|z|)    
#> gamma        0.70558    0.04446   15.87   <2e-16 ***
#> (Intercept) -0.43937    0.03749  -11.72   <2e-16 ***
#> x1           0.71102    0.02020   35.20   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> --- BCM ---
#>             Estimate Std. Error z value Pr(>|z|)    
#> gamma        0.75378    0.04485   16.80   <2e-16 ***
#> (Intercept) -0.46289    0.03831  -12.08   <2e-16 ***
#> x1           0.71065    0.02019   35.20   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> --- CS ---
#>             Estimate Std. Error z value Pr(>|z|)    
#> gamma        0.75410    0.05933   12.71   <2e-16 ***
#> (Intercept) -0.46751    0.04282  -10.92   <2e-16 ***
#> x1           0.71067    0.02035   34.92   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual deviance (naive): 2265  on 1997 degrees of freedom
#> AIC (naive): 5026
#> 
#> Bias correction (difference from naive):
#>              bca        bcm        cs       
#> gamma         0.140293   0.188490   0.188810
#> (Intercept)  -0.068454  -0.091968  -0.096587
#> x1           -0.001052  -0.001418  -0.001403
confint(fit, method = "cs")
#>                  2.5 %     97.5 %
#> gamma        0.6378256  0.8703766
#> (Intercept) -0.5514243 -0.3835923
#> x1           0.6707755  0.7505572

# Same problem with the AKN98 unbiased-surrogate corrected score
# (needs only Pi, not pi_z). The two CS variants give population-
# equivalent estimators and agree closely in large samples.
fit_akn <- mcglm(y ~ mc(z, Pi) + x1, data = df, family = "poisson",
                 method = c("cs", "cs_akn"))
coef(fit_akn, method = "cs_akn")
#>       gamma (Intercept)          x1 
#>   0.7533845  -0.4645779   0.7088654 

# --- Matrix interface, supplying p01/p10/pi_z directly ---
x_mat <- cbind(1, x1)
fit2 <- mcglm(y, z_hat = as.integer(z_hat), x = x_mat,
              family = "poisson", method = c("naive", "cs"),
              p01 = 0.10, p10 = 0.15, pi_z = 0.4)

# --- Multicategory misclassification (K = 3) ---
# \donttest{
set.seed(3)
n   <- 1000
Pi3 <- matrix(c(0.8, 0.1, 0.1,
                0.1, 0.8, 0.1,
                0.1, 0.1, 0.8), 3, 3, byrow = TRUE)
pi3 <- c(0.5, 0.3, 0.2)
Z   <- sample(0:2, n, replace = TRUE, prob = pi3)
Zh  <- vapply(Z, function(z) sample(0:2, 1, prob = Pi3[, z + 1]), 0L)
x1  <- rnorm(n)
y   <- rpois(n, exp(c(0, 1.0, -0.9)[Z + 1] + 0.8 - 0.7 * x1))
fit3 <- mcglm(y, z_hat = Zh, x = cbind(1, x1), family = "poisson",
              method = c("naive", "cs"), Pi = Pi3, pi_z = pi3)
coef(fit3, method = "cs")
#>     gamma1     gamma2     alpha0     alpha1 
#>  1.0131271 -1.1437257  0.8162136 -0.7088909 
# }
```
