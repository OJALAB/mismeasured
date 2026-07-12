# nonprobsvyMC

Misclassification-corrected mass imputation for non-probability
samples: bridges [nonprobsvy](https://github.com/ncn-foreigners/nonprobsvy)
with [mismeasured](https://github.com/OJALAB/mismeasured).

When the outcome model of a mass-imputation estimator contains a
discrete covariate observed with misclassification (e.g. an AI/ML-coded
label in a big non-probability source), the naive estimator is biased —
most visibly for **domain estimation**, where the attenuated
coefficient of the misclassified covariate propagates directly into the
domain means. `nonprob_mc()` corrects the outcome model with the
BCA / BCM / corrected-score estimators of `mismeasured::mcglm()` and
returns a regular `nonprob` object.

## Usage

```r
library(nonprobsvyMC)
library(survey)

Pi <- matrix(c(0.9, 0.1,     # Pi[j, l] = P(z_hat = j - 1 | z = l - 1)
               0.15, 0.85), 2, 2)

fit <- nonprob_mc(
  data           = nonprob_sample,        # y, proxy z, covariates
  outcome        = y ~ mc(z, Pi) + x1 + x2,
  svydesign      = prob_design,           # TRUE z, no y needed
  family_outcome = "poisson",
  mc_method      = "cs"                   # naive | bca | bcm | cs
)

fit                # print.nonprob
summary(fit)
confint(fit)
fit$outcome[[1]]$model_fitted   # the underlying mcglm fit
```

The probability sample (`svydesign`) must carry the **true** category
of the `mc()` variable; the non-probability sample carries the proxy.
Misclassification parameters are passed as in `mcglm()`: `Pi` (inside
`mc()` or as an argument), or `p01`/`p10`/`pi_z` for a binary `z`.

The variance is `var_prob` (design-based, from **survey**) plus
`var_nonprob` (corrected-score linearization of the outcome-model
uncertainty, assembled from `sandwich::bread()`/`meat()` of the
`mcglm` fit).

For domain estimation, the returned object's `svydesign` element is the
updated design with predictions in `.y_mc_hat`:

```r
svyby(~.y_mc_hat, ~region, fit$svydesign, svymean)   # design variance only
```

## Installation

```r
# install.packages("remotes")
remotes::install_github("OJALAB/mismeasured")   # >= 0.7.0
install.packages("nonprobsvy")                  # >= 0.3.0
remotes::install_github("OJALAB/nonprobsvyMC")
```

## Scope and roadmap

Mass imputation only (no IPW/DR), canonical links
(gaussian/binomial/poisson), K >= 2 categories, known
misclassification parameters. Planned: validation-sample estimated
misclassification with the extra variance terms, proxy-only probability
samples, bootstrap variance. See `PLAN.md` in the parent folder of the
development repo for the design rationale and theory references.
