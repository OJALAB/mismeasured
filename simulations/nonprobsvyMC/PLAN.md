# nonprobsvyMC — plan

Bridge package: `nonprobsvy`'s non-probability-survey API with
`mismeasured`'s misclassification-corrected GLMs, so that

```r
nonprob_mc(
  outcome        = y ~ mc(z, Pi) + x1 + x2,
  data           = nonprob_sample,
  svydesign      = prob_design,
  method_outcome = "mcglm",
  family_outcome = "poisson",
  mc_method      = "cs"          # naive | bca | bcm | cs
)
```

returns a `nonprob`-class object whose mass-imputation outcome model is
bias-corrected for misclassification in `z`.

## 1. Verified facts (nonprobsvy 0.3.0, mismeasured 0.7.0)

Checked against the installed `nonprobsvy` 0.3.0 namespace and
`mismeasured` 0.7.0 as on GitHub `main`
(`OJALAB/mismeasured@59c7540..db8a43e`, 2026-07-12), which ships the
ecosystem API this plan depends on (Section 4).

1. **`nonprob()` dispatch is a hard-coded `switch`.** Inside
   `nonprobsvy:::nonprob_mi`:
   `switch(method_outcome, glm = method_glm, nn = method_nn, pmm = method_pmm, npar = method_npar)`.
   A third-party package **cannot** register a new `method_outcome`
   string against the CRAN version. Any short-term bridge must wrap, not
   plug in.

2. **The outcome-method contract is clean and small.**
   `method_glm(y_nons, X_nons, X_rand, svydesign, weights, family_outcome,
   start_outcome, vars_selection, pop_totals, pop_size, control_outcome,
   control_inference, verbose, se)` returns a `"nonprob_method"` object:
   `model_fitted, y_nons_pred, y_rand_pred, coefficients, svydesign
   (updated with predictions), y_mi_hat, var_prob, var_nonprob,
   var_total, model, family`. A `method_mcglm()` against this contract
   would be straightforward, but the dispatch is closed upstream — hence
   the standalone-wrapper architecture (Section 3). The contract remains
   the blueprint for what `nonprob_mc()`'s `outcome` slot must mimic.

3. **`method_glm`'s variance is exactly the Kim et al. (2021)
   linearization** that Theorem `thm:mass-imp-cs` of
   `theory/ml-glm/paper.tex` generalizes: `var_prob` comes from
   `survey::svymean` on the predicted values; `var_nonprob` is
   `n^{-2} * sum(resid^2 * (X c)^2)` with
   `c = ginv(I_naive) %*% G_P`. The corrected-score analogue replaces
   `I_naive` with `I_hat + M_hat` and the residual score with
   `phi_i = xi_i(y_i - mu_i) - m_i` — already implemented and
   Monte-Carlo-validated in
   `simulations/sim-paper-nonprob-misclass/functions.R::mi_variance()`
   (coverage 0.94–0.95 at the nominal 0.95).

4. **The `nonprob` object fields** (constructor at the end of
   `nonprob_mi`): `data, X, y, R, ps_scores, case_weights, ipw_weights,
   control, output, SE, confidence_interval, nonprob_size, prob_size,
   pop_size, pop_size_fixed, pop_totals, pop_means, outcome (list of
   nonprob_method), selection, boot_sample, svydesign, ys_rand_pred,
   ys_nons_pred, ys_resid`. Populating these lets `print`, `summary`,
   `confint`, `nobs`, `weights`, `extract` work unchanged.

5. **`mismeasured`'s formula machinery**: `mc(z, Pi)` terms are parsed
   by the exported `mc_parse_formula(formula, data, env)` (since 0.7.0),
   which returns `list(y, z_hat (0-based), x (model matrix of the
   remaining terms), Pi, K)`; exactly one `mc()` term is required;
   factor `z` is re-coded `as.integer - 1`; the response may be absent
   (prediction data). The `mc()` call must be stripped before
   `model.matrix()` is ever called — `nonprobsvy`'s own
   `make_model_frame`/`merge_formulas` would choke on it, so the bridge
   must intercept the outcome formula *before* handing anything to
   `nonprobsvy`.

6. **Constraints inherited from `mismeasured`**: canonical links only
   (poisson/binomial/gaussian — matching `family_outcome` options);
   no offsets; frequency weights supported (`weights` → `case_weights`);
   `mc_method` ∈ {naive, bca, bcm, cs}; K ≥ 2 categories; multinomial
   responses excluded (BCA/BCM/CS unsupported there). Since 0.7.0,
   `coef()`/`vcov()`/`predict()` all default to the same (last-fit)
   method, and `predict(newdata)` accepts a data frame (formula fits)
   or `list(z_hat, x)` (matrix fits).

7. **Estimand caveat**: the theory (and the bridge, initially) assumes
   the **probability sample carries the true `z`** while the
   non-probability sample has the misclassified proxy. If `S_P` also
   only has the proxy, predicting requires marginalizing
   `E[mu | z_hat]` over `P(z | z_hat)` (needs `Pi` *and* `pi_z`) — a
   well-defined extension, but not v0.1 (Section 5).

### Where the correction matters (verified in `poc.R`)

For the **overall** population mean, the intercept of the naive fit
absorbs most of the misclassification bias (Monte Carlo, 200 reps,
p01 = 0.15 / p10 = 0.20, Poisson: naive bias −0.005 vs CS +0.002).
The headline gain is **domain estimation**: predicting the `z = 1`
subpopulation with the attenuated naive `gamma` (0.52 vs true 0.80)
gives a domain-mean bias of **−0.21 (~11%)**, while CS is unbiased
(+0.002). Domain/subset estimation should therefore be a first-class
feature of `nonprob_mc()` from M1, not an afterthought.

## 2. What the package does (scope)

- **v0.x scope: mass imputation (MI) only.** `nonprob(outcome=...)`
  without `selection=` — the estimator our Theorem `thm:mass-imp-cs`
  covers. IPW is untouched by outcome misclassification; DR needs a
  corrected outcome model *plus* a selection model possibly involving
  the misclassified covariate — deferred (Section 5).
- Point estimator: `mu_hat = svymean` of `mu(psi_hat' xi_i)` over the
  probability sample (this is the Hájek form — `survey`'s default —
  consistent with what `method_glm` already does via
  `svydesign_updated`).
- Variance: `var_prob` from `survey` (works for any design, not just
  the Poisson sampling of our simulation) + `var_nonprob` from the
  corrected-score linearization (port of `mi_variance()`, with
  `G_P = N_hat^{-1} sum w_i mudot_i xi_i` and `A = I_hat + M_hat` for
  CS, `A = I_hat` for naive/BCA/BCM per Theorem `thm:mis-bc`).

## 3. Architecture: standalone `nonprob_mc()` wrapper (DECIDED)

The package is a standalone wrapper — the previously sketched
alternatives (an upstream `method_outcome` extension point in
`nonprobsvy`; masking `nonprob()`) are dropped from the plan: the
wrapper needs no upstream changes, ships independently of `nonprobsvy`
release cycles, and nothing in it is wasted if an upstream hook ever
appears later.

New package `nonprobsvyMC`, `Imports: nonprobsvy, mismeasured
(>= 0.7.0), survey, stats`. One user-facing function `nonprob_mc()`
mirroring the `nonprob()` signature (MI arguments only) plus
`mc_method` (which correction: `"naive"`, `"bca"`, `"bcm"`, `"cs"`)
and the misclassification arguments (`Pi`/`pi_z`/`p01`/`p10`,
validation data — Section 4a).

Flow (public mismeasured 0.7.0 API only, zero `:::` — see the CRAN
constraint below):

1. Parse `outcome` with `mismeasured::mc_parse_formula()`; obtain
   `(y, z_hat, x, Pi, K)` and the *sanitized* formula
   (`y ~ z + x1 + x2`).
2. Fit `mismeasured::mcglm()` on the non-probability sample
   (all four methods cheaply; report the one in `mc_method`, keep the
   rest for diagnostics like the existing `mcglm` output).
3. Predict on the probability sample with
   `predict(fit, newdata = svydesign$variables, type = "response")`
   (true `z` there), `update()` the svydesign, `svymean()` → point
   estimate and `var_prob`.
4. `var_nonprob` from the corrected-score linearization via
   `sandwich::bread()`/`sandwich::meat()` on the fit
   (`t(G_P) %*% B %*% S %*% t(B) %*% G_P / n`, with `G_P` from the
   design weights and `family()$mu.eta`).
5. Assemble the `nonprob`-class object with the fields from Fact 4
   (`selection = NULL`, `outcome` = a `nonprob_method`-shaped list with
   `model = "mcglm"`), so all of `nonprobsvy`'s S3 methods work.

Steps 1–4 are demonstrated end-to-end in `poc.R` in this folder
(naive `nonprob()` vs the corrected pipeline on a DGP with a
misclassified binary covariate), running entirely on exported API.

### CRAN constraint: no `:::` in the package

`R CMD check --as-cran` NOTEs every cross-package `:::` call
("Unexported objects imported by ':::' calls"); `?":::"` calls it a
design mistake and CRAN expects it fixed at submission. The
same-maintainer tolerance technically applies here (both upstreams are
ours), but relying on it just trades a check NOTE for silent breakage
on the next internal refactor. `getFromNamespace()`-style workarounds
evade the check, not the policy. Hard rule for nonprobsvyMC: **every
upstream symbol it uses must be exported** — from `mismeasured` via the
Section 4 API (shipped in 0.7.0), from `nonprobsvy` via the exported
`nonprob()` only. `poc.R` in this folder is the acceptance check: it
runs the whole pipeline with zero `:::` calls.

## 4. Changes in mismeasured — SHIPPED in 0.7.0 (GitHub `main`)

Design decision: do NOT export the raw internals
(`.mcglm_compute_m_bin`, `.mcglm_compute_Ihat`, ... — that freezes
implementation-level signatures forever) and do NOT vendor the math
into nonprobsvyMC (two copies of verified code will drift). Instead an
**object-level accessor API**, registered on established ecosystem
generics rather than package-specific names. All of the following is
now implemented in `mismeasured` (see `NEWS.md` 0.7.0, files
`R/mcglm-sandwich.R`, `R/mcglm-marginaleffects.R`):

1. `mc_parse_formula(formula, data, env)` — exported `mc()` parser
   (single source of truth for `mc()` semantics).
2. `estfun.mcglm` — delayed S3 registration on `sandwich::estfun`
   (sandwich in Suggests): n x p matrix of per-observation estimating
   functions at the chosen method's coefficients (`phi_i` for CS, the
   proxy score for naive/BCA/BCM). Gives `S` via `sandwich::meat()`.
3. `bread.mcglm` — on `sandwich::bread`: `(I_hat + M_hat)^{-1}` for CS,
   `I_hat^{-1}` otherwise. Gives `A^{-1}`. **Caveats** (documented in
   `?estfun.mcglm`): `sandwich::sandwich()` does not forward `method`
   to the bread and assembles `B %*% meat %*% B` without a transpose,
   which is wrong for the asymmetric CS bread — assemble
   `B S B' / n` manually (what nonprobsvyMC's variance does anyway) or
   use `vcov(fit, method =)`, which is already the correct sandwich.
4. `predict.mcglm(object, newdata)` — data frame for formula fits
   (predicts at whatever category `newdata` carries, i.e. the true `z`
   on the probability sample), `list(z_hat =, x =)` for matrix fits.
   Also covers the M3 validation-arm D-matrix without further exports:
   the toggle gap `delta_i` is
   `predict(newdata, z = 1) - predict(newdata, z = 0)`.
5. `coef.mcglm()` default changed to return the default method's
   vector (`method = "all"` for the old list) — required by
   `lmtest::coeftest` and marginaleffects, and consistent with
   `vcov()`/`predict()`. `lmtest::coeftest(fit)` now works as-is.
6. marginaleffects support (formula fits): `get_coef`/`set_coef`/
   `get_vcov`/`get_predict` methods registered on the marginaleffects
   generics plus the `marginaleffects_model_classes` opt-in set in
   `.onLoad()` — `avg_comparisons()`/`avg_slopes()` etc. work with
   delta-method SEs propagated from `vcov(fit)`.
7. The fitted object stores `$misclass` (`c1`, `c2`, `p01`, `p10`,
   `pi_z`, `Pi`) so the accessors are self-contained.

With these, nonprobsvyMC needs only `stats::<family>()$mu.eta` (for
`G_P`), `sandwich::` accessors, and `nonprobsvy::nonprob()` — zero
`:::` calls. `poc.R` in this folder has been rewritten against the
public API only and is the acceptance check. Rejected alternatives:
exporting raw internals (wrong altitude), vendoring (verified-code
drift), storing component matrices in the fitted object (bloat, fixes
the evaluation point), a shared core package (overkill).

## 4a. Estimated misclassification probabilities

When `Pi` is not known, the M3 milestone needs more than a point
estimate of the misclassification parameters — the extra sampling
uncertainty enters both `psi_hat` and the mass-imputation variance
(paper subsection "Corrected-score under estimated misclassification
probabilities"; prototyped and Monte-Carlo-validated in
`simulations/sim-paper-nonprob-misclass/functions.R::mi_variance()`).

**Additional information the API must accept** (one of):

1. **A validation sample** carrying paired `(z, z_hat)` — the richest
   input. Proposed argument:
   `validation = list(data = <df>, true = "z", proxy = "z_hat")`, plus
   a flag for its provenance:
   - `subset_of = "nonprob"` (S_V drawn from S_NP): the variance needs
     the extra `Cov(phi_i, s_i)` term, estimable because `phi_i` and
     `s_i` are jointly observed for i in S_V;
   - independent external sample: covariance term drops, only `V_V`
     (with `V(s_i)` estimated from S_V) is added.
   From S_V the bridge computes `p01_hat`, `p10_hat` (binary) or
   `Pi_hat` (multicategory, cell frequencies), `pi_hat`, and the
   influence functions `s_i`. The `D = -d mbar / d eta'` matrix comes
   from the toggle-gap trick via `predict()` (binary); for K > 2 the
   `d m / d Pi_{jl}` derivatives are linear in the `mu` differences and
   computable from K `predict()` calls per category.
2. **An externally estimated `Pi_hat` with its covariance matrix**
   (e.g. from a published confusion matrix or an audit study):
   `Pi = Pi_hat, Pi_vcov = <matrix>`. Then `V_V` is the delta-method
   term `G' H^{-1} D Pi_vcov D' H^{-T} G` and no covariance term
   applies. If only `Pi_hat` is supplied without `Pi_vcov`, treat as
   known and *warn* that the reported variance conditions on `Pi_hat`
   (anti-conservative).

**Theory caveats to surface in the docs**: consistency needs
`sup ||Phi_n* - Phi_n|| = o_P(1)` and asymptotic unbiasedness the
stronger `o_P(n^{-1/2})` rate (Prop. `prop:est_prob`); with rare
misclassification the effective validation information is the *number
of misclassified units observed* in S_V, not `n_V` — with small
`p01`/`p10` the estimated-`Pi` variance term can dominate `V_NP`
(demonstrated in the drifting-regime simulation, where `n_V = 0.1
n_NP` sits at the boundary of the required rate).

**Sample-size guidance** (to include in the vignette): for the binary
case, `s_i` has variance `p01(1-p01)/(1-pi)` and `p10(1-p10)/pi`
divided by `n_V` — a quick pre-check `n_V * p_hat >= ~20` per error
type should gate a warning.

## 5. Roadmap

| milestone | content |
|---|---|
| M0 (done) | POC in this folder: corrected MI pipeline against `nonprobsvy` naive MI on simulated data |
| M1 | Package skeleton; `nonprob_mc()` point estimation, known `Pi`, binary `z`, MI only; `nonprob`-object compatibility (`print`/`summary`/`confint`) |
| M2 | Variance (`var_prob` via survey + corrected `var_nonprob`); tests reusing `simulations/sim-paper-nonprob-misclass` DGP as integration test (coverage check vs known truth) |
| M3 | Estimated misclassification (Section 4a): `validation =` / `Pi_vcov =` inputs, `V_V` + `Cov(phi, s)` terms (already prototyped in `mi_variance()`), small-validation-sample warning; multicategory `z` (K > 2) |
| M4 | Proxy-only probability sample (`E[mu | z_hat]` marginalization using `Pi`, `pi_z`); bootstrap variance option (mirror `control_inf(var_method = "bootstrap")`) |

## 6. Risks / open questions

- **`nonprob` S3 internals may assume selection-side fields** even for
  MI-only objects (e.g. `summary` touching `ps_scores`). M1 must test
  every exported S3 method against the assembled object; fall back to
  our own `print` method if any assumption is baked in too deep.
- **`survey` design generality**: `var_prob` from `svymean` covers
  stratified/clustered designs automatically, but the paper's variance
  theorem is stated for Poisson sampling; document that `var_prob` is
  design-based and exact-form only where `survey` is.
- **Naming (RESOLVED)**: `method_outcome = "mcglm"` + a separate
  `mc_method` argument for the correction (`"naive"`, `"bca"`,
  `"bcm"`, `"cs"`) — the `mc_` prefix underlines the misclassification
  correction and matches the `mc()` / `mc_parse_formula()` family;
  no overloaded `"mcglm-cs"` strings.
- **Weights semantics**: `nonprobsvy::nonprob` `case_weights` are
  frequency weights — matches `mcglm(weights=)`; verify no rescaling
  happens upstream.
