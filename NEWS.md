# leftcens 0.4.1

* New vignette "Preparing your own data for leftcens" --- a step-by-step
  on-ramp from a wide laboratory export through per-analyte `cens_data`
  construction, descriptive summaries, and imputation, with recipes for the
  common input formats: pre-split `left`/`right` columns, a value column plus a
  qualifier-flag column per analyte (e.g. `A`, `A.flag`), and single-column
  value-plus-censored-flag data.

# leftcens 0.4.0

* `as.data.frame()` methods for `desc_np()` and `desc_sp()` fits return tidy
  data frames (a quantile table and a coefficient table, respectively), so a
  descriptive summary can be passed directly to a table-formatting package such
  as `gt`, `flextable`, `knitr::kable()`, or `DT`.
* broom-style `tidy()` and `glance()` methods for `desc_np()` and `desc_sp()`
  fits (via the lightweight `generics` package), for tidyverse-ecosystem
  integration (e.g. `gtsummary`, `modelsummary`). `tidy()` gives the
  component-level table; `glance()` a one-row model summary.
* `desc_np()` now records the `detection_limit` and `quantitation_limit` of the
  data, and `quantile.cens_np_fit()` reports any quantile that falls below the
  quantitation limit as `NA` (an estimate there is an extrapolation into the
  censored region, not a reliably quantified value). Pass `ql = 0` for the raw
  NPMLE estimates, or `ql = fit$detection_limit` to mask below the detection
  limit instead. The `print()` method shows the limits.

# leftcens 0.3.0

* New `preflight_reliability()` and `preflight_from_data()`: a Monte Carlo
  pre-flight that simulates a study design (or reads it from your `cens_data`),
  imputes it, and reports the per-analyte bias and multiple-imputation coverage
  you can expect --- flagging which analytes are reliably imputable (small bias,
  calibrated coverage, and non-detects below 50%) before you trust imputation on
  real data.
* New vignette "How much can you trust the imputation?" with practical guidance
  (model choice, `iters_all`/`M`, the ND < 50% reliability target, and caveats),
  synthesising the Monte Carlo validation study.

# leftcens 0.2.0

* **`gsimp_impute()` now uses `imp_model = "tobit"` by default** (previously
  `"ridge"`). The Tobit model is an interval-censored Gaussian conditional model
  (via `survival::survreg()`) fit on observed *and* censored rows, so it is not
  subject to the selection bias of the observed-only models that fit only on the
  (upper-truncated) detected values. In a Monte Carlo validation it holds
  mean-recovery bias near zero and multiple-imputation coverage near nominal
  down to ~55% non-detects, where ridge (reliable only to ~25% non-detects) is
  badly biased and its intervals fail to cover. **This changes default output:**
  imputations and any downstream summaries will differ from 0.1.0 and be less
  biased. Pass `imp_model = "ridge"` to recover the old behaviour.
* The Tobit model reduces the predictor dimension by PCA when analytes outnumber
  samples, so it applies to wide (e.g. metabolomics) data rather than falling
  back to ridge. It is several times slower than ridge.

# leftcens 0.1.0

First release. Descriptive statistics and imputation for left- and
interval-censored measurement data, built on a single three-tier interval model
(non-detect / detected-not-quantified / quantified).

## Data model

* `as_interval_data()` and `as_interval_data_from_single()` construct `cens_data`
  objects; single-limit (Helsel-style) left censoring is handled as a degenerate
  two-tier case of the same model.

## Descriptive statistics

* `desc_np()` fits the non-parametric NPMLE (Turnbull via `survival`, or Wang via
  `icenReg`), with `quantile()` and `plot()` methods.
* `desc_sp()` fits semi-parametric proportional-hazards / proportional-odds
  baselines (via `icenReg`), with `coef()`, `summary()`, and `plot()` methods.

## Imputation

* `build_bounds()` assembles per-cell `(lo, hi)` bound matrices from `cens_data`
  columns, so non-detect and detected-not-quantified cells of the same analyte
  are imputed within their respective bands.
* `gsimp_impute()` performs Gibbs-sampler imputation constrained to those bounds,
  with a pluggable conditional model (`ridge`, `lm`, `glmnet`, or a user
  function). `rnorm_trunc()` provides the underlying truncated-normal draw.

## Provenance

* Released under the MIT license. The imputation is an independent, clean-room
  implementation of the GSimp method (Wei et al. 2018,
  <doi:10.1371/journal.pcbi.1005973>), written from the published algorithm
  rather than ported from the CC BY-NC-SA GSimp source.
