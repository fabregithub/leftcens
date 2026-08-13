# leftcens 0.9.0

* **New `impute_censored_conditional()`**: multiple imputation of a single
  left-/interval-censored variable conditional on a *general* design matrix (an
  outcome, mixed-type covariates, other analytes) — the general-predictor
  complement to `gsimp_impute()`, which conditions on other analytes only. Handles
  skew via the sinh-arcsinh margin, respects per-cell censoring bounds, and does
  proper MI (draws the margin and regression parameters). This is the reusable
  core for fully-conditional-specification loops and congenial imputation of a
  censored *exposure* using the outcome.
* **Exported the sinh-arcsinh margin API** (previously internal): `fit_shash_margin()`
  (interval-censored MLE), `draw_margin()` (proper-MI parameter draw), and the
  latent-scale transforms `x_to_z()` / `z_to_x()`.

# leftcens 0.8.0

* **`imp_model = "copula"` is now the default** for `gsimp_impute()`,
  `gsimp_mi()`, and the `preflight_*()` functions (previously `"tobit"`). Real
  measurement data are often more skewed than log-normal, and the shape is usually
  undiagnosable under censoring; the copula is near-unbiased across skew and
  matches `"tobit"` on log-normal data, so it is the safer default. Set
  `imp_model = "tobit"` (or `"ridge"`) to restore the previous, faster behaviour
  for known-log-normal or very large/wide problems.
* New `gsimp_mi()`: multiple imputation with Rubin's-rules pooling of the
  per-analyte means, reporting the fraction of missing information (FMI). Draws
  imputations in parallel (`n_cores`), uses the copula hybrid automatically
  (plug-in point estimate + drawn-margin variance), can choose the number of
  imputations adaptively (`adaptive = TRUE`), and returns the completed datasets
  (`return_imputations = TRUE`) for a downstream analysis.
* Guidance vignette and README reframed around the copula default and a
  multiple-imputation workflow (a single imputation under-covers; choose the
  number of imputations by FMI --- about 30 for the copula, or adaptive). The
  README states the package's scope: descriptive statistics and marginal-
  distribution reconstruction, not association modelling with a censored variable.
* `imp_model = "copula"` now handles an entirely-censored analyte gracefully
  (with no observed values the margin falls back to the censoring bounds).

# leftcens 0.7.0

* New skew-robust imputation model `imp_model = "copula"` for `gsimp_impute()`.
  Right-skewed analytes bias the default `"tobit"` mean-imputation downward under
  censoring; `"copula"` fits a flexible (sinh-arcsinh) margin per analyte by
  interval-censored MLE, imputes in latent-Gaussian space with the same censored
  sampler, and back-transforms --- keeping the mean near-unbiased and its
  multiple-imputation intervals calibrated under strong skew, across a range of
  dependence structures. Opt-in; the log-normal default `"tobit"` is unchanged.
* `gsimp_impute()` gains a `margin_draw` argument (copula only). The recommended
  multiple-imputation workflow uses the plug-in margin for the point estimate and
  drawn-margin imputations for the between-imputation variance, so Rubin's-rules
  intervals stay calibrated at heavy censoring without a point-estimate bias.
* `preflight_from_data()` now estimates each analyte's marginal skewness (by an
  interval-censored fit, not the imputed values) and reports the skew-induced mean
  bias instead of assuming log-normality, recommending `imp_model = "copula"` when
  it detects material right-skew. `preflight_reliability()` accepts a per-analyte
  `skew` vector.
* The imputation-guidance vignette gains a skew-tolerance rule (the mean is
  reliable while marginal skewness times the non-detect fraction stays below
  ~0.18) and the copula workflow.

# leftcens 0.6.1

* The "Reporting censored-data results" vignette's summary table now reports a
  range of quantiles (10th, 25th, median, 75th, 90th); those below the
  quantitation limit come back as `NA` (read as `< LCMRL`), making the effect of
  censoring on the lower tail explicit, with a note that "detected" (>= MDL)
  differs from "quantified" (>= LCMRL).

# leftcens 0.6.0

* New vignette "Reporting censored-data results" --- how to present
  censored-measurement results honestly (detection frequency and limits, when the
  median/mean are estimable, `< LCMRL` reporting, imputation disclosure with
  Rubin pooling, and citation), complementing the reliability guidance.

# leftcens 0.5.0

* Two example data sets, `groundwater` (three-tier `.left`/`.right` layout) and
  `surfacewater` (value + qualifier-flag layout), so the two common lab-export
  shapes can be tried directly (`data(groundwater)`).
* Descriptive statistics are now data-frame-first, mirroring the imputation
  side. `as_cens_list()` turns a wide laboratory table into a named list of
  `cens_data` (one per analyte) --- supporting both `.left`/`.right` column pairs
  and value-plus-qualifier-flag columns --- and that list feeds `desc_np()`,
  `desc_sp()`, and `build_bounds()` directly.
* `desc_np()` and `desc_sp()` now accept a list (or `cens_list`) of analytes and
  return a `cens_np_fits` / `cens_sp_fits` collection. `print()`,
  `as.data.frame()`, `tidy()`, `glance()`, `quantile()` (non-parametric), and
  `plot()` (non-parametric) combine the analytes into one table or overlay ---
  no more `lapply()`.

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
