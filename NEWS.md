# leftcens (development version)

* `gsimp_impute()` gains an experimental `imp_model = "tobit"`: an
  interval-censored Gaussian conditional model (via `survival::survreg()`) fit
  on observed *and* censored rows. Unlike the observed-only models it is not
  subject to the selection bias of fitting on the (upper-truncated) detected
  values, and in simulation it removes essentially all of the mean-recovery bias
  even at heavy censoring. It falls back to ridge for wide or thin designs.

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
