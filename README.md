
<!-- README.md is generated from README.Rmd. Please edit that file -->

# leftcens

<!-- badges: start -->

[![R-CMD-check](https://github.com/fabregithub/leftcens/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fabregithub/leftcens/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/fabregithub/leftcens/actions/workflows/pkgdown.yaml/badge.svg)](https://fabregithub.github.io/leftcens/)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License:
MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**leftcens** provides descriptive statistics and imputation for left-
and interval-censored measurement data — the kind produced by analytical
and environmental laboratories, where a result may fall below a
reporting threshold rather than being a clean number.

A measurement is reported in one of three tiers, and each is really an
*interval* the true value lies inside:

| Report                   | What is known             | Interval         |
|--------------------------|---------------------------|------------------|
| Non-detect               | below the detection limit | `(0, MDL)`       |
| Detected, not quantified | between MDL and LCMRL     | `(MDL, LCMRL)`   |
| Quantified               | an actual measurement     | `(value, value)` |

One data representation (`cens_data`) underlies everything, so the
interval semantics are defined exactly once and shared by two capability
modules:

- **Descriptive statistics** — non-parametric (Turnbull NPMLE) and
  semi-parametric (proportional-hazards / proportional-odds) summaries
  of interval-censored data.
- **Imputation** — Gibbs-sampling imputation that draws each censored
  cell from within its own interval, extended to accept *per-cell*
  bounds so non-detect and detected-not-quantified values in the same
  analyte are handled correctly.

Classic single-limit (Helsel-style) left censoring is treated as a
degenerate two-tier case of the same model, not a separate mode.

## Installation

leftcens is not yet on CRAN. Install the development version from
GitHub:

``` r
# install.packages("remotes")
remotes::install_github("fabregithub/leftcens")
```

## Getting started

`leftcens` ships two example lab exports; `as_cens_list()` turns either
into one `cens_data` per analyte, which flows into the descriptive and
imputation tools:

``` r
library(leftcens)

# `groundwater` is bundled with the package (a wide ".left"/".right" lab export)
cl <- as_cens_list(groundwater)               # one cens_data per analyte

quantile(desc_np(cl), probs = c(0.5, 0.9))    # per-analyte quantiles (NA below the QL)
#>        50%    90%
#> As      NA  7.514
#> Cd      NA  0.300
#> Pb      NA  1.710
#> Cu  6.6005 19.694
#> Zn 10.8895 58.264

filled <- gsimp_impute(build_bounds(cl))      # impute the censored cells (tobit model)
```

Descriptive summaries are table-ready via `as.data.frame()` (`gt`,
`flextable`, `knitr::kable()`, …) and broom-style `tidy()` / `glance()`.

The four vignettes (also on the [package
website](https://fabregithub.github.io/leftcens/)) carry the detail:

- **Getting started** — the interval data model and the descriptive and
  imputation tools.
- **Preparing your own data** — from a wide lab export (interval *or*
  flagged) to results, using the bundled `groundwater` / `surfacewater`
  datasets.
- **Per-cell imputation bounds** — the Gibbs sampler and a recovery
  demonstration.
- **How much can you trust the imputation?** — practical guidance and
  the reliability limit (imputation stays trustworthy to roughly 50%
  non-detects; the full Monte Carlo study is in `validation/`).

## License and provenance

leftcens is released under the **MIT license**.

The imputation is an independent, clean-room implementation of the GSimp
*method* (Wei R. *et al.* 2018, *PLoS Comput. Biol.* 14(1):e1005973,
[doi:10.1371/journal.pcbi.1005973](https://doi.org/10.1371/journal.pcbi.1005973)),
written from the published algorithm rather than ported from the (CC
BY-NC-SA) GSimp source. Please cite the GSimp paper when you use the
imputation.
