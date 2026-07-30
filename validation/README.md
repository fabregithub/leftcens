# Monte Carlo validation study

A standalone, reproducible study of `leftcens::gsimp_impute()`. It is **not**
part of the package build or CI test suite — it is a heavier simulation study
intended as the basis for a methods/validation write-up.

## What it does

For a grid of data-generating conditions, and many replications each, it:

1. simulates ground-truth log-concentrations (multivariate, exchangeable
   correlation `rho`, optional sinh-arcsinh skew for model misspecification);
2. censors them into the three-tier structure (non-detect / detected-not-
   quantified / quantified) at analyte-specific limits;
3. imputes with `gsimp_impute()` and, for comparison, naive detection-limit
   substitution;
4. measures recovery against the known truth.

### Metrics

| Metric | Meaning |
|---|---|
| `gs_bias`, `naive_bias` | mean bias of recovered analyte means (log scale) |
| `gs_rmse`, `naive_rmse` | RMSE of the recovered mean across replications |
| `*_cor_mae` | mean absolute error of recovered between-analyte correlations |
| `mi_coverage` | 95% multiple-imputation (Rubin's rules) CI coverage of the true mean — the key calibration metric; ~0.95 is well-calibrated |

Conditions varied: sample size `n`, analyte count `p`, correlation `rho`,
censoring fraction (`nd_frac` + `dnq_frac`), and marginal `skew`.

## Running it

From the package root, with `leftcens` installed:

```bash
# quick sanity config
Rscript validation/mc_validation.R

# full publication grid (slower)
CONFIG=full N_REP=500 M=50 Rscript validation/mc_validation.R
```

Environment variables: `CONFIG` (`quick`|`full`), `N_REP`, `M`, `IMP_MODEL`
(`ridge`|`glmnet`|`lm`).

Results are written to `validation/results/` (`latest.rds`, a timestamped
`.rds`, and two `.csv` files). These are git-ignored — regenerate them by
running the script.

Then render the report:

```r
rmarkdown::render("validation/report.Rmd")
```

## Files

- `mc_validation.R` — the engine: simulation, censoring (per-analyte limits
  supported), metrics, runner, summariser, and a `main` block. Source it to
  reuse the functions (`run_validation()`, `summarise_validation()`, ...).
- `detection_rate_standard.R` — sweeps the detection rate x correlation and
  derives the lowest detection rate that is still "reliable" (bias and MI
  coverage criteria).
- `heterogeneous_limits.R` — tests whether well-observed "anchor" analytes
  rescue heavily-censored "target" analytes as correlation increases.
- `run_full.R` — runs the full 72-scenario grid in batches, accumulating
  replications across calls (`MODE=batch` ... then `MODE=summarise`). Use when a
  single high-replication pass exceeds the available run time.
- `tobit_vs_ridge.R` — paired comparison of the ridge and tobit conditional
  models across the detection-rate sweep (same data, only the model differs).
- `report.Rmd` — reads `results/latest.rds` and produces tables and figures.
- `results/` — generated outputs (git-ignored).

## Findings so far

**A rough detection-rate standard (uniform limits, n=200, p=6, log-normal).**
Reliability requires |mean bias (log)| <= 0.10 AND MI 95% coverage >= 0.90.

| Detection rate | Non-detects | Reliable? | Note |
|---|---|---|---|
| >= 80% | <= 20% | yes | point estimates and calibrated MI inference both hold |
| 70-80% | 20-30% | point only | bias < 0.10 but MI coverage already below nominal |
| < 70% | > 30% | no | bias grows, coverage collapses |

Two structural results:

1. **Coverage is the binding constraint, not bias.** Point recovery degrades
   gracefully; MI coverage falls off a cliff, because the small residual upward
   bias exceeds the (correctly narrow) MI standard error. The reliable limit is
   stricter for *inference* (>= 80% detection) than for *point summaries*
   (>= 70%).

2. **Correlation and auxiliary well-observed analytes do NOT extend the reliable
   range** (`heterogeneous_limits.R`). Even with anchor analytes at 5% ND and
   correlation up to 0.8, heavily-censored (40% ND) target analytes stay
   unreliable (coverage ~0.04 -> ~0.17, never near 0.95). The reason: each
   analyte's conditional imputation model is fit on that analyte's *observed*
   (upper-truncated) values, so it inherits selection bias no matter how well
   the neighbours are observed. **Reliability is governed by each analyte's own
   detection rate.** The implied improvement lever is the conditional *model*
   (a censored/Tobit regression that accounts for the truncated training
   sample), not more auxiliary analytes.

### Full-grid results (72 scenarios, 50 reps, M=10; n in {100,300}, p in {5,15})

The full grid (`run_full.R`) confirms the boundary and adds two things the
fixed-`n` sweep could not show:

| Censored | gsimp bias | MI coverage | Reliable? |
|---|---|---|---|
| 25% (15% ND) | ~0.04 | 1.00 (all n, p, rho, skew) | yes |
| 45% (35% ND), symmetric | ~0.15 | 0.57 @ n=100 -> 0.01 @ n=300 | no |
| 45% (35% ND), right-skew | ~0.01-0.04 | ~1.00 | yes |
| 65% (55% ND) | 0.35 (0.22 skewed) | ~0.00 | no |

3. **Coverage failure worsens with sample size.** At 45% censoring the point
   bias is unchanged across n (~0.15), but coverage falls from ~0.5 (n=100) to
   ~0.01 (n=300): the bias is systematic, so a larger n only shrinks the CI
   around the wrong centre. **More samples cannot fix censoring bias -- they make
   the inference worse.** The reliable detection-rate threshold is therefore
   stricter for larger studies.

4. **Distribution shape matters.** Right-skew (sinh-arcsinh) sharply reduced bias
   and restored coverage at 45% censoring (bias 0.15 -> ~0.02), so the standard
   is worst for symmetric log-scale data and more forgiving for right-skewed
   data. `p` (5 vs 15) is a minor effect; gsimp ~halves naive-substitution bias
   throughout.

### The tobit conditional model fixes it (`tobit_vs_ridge.R`)

The observed-only ridge model fails above ~25-35% censoring because it is fit on
the upper-truncated *detected* values. The `imp_model = "tobit"` censored
Gaussian model, fit on observed **and** censored rows, removes that bias
entirely. Paired comparison (same data, detection sweep, n in {100, 300}):

| Detection rate | ridge bias / coverage | tobit bias / coverage |
|---|---|---|
| 85% | 0.03 / 1.00 | ~0.00 / 1.00 |
| 75% | 0.07 / 0.96-1.00 | ~0.00 / 1.00 |
| 65% | 0.13 / 0.04-0.71 | ~0.00 / 1.00 |
| 55% | 0.21 / ~0.05 | ~0.00 / 1.00 |
| 45% | 0.31 / ~0.01 | ~-0.01 / ~1.00 |

**Lowest reliable detection rate: ridge 75%, tobit at least 45%** (the sweep
floor -- tobit may extend further), at *both* sample sizes. tobit also removes
the sample-size penalty (it works at n=300 where ridge fails hardest). Its
coverage sits at ~1.00, i.e. mildly conservative rather than under-covering --
the safe direction. This makes the censored conditional model the clear path to
extending reliable imputation into heavy-censoring regimes.

**Heterogeneous limits, revisited under tobit (`heterogeneous_limits.R`,
MODELS="ridge,tobit").** For ridge, well-observed anchor analytes did not rescue
heavily-censored (40% ND) targets at any correlation. tobit makes those targets
reliable (bias ~0, coverage 1.00) at **every** correlation *including rho = 0* --
so the rescue is not via the anchors but via tobit's per-analyte censored
likelihood. Correlation is expected to affect tobit's *efficiency* (interval
width), not its bias; that is the open follow-up.

**Effect of the number of variables (`effect_of_p.R`, at 40% ND).** tobit meets
the ND<50% target at every `p` from 3 to 200 (`p/n` up to 2, into the PCA
regime), while ridge fails at all `p` and worsens as `p` grows. Two nuances for
tobit: a mild *negative* bias drift as `p >> n` (-0.02 -> -0.045, from PCA signal
loss -- within target but eroding, so `max_pc` is a tunable knob), and RMSE that
is U-shaped in `p`, minimised around `p = 12-50`. Median bias ~0 throughout.

These are the substantive results for a methods write-up. Grids used moderate
replication (transition-zone coverage carries Monte Carlo SE ~+/-0.03-0.05). For
the definitive run use higher replication on capable hardware:

```bash
CONFIG=full N_REP=500 M=50 Rscript validation/mc_validation.R
N_REP=200 M=25 Rscript validation/tobit_vs_ridge.R
```
