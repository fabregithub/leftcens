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

- `mc_validation.R` — simulation, censoring, metrics, runner, summariser, and a
  `main` block that runs a config and saves results. Source it to reuse the
  functions (`run_validation()`, `summarise_validation()`, ...).
- `report.Rmd` — reads `results/latest.rds` and produces tables and figures.
- `results/` — generated outputs (git-ignored).

## What to expect

gsimp reliably reduces the upward bias of naive substitution. The limiting
regime is heavy censoring: as the censored fraction rises, residual bias
persists and MI coverage degrades (intervals correctly narrow but centred
slightly high). That boundary — where the method helps most and where it starts
to strain — is the substantive result to report.
