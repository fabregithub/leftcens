# Research plan: sensitivity analyses, guidance, and a pre-flight tool

Planning doc for the next phase of `leftcens` validation work. Committed so it
travels across machines. Not part of the package build (`validation/` is
Rbuildignored).

## ▶ Start here (resuming on another machine)

1. `git pull`, then install so the scripts' `library(leftcens)` resolves:
   `R CMD INSTALL .` (or `devtools::install()`).
2. **E1 is DONE** (tobit arm added to `heterogeneous_limits.R`; see the E1
   result below). **Next task: E2 — number of variables / PCA regime.**
3. Then follow **Recommended order** at the bottom; every driver is
   env-configurable (`N_REP`, `M`, `ITERS_ALL`, grid params) for scaling up.
4. Package is at **0.2.0** (tobit is the default `imp_model`); tests green,
   `R CMD check` 0/0/0. Heavy runs go here on the workstation — the 10-min /
   library-path limits from the Claude Code environment do not apply.

## Where things stand (start-of-plan snapshot)

- **Package: leftcens 0.2.0.** `gsimp_impute()` defaults to `imp_model = "tobit"`
  (interval-censored Gaussian regression via `survival::survreg`, fit on
  observed *and* censored rows). Wide data (`p >= n`) handled by PCA reduction
  (cap `k = n/3`). Alternatives: `ridge` (old default), `lm`, `glmnet`, or a
  user `function(yo, Xo, Xm)`.
- **Harness (`validation/`):** `mc_validation.R` (engine: `simulate_truth()`,
  `censor_three_tier()` [per-analyte fractions], `run_validation(imp_model=)`,
  `mi_coverage()`, `summarise_validation()`), `detection_rate_standard.R`,
  `heterogeneous_limits.R`, `tobit_vs_ridge.R`, `run_full.R` (batch/accumulate),
  `report.Rmd`.
- **Findings so far** (see `validation/README.md`):
  - Observed-only (ridge): reliable to ~25% non-detects; **coverage fails before
    bias**; **larger n makes inference worse**; **correlation does not help**;
    **heterogeneous anchors do not rescue** heavily-censored targets;
    right-skew is more forgiving. Common cause: selection bias from fitting on
    upper-truncated observed values.
  - **tobit** (new default) removes that bias: near-zero bias and near-nominal
    coverage down to at least ~55% non-detects, at both small and large n; PCA
    path extends it to wide data.

## Framing

Most of the sensitivity questions below were **already answered for ridge**, and
those answers were driven by ridge's selection bias. The scientific narrative for
this phase is **how the censored (tobit) model changes each conclusion.** Prefer
paired ridge-vs-tobit designs (same data, only the model differs), as in
`tobit_vs_ridge.R`.

## Experiments (ranked by novelty / value)

### E1 — Heterogeneous detection limits under tobit  (was point 5)  [DONE]
- **Question:** do well-observed "anchor" analytes now rescue heavily-censored
  "target" analytes, given tobit removes the selection bias that blocked this
  for ridge?
- **Prior (ridge):** no rescue even at rho = 0.8 (`heterogeneous_limits.R`).
- **Result (tobit arm, MODELS="ridge,tobit"; 40% ND targets, 5% ND anchors, 20
  reps, M=10):** tobit makes the target analytes reliable (bias ~0, coverage
  1.00) at **every** rho, *including rho = 0*, while ridge stays biased
  (~+0.16) with coverage ~0.05-0.12 throughout. So the rescue is **not** via the
  anchors/correlation — it is tobit's *per-analyte* censored likelihood removing
  the selection bias directly (even an intercept-only censored fit recovers the
  mean). Overturns the ridge-era "targets unreliable" conclusion, but the
  mechanism is "the censored model makes correlation unnecessary", not
  "correlation now helps".
- **Refinement this opened:** correlation may improve tobit *efficiency* (tighter
  intervals — its coverage sits at a conservative 1.00) even though it does not
  affect *bias*. Test with an RMSE / interval-width metric in E3/E4.
- **Run it (higher rep):** `N_REP=200 M=25 Rscript validation/heterogeneous_limits.R`

### E2 — Number of variables k (= p) imputed together  (point 2)  [NEW]
- **Question:** how does recovery behave as `p` grows and crosses `p ≈ n`, where
  the PCA reduction engages? Does it degrade gracefully at `p >> n`? Do more
  (informative) analytes help tobit (unlike ridge)?
- **Prior:** only `p in {5,15}` tested, ridge only (minor effect).
- **How:** grid over `p` (e.g. 3, 6, 12, 25, 50, 100, 200) at fixed `n`, both
  models; watch the narrow -> PCA transition and the PCA `max_pc`/`k` behaviour.
  New driver, reuses the engine.

### E3 — Correlation: parametric vs non-parametric dependence  (point 3)  [NOVEL ROBUSTNESS]
- **Question:** tobit's conditional model is *linear*. Does non-linear / copula
  dependence break it, where linear-Gaussian dependence does not?
- **Prior:** all sims used Gaussian (linear, exchangeable) correlation.
- **How:** needs **new DGP code** — a copula generator and/or a non-linear
  dependence generator, added to `mc_validation.R`. Compare tobit under linear
  vs non-linear dependence.
- **If tobit degrades under non-linearity:** motivates a future non-linear
  censored conditional model (spline/GAM or forest-based) — a Phase-3 extension.

### E4 — Censoring rate x correlation interaction under tobit  (point 4)
- **Question:** does correlation help tobit at high censoring (where it was
  inert for ridge)? Per E1, expect correlation to affect *efficiency*, not
  *bias*, for tobit.
- **How:** extend `detection_rate_standard.R` to run a tobit arm across the
  rho x detection-rate grid. **Record RMSE and mean MI interval width**, not
  just bias/coverage, since for tobit the interesting effect of correlation is
  narrower intervals / lower variance rather than reduced bias. Partly overlaps
  E3.

### Gaps to fold in
- **E5 — Three-tier (DNQ) band.** Almost everything used pure left-censoring
  (`dnq_frac = 0`). The middle tier is the package's raison d'etre; give it a
  sensitivity slice (`dnq_frac > 0`).
- **E6 — Convergence guidance.** How large must `M` (imputations) and
  `iters_all` (sweeps) be for stable bias/coverage? Cheap; needed for E-guidance.
- **E7 — Runtime / scalability.** Profile tobit vs ridge vs `p`/`n` (tobit is
  ~8x ridge; PCA adds cost). Feeds the guidance's compute recommendations.

## Deliverables

### D1 — Definitive high-replication run  (point 1)  [WORKSTATION]
Capstone validation once the axes above settle the design. On the work machine:
```bash
CONFIG=full N_REP=500 M=50 Rscript validation/mc_validation.R      # ~50 hrs
N_REP=200 M=25 Rscript validation/tobit_vs_ridge.R
```

### D2 — Guidance for the best available approach  (point 6)
Synthesise E1-E7 into a practical guide (a vignette or `validation/` doc): when
to use tobit vs ridge; expected reliability by censoring rate, correlation
strength/shape, `p`, skew, and tier structure; recommended `M` / `iters_all`;
and the compute cost. This is the user-facing decision aid.

### D3 — Pre-flight reliability tool  (point 7, reading A)
Package the harness so a user can **simulate data matching their own study
design** (their `n`, `p`, correlation structure, per-analyte censoring pattern,
skew) and get a reliability read (bias / MI coverage) **before** trusting
imputation on the real data. E.g. a `leftcens` function like
`preflight_reliability(design_spec, ...)` or a templated validation vignette the
user fills in with their design. Turns the internal study into a
"test-it-on-your-design-first" workflow.

## Division of labour

- **In Claude Code (this repo):** build the experiment drivers (E1-E7 are mostly
  new grids + drivers reusing the engine; E3 needs new DGP code) and run modest
  foreground confirmations to check each works and shows the expected direction.
- **On the workstation:** the heavy / definitive runs (D1, and high-replication
  versions of E1-E7). All drivers are env-configurable (`N_REP`, `M`,
  `ITERS_ALL`, grid params) for scaling up.

## Recommended order

E1 -> E2 -> E3 -> E4 -> (E5, E6, E7) -> D1 (definitive) -> D2 (guidance)
-> D3 (pre-flight tool).

## Quick start when resuming (work machine)

```bash
# from the repo root, install so validation scripts' library(leftcens) resolves
R CMD INSTALL .            # or devtools::install()

# run any driver (env-configurable), e.g. the ridge-vs-tobit comparison:
N_REP=200 M=25 Rscript validation/tobit_vs_ridge.R

# render the report from the latest results:
Rscript -e 'rmarkdown::render("validation/report.Rmd")'
```

Notes: install `leftcens` into a library on your `.libPaths()` (a plain
`R CMD INSTALL .` does this). The 10-minute foreground cap and the detached-R
library-path quirks seen inside Claude Code do **not** apply on a normal
workstation — run the heavy jobs directly there.
