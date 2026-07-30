# Research plan: sensitivity analyses, guidance, and a pre-flight tool

Planning doc for the next phase of `leftcens` validation work. Committed so it
travels across machines. Not part of the package build (`validation/` is
Rbuildignored).

## ▶ Start here (resuming on another machine)

1. `git pull`, then install so the scripts' `library(leftcens)` resolves:
   `R CMD INSTALL .` (or `devtools::install()`).
2. **E1, E2, E3 are DONE** (see their result blocks below). **Next task: E5 —
   the three-tier DNQ band** (the untested gap and the package's raison d'etre;
   E4's correlation question is largely answered by E1+E3, so it drops in
   priority). Then E4 (efficiency/interval-width), E6, E7.
3. Then follow **Recommended order** at the bottom; every driver is
   env-configurable (`N_REP`, `M`, `ITERS_ALL`, grid params) for scaling up.
4. Package is at **0.2.0** (tobit is the default `imp_model`); tests green,
   `R CMD check` 0/0/0. Heavy runs go here on the workstation — the 10-min /
   library-path limits from the Claude Code environment do not apply.

## Target (north star)

**Rescue analytes with ND < 50%.** Rationale: below 50% non-detects the median
sits at or above the detection limit, so a recoverable center always exists --
that makes ND < 50% the natural line for "worth rescuing." Success = the default
(tobit) meets the reliability criteria (|mean bias| <= 0.10 AND MI coverage
>= 0.90, with the median recovered) across the whole 0 - <50% ND range, under
realistic conditions. The 40-49% band is the hard edge; ~50% is a fragile
boundary (the median lands on the limit), so the target is strictly < 50%.

Every experiment below is measured against this target: the question is not just
"what happens" but "does rescue still hold up to ~49% ND, and if not, under which
conditions does it fail first?" Report the **median** (and other quantiles), not
just the mean --- the median is the user's motivating estimand.

**Baseline boundary check** (n=200, p=6, rho=0.5, symmetric, pure left-censoring,
25 reps, M=12, tobit):

| ND  | mean bias | median bias | MI coverage |
|-----|-----------|-------------|-------------|
| 40% | +0.006    | 0.000       | 1.00        |
| 45% | -0.007    | 0.000       | 1.00        |
| 49% | +0.004    | 0.000       | 0.98        |

=> Target looks **achievable at baseline**. E2-E5 test whether it survives more
variables, non-parametric correlation, skew, the DNQ tier, and heterogeneous
limits.

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

### E2 — Number of variables k (= p) imputed together  (point 2)  [DONE]
- **Question:** how does recovery behave as `p` grows and crosses `p ≈ n`, where
  the PCA reduction engages? Does it degrade gracefully at `p >> n`?
- **Result** (`effect_of_p.R`; n=100, rho=0.5, 40% ND, 12 reps, p=3..200):
  **the rescue is robust to variable count — tobit meets the target
  (|mean bias| <= 0.10) at every p, including p=200 (p/n=2, PCA regime)**, while
  ridge fails at all p and worsens with p (0.16 -> 0.28). Median bias ~0
  throughout (< 50% ND). Two nuances: (a) tobit shows a mild *negative* bias
  drift as p >> n (-0.02 at p=3 -> -0.045 at p=200) from PCA signal loss --
  within target but eroding; (b) **RMSE is U-shaped in p, minimised ~p=12-50**
  (0.014) -- too few analytes = weak model, too many = PCA truncation.
- **Follow-ups:** tune the PCA `max_pc` cap (currently 20) to reduce the wide-p
  overcorrection; add MI coverage at a few large-p points (this sweep was
  point-metrics only for speed).
- **Run it (higher rep):** `REPS=50 ITERS=20 Rscript validation/effect_of_p.R`

### E3 — Correlation: parametric vs non-parametric dependence  (point 3)  [DONE]
- **Question:** tobit's conditional model is *linear*. Does non-linear
  dependence break it, where linear-Gaussian dependence does not?
- **Generator:** `simulate_truth_nl()` in `mc_validation.R` — a shared latent
  through orthogonal non-linear bases (Hermite H2/H3, sine, abs), giving strong
  dependence with low linear correlation (verified: linear |corr| ~0.29 vs
  ~0.60 for the Gaussian generator, while 2nd-moment dependence is higher).
- **Result** (`effect_of_dependence.R`; n=150, p=6, 40% ND, 30 reps): **tobit
  meets the target under BOTH regimes** (linear: bias -0.005, RMSE 0.026;
  non-linear: bias -0.036, RMSE 0.044); ridge fails both. Non-linear dependence
  does **not** break the rescue — it costs *efficiency* (RMSE 0.026 -> 0.044)
  and a mild bias drift, not validity. Confirms tobit rescues via the
  per-analyte censored likelihood (E1), robust to dependence *shape*. Median
  bias 0 throughout.
- **Follow-up / extension:** a **non-linear censored conditional model**
  (spline/GAM or forest-based, interval-censored) could reclaim the efficiency
  lost under non-linear dependence — a Phase-3 methodological extension.
- **Run it (higher rep):** `REPS=100 ITERS=25 Rscript validation/effect_of_dependence.R`

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
