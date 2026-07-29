# Design Note: `leftcens` — Left-Censored / Interval-Censored Descriptive Stats and Imputation Package

**Status:** Naming decided — package name is **`leftcens`**.

### Naming history

Candidates considered, roughly in order: `limitstats`, `mdlstats`, `gsimpIC`,
`icimpute`, `icdstats`, `censIC` (interval-centric framing) → `detectstats`,
`censdetect`, `cenStats`/`censtats`, `censimp` (broadened once interval censoring was
recognized as one input convention rather than the organizing concept) →
`lftcenstats`, **`leftcens`** (final — anchors the name on left-censoring, the actual
statistical problem, rather than on "interval" or "stats," which each only described
one part of the package).

`ic_data`/`ic_np_fit`/`ic_sp_fit`/`gsimp_ic()` in earlier drafts of this note are
renamed below to `cens_data`/`cens_np_fit`/`cens_sp_fit`/`gsimp_impute()` to match —
this resolves the open naming question from the previous draft (Section 8).

Before first `library()`/`devtools::load_all()` use, confirm with
`available::available("leftcens")` that there's no CRAN/Bioconductor collision or
other conflict.

---

## 1. Purpose and scope

A single R package providing tools for analytical/environmental measurement data with a
three-tier interval-censoring structure:

| Category | Interval coding |
|---|---|
| Non-detect | `(0, MDL)` |
| Detected, not quantified | `(MDL, LCMRL)` |
| Quantified | `(value, value)` |

Two capability modules, sharing one data representation:

1. **Descriptive statistics** — non-parametric (Turnbull NPMLE) and semi-parametric
   (proportional-hazards/odds baseline) summaries of interval-censored data.
2. **Imputation** — Gibbs-sampling-based imputation (GSimp lineage) extended to accept
   per-observation `(lo, hi)` bounds rather than a single per-column bound, so it can
   impute the non-detect and detected-not-quantified bands correctly within the same
   variable.

**Explicitly out of scope (for now):**
- Right-censored / above-quantitation-limit data.
- General n-way interval schemes beyond the three categories above.
- Model fitting (regression / brms) — that stays in `generic-mi-brms-pipeline`, which
  consumes this package's imputed output.

---

## 2. Provenance

- Imputation core is carried over from `NDim4jecs` (itself a package wrapper around
  Rum Wei et al.'s **GSimp**, https://github.com/WandeRum/GSimp,
  doi:10.1371/journal.pcbi.1005973), specifically `multi_impute()` / `GS_impute()` and
  the truncated-normal Gibbs draw (`rnorm_trunc`).
- Descriptive-stats module is new, built on `survival::survfit(type="interval2")` /
  `icenReg` rather than on NADA/NADA2, since those assume single-limit left censoring.

### 2.1 Licensing — needs resolution before any code is copied

This is the one item that should be settled **before** carrying code over, not after:

- The original GSimp scripts are licensed **CC BY-NC-SA 4.0** (non-commercial,
  share-alike).
- `NDim4jecs`'s own repository license is listed as **CC BY 4.0**, which does not
  obviously reconcile with wrapping NC-SA-licensed source underneath it.
- CC BY-NC-SA is **not** compatible with the standard R package licenses
  (MIT/GPL/CC0) typically expected for CRAN, and the NC clause restricts commercial
  use/redistribution outright regardless of what license wraps it.

Action items:
- Confirm with the original GSimp author whether relicensing or explicit permission is
  available, **or**
- Treat the ported Gibbs-sampling code as CC BY-NC-SA-encumbered, keep the new package
  itself unpublished/internal-only (or clearly marked NC), and do **not** submit to
  CRAN under a conflicting license, **or**
- Reimplement the truncated-normal Gibbs step independently from the published
  algorithm description (PLOS Comp Biol paper) rather than copying source, which
  sidesteps the code-license question (the statistical method itself is not
  copyrightable, only the specific implementation).

This determines how the package can legally be distributed, so it should be the first
open question closed out, not a footnote.

---

## 3. Shared data representation

Single constructor/validator used by both modules, so the interval semantics are
defined exactly once:

```r
as_interval_data(left, right, log_transform = TRUE)
```

Responsibilities:
- Validate `left <= right` elementwise (error otherwise).
- Classify each row into `non_detect` / `detected_not_quantified` / `quantified`
  (derived, not user-supplied, to avoid the two representations drifting out of sync).
- Handle `left = 0` correctly under `log_transform = TRUE` (maps to `-Inf` on log
  scale — this is the expected non-detect lower bound, not an error).
- Return a lightweight S3 object (class `cens_data`) with `$left`, `$right`,
  `$category`, `$log_transform`, so downstream functions can dispatch on it and print/
  summarize it consistently.

Both the descriptive-stats functions and the imputation bound-matrix builder consume
`cens_data` objects — this is the thing that justifies keeping both modules in one
package rather than splitting them.

### 3.1 Single (Helsel-style) left-censoring is a degenerate case, not a second mode

Classic single-limit left-censoring — one detection limit per observation, no
detected-not-quantified band (the format NADA/NADA2 work with) — is just the two-tier
collapse of the same interval representation: `(0, DL)` for non-detect and
`(value, value)` for detected. Nothing downstream needs to know whether the middle
tier is present or absent; `desc_np()`, `desc_sp()`, and `gsimp_impute()` all operate on
`cens_data` objects regardless of how many distinct interval shapes occur in the data.

Consequence: the package should target **both** conventions through one shared data
model, not through two separate modes or two packages. The only thing that needs
adding is a thin, ergonomic constructor so users coming from single-limit workflows
don't have to hand-expand `(value, censored)` vectors into `(left, right)` pairs
themselves:

```r
as_interval_data_from_single(value, censored, log_transform = TRUE)
# censored == TRUE  -> left = 0,     right = value
# censored == FALSE -> left = value, right = value
# internally just calls as_interval_data(); no separate logic path
```

A similar adapter could accept a `NADA::Cen()`-style object directly, for users
migrating existing single-limit datasets, again as a wrapper around
`as_interval_data()` rather than new machinery.

---

## 4. Package architecture

```
R/
  core-cens.R          # as_interval_data(), as_interval_data_from_single(),
                        # validation, category classification -> cens_data objects
  desc-nonparam.R      # Turnbull NPMLE wrappers (survival::survfit / icenReg::ic_np),
                        # quantile extraction, plotting
  desc-semiparam.R     # icenReg::ic_sp() wrappers (baseline PH/PO curves)
  impute-bounds.R      # build lo_mat / hi_mat (per-cell bounds) from a cens_data object
  impute-gsimp.R       # gsimp_impute() (patched multi_impute()/GS_impute()),
                        # matrix-valued lo/hi
  utils.R              # shared helpers (log/exp round-trip, printing, etc.)
tests/
  testthat/
    test-core-cens.R
    test-desc-nonparam.R
    test-desc-semiparam.R
    test-impute-bounds.R
    test-impute-gsimp.R      # includes a recovery test against simulated ground truth
vignettes/
  getting-started.Rmd
  imputation-bounds.Rmd      # documents the per-cell lo/hi patch and its assumptions
```

### 4.1 Dependency strategy

| Package | Role | `Imports` or `Suggests` |
|---|---|---|
| `survival` | Turnbull NPMLE | Imports (light, base-adjacent) |
| `icenReg` | NPMLE + semi-parametric IC regression | Imports |
| `missForest`, `imputeLCMD`, `glmnet` | GSimp prediction models / init | Suggests |
| `doParallel`, `foreach` | Parallel Gibbs iterations | Suggests |
| `impute`, `pcaMethods`, `ropls` (Bioconductor) | GSimp comparison wrappers | Suggests |

Imputation functions should `requireNamespace(..., quietly = TRUE)` and fail with an
informative message (not a silent `NA`) if a `Suggests`-only package is missing. This
keeps `library(leftcens)` cheap for a user who only wants descriptive statistics.

---

## 5. Function API sketch

```r
## --- shared ---
as_interval_data(left, right, log_transform = TRUE)                   # -> cens_data
as_interval_data_from_single(value, censored, log_transform = TRUE)   # convenience wrapper
as_interval_data_from_cen(x, log_transform = TRUE)                     # optional: NADA::Cen() adapter

## --- descriptive: non-parametric ---
desc_np(x, method = c("turnbull", "wang"))       # -> cens_np_fit
quantile.cens_np_fit(object, probs = c(.1,.25,.5,.75,.9))
plot.cens_np_fit(x, ...)

## --- descriptive: semi-parametric ---
desc_sp(x, covariates = NULL, model = c("ph", "po"))  # -> cens_sp_fit
plot.cens_sp_fit(x, ...)

## --- imputation ---
build_bounds(x)                     # cens_data -> list(lo_mat, hi_mat)
gsimp_impute(data_wide, bounds, iters_each = 100, iters_all = 20,
             initial = "qrilc", imp_model = "glmnet_pred", n_cores = 1)
```

`gsimp_impute()` is the patched `multi_impute()`: accepts `lo`/`hi` as matrices matching
`dim(data_wide)` (one bound pair per cell, `NA` for non-missing cells), and slices
`lo_mat[is.na(y_miss), j]` / `hi_mat[is.na(y_miss), j]` per column before calling the
existing (unmodified) `rnorm_trunc()` draw, which already vectorizes over per-row
bounds. This is the minimal patch identified against the upstream `GSimp.R` source —
no change needed to the Gibbs sampler itself, only to how bounds are constructed and
passed down from `multi_impute()`.

---

## 6. Testing / validation plan

- **Unit tests** for `as_interval_data()`: boundary cases (`left = 0`, `left = right`,
  invalid `left > right`), correct category classification.
- **Recovery test** for `gsimp_impute()`: simulate a complete dataset, artificially censor
  it into the three-tier interval structure with known bounds, impute, and check that
  imputed values fall within their assigned `(lo, hi)` and that aggregate statistics
  (mean/sd/correlation structure) recover the uncensored ground truth within tolerance.
  This is the equivalent of GSimp's own MNAR simulation/evaluation harness, adapted to
  the interval case.
- **Cross-check** `desc_np()` Turnbull output against `survival::survfit()` directly on
  a small hand-computable example, to catch wrapper bugs independent of the underlying
  algorithm.
- **Regression check**: verify `gsimp_impute()` reduces to the original `multi_impute()`
  behavior when `lo`/`hi` are supplied as scalars/vectors (backward compatibility with
  the single-limit, single-bound use case).

---

## 7. Roadmap

1. Resolve licensing (Section 2.1) — blocking.
2. Set up the `leftcens` package skeleton (`usethis::create_package()`,
   `DESCRIPTION`, `NAMESPACE`, roxygen scaffolding) and implement `core-cens.R` + tests.
3. Implement `desc-nonparam.R` (Turnbull wrapper) — lowest-risk, no license entanglement.
4. Implement `desc-semiparam.R`.
5. Implement `impute-bounds.R` + `impute-gsimp.R`, with the recovery-test
   simulation as the acceptance criterion before this is used on real study data.
6. Vignettes; run `available::available("leftcens")` before any public release/CRAN
   submission if not already done in Section 0.
7. **Keep `leftcens` independent of `generic-mi-brms-pipeline`** (decision
   2026-07-29 — supersedes the earlier plan to add a `gsimp_impute()` strategy
   inside the pipeline). Rationale: the pipeline's imputer (`miceRanger`) assumes
   **MAR** — it treats missingness as uninformative about magnitude. Left-censored
   data is **MNAR** (a value is missing precisely because it falls below a
   detection limit), so routing it through a MAR imputer biases estimates upward
   (the naive-substitution bias). The pipeline's own docs already say to handle
   such data appropriately upstream; `leftcens` *is* that upstream step
   (`MNAR-censored data → completed data → pipeline`). Independence honours the
   boundary the pipeline declares, avoids forcing a censoring representation into
   a mature MAR-oriented pipeline, and avoids a local non-CRAN dependency.
   Investigated and confirmed: the pipeline has no censoring representation today.

   If interoperability is ever wanted, the clean seam — *without* coupling the
   packages — is `leftcens` emitting **M** completed datasets that the pipeline
   consumes as its M imputations (pipeline imputation set to pass-through), so the
   censoring-imputation uncertainty propagates through the pipeline's Rubin
   pooling, rather than feeding a single completed dataset as if it were ground
   truth.

---

## 8. Open questions

- Does the semi-parametric module need to support repeated-measures / clustered data
  (mirroring the `generic-mi-brms-pipeline` repeated-outcome support), or is a single
  cross-sectional interval variable sufficient for now?
- Should `gsimp_impute()`'s per-cell bounds also allow *analyte-and-sample*-varying
  MDL/LCMRL (i.e., a full matrix supplied by the user), or is
  column-constant-with-row-exceptions sufficient for the current data? (The patch as
  designed supports the fully general case at no extra cost, so this is a
  documentation/API-surface question, not an implementation one.)
