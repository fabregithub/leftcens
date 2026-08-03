# Phase 3 plan: the skewed-mean problem — diagnosis, candidate solutions, and a test program

Planning doc for the next methodological phase of `leftcens`, triggered by the
D1 capstone finding (2026-08-02): under substantial right-skew the default
tobit imputation gives a **downward-biased log-scale mean** that breaches the
|bias| ≤ 0.10 target by ~35% ND, even though the **median is recovered exactly**
across the whole ND < 50% range. This doc works the mathematics, enumerates
every solution family, and lays out experiments to decide among them. Committed
so it travels across machines (like `RESEARCH_PLAN.md`; `validation/` is
Rbuildignored).

Two questions to answer, from the brief:
1. **Is the Gibbs sampler the right engine to stick with?**
2. **What if the right-skewed variables to be imputed are correlated?**

---

## 1. Diagnosis — where the bias actually comes from (the maths)

### 1.1 The estimand is tail-sensitive under left-censoring
The target is the marginal log-scale mean `θ_j = E[X_j]`. Split it at the
detection limit `c_j = log(MDL_j)`:

```
E[X_j] = E[X_j | X_j ≥ c_j]·P(X_j ≥ c_j)  +  E[X_j | X_j < c_j]·P(X_j < c_j)
                 (observed part)                     (imputed part)
```

The second term is **entirely model-extrapolated** — no data live below `c_j`,
only the count. So the mean's accuracy is exactly the accuracy of the modelled
left-tail expectation, weighted by the non-detect fraction `π_j = P(X_j < c_j)`.
By contrast the **median** for `π_j < 0.5` sits at/above `c_j`, so it is a
functional of the *observed* part only — which is why D1 shows median bias ≈ 0
at every skew. **The mean fails where the median cannot: in the unobserved tail.**

### 1.2 The misspecification-bias formula
The current conditional (`tobit_predict` → `survreg` Gaussian) draws each
non-detect from a **homoscedastic truncated normal** `N(μ_ij, σ)` on `(−∞, c_j)`,
with a single scale `σ = fit$scale` (see `R/impute-gsimp.R:412`). The resulting
marginal-mean bias is, to first order,

```
bias(θ_j) ≈ π_j · ( E_model[X_j | X_j < c_j] − E_true[X_j | X_j < c_j] )
          = π_j · Δ_tail
```

`Δ_tail` is the gap between the modelled and true left-tail means. For a Gaussian
fit whose σ is inflated by a long *right* tail, the modelled left tail is too
heavy → `Δ_tail < 0` → **downward bias, scaling linearly with π_j.** This predicts
the D1 pattern (−0.07 → −0.21 → −0.38 as ND = 15 → 35 → 55%).

**Numerically verified** (sinh-arcsinh skew δ = 0.75, matches the DGP): skew(Y) =
+1.0; below the 35% MDL, `E_true = −0.14` but the moment-matched Gaussian gives
`E_model = −0.46`, so `Δ_tail = −0.32`; `π·Δ_tail = 0.35·(−0.32) = −0.11`. The
full-model bias (−0.21) is larger because the *conditional* linear predictor and
inflated conditional scale add to the marginal effect. Sign, scaling, and order
of magnitude all confirm the mechanism.

### 1.3 Why sinh-arcsinh skew does this
For `Y = sinh(asinh(Z) + δ)`, `Z ~ N(0,1)`, `δ > 0`: the right tail is stretched
by `≈ e^{+δ}` and the left tail **compressed** by `≈ e^{−δ}` (from
`sinh(asinh(z)+δ) ≈ z·e^{±δ}` as `z → ±∞`). So the true left tail is *thin*
(verified sd = 0.36 vs the Gaussian's 1.34). Left-censoring samples precisely
this thin tail, while σ is estimated from the whole (right-tail-inflated)
distribution. Any generator with a lighter-than-Gaussian left tail (log-gamma,
generalized-gamma, GEV) will do the same; sinh-arcsinh is one representative.

### 1.4 What is and is not identified
- **Observed-range shape** (skew above `c_j`) is identified and can be estimated.
- **Below-`c_j` marginal tail** is *not* nonparametrically identified from analyte
  `j` alone — we only know its mass `π_j`. Recovering `E_true[X_j | X_j<c_j]`
  therefore requires either (a) a parametric/semiparametric tail **assumption**,
  or (b) **information borrowed from correlated analytes** (see §3). This is the
  crux of the whole problem, and it is exactly where correlation earns its keep.

### 1.5 Is Gibbs the culprit? (answering question 1)
**No.** Separate the two layers:
- **Engine:** Gibbs / fully-conditional-specification (FCS) — a scheme for
  sampling the joint from per-analyte conditionals.
- **Model:** the linear-Gaussian interval-censored conditional.

The bias is 100% in the *model*: Gibbs is faithfully sampling the wrong
conditional. Moreover the FCS conditionals here (each analyte linear-Gaussian on
the rest) are **compatible with a joint multivariate normal**, so on a correctly
specified latent (see S1) Gibbs converges to the right joint — E6 already showed
fast mixing (`iters_all ~ 3–5`). **Conclusion: keep the Gibbs/FCS engine; fix the
conditional distribution.** We will still benchmark a non-Gibbs joint-EM
alternative (S3) to *prove* the engine is not the limiting factor, but the prior
is strongly that the engine stays.

---

## 2. The key constructive result — Gaussian-copula imputation

Model the joint as a **Gaussian copula with arbitrary margins**:
`X_j = F_j^{-1}(Φ(Z_j))`, with `Z = (Z_1,…,Z_p) ~ N(0, Σ)`. Then:

1. **Skew lives entirely in the margins** `F_j`; the dependence lives in `Σ`.
2. Map observed values to latent normal scores `Z_j = Φ^{-1}(F_j(X_j))`.
   Left-censoring `X_j < c_j` becomes a **latent threshold** `Z_j < Φ^{-1}(F_j(c_j))`.
3. In latent space the conditional is **exactly linear-Gaussian and
   homoscedastic**:
   ```
   Z_j | Z_{-j} ~ N( Σ_{j,-j} Σ_{-j,-j}^{-1} Z_{-j},  1 − Σ_{j,-j}Σ_{-j,-j}^{-1}Σ_{-j,j} )
   ```
   — precisely the assumption `survreg` + `rnorm_trunc` already encode. So the
   **existing machinery becomes correctly specified** once we impute in `Z`-space.
4. Back-transform `X_j = F_j^{-1}(Φ(Z_j))` and average. The mean is now unbiased
   **given `F_j`**, so the whole problem collapses onto *marginal tail estimation*
   (§1.4) — a 1-D problem per analyte, not a joint one.

**Correlation and question 2 (hypothesis — see G3 for the verdict).** The latent
conditional variance `1 − Σ_{j,-j}Σ_{-j,-j}^{-1}Σ_{-j,j}` shrinks as correlation
grows, so *per-cell* a censored `Z_j` is pinned into a narrower band by observed
correlated analytes. This suggested correlation might now **help the marginal
mean** under the copula (the predicted "flip").

> **G3 REFUTED this for the marginal mean.** Correlation does NOT reduce the
> copula's mean bias/RMSE — RMSE mildly *rises* with rho, because (i) the marginal
> mean's variance is a *marginal* quantity that per-cell conditional shrinkage
> does not lower, and (ii) *concordant* censoring at high rho means correlated
> analytes are censored together → *less independent* information. The pinning
> intuition holds only when the neighbour is *observed*, which fails under heavy,
> concordant censoring. Correlation's real payoff is the **joint** estimand
> (recovering the dependence structure), not the mean — confirming E4. So the
> copula's rescue is the censored likelihood + transform (per §1), NOT correlation.

Remaining risk: estimating `F_j` **below `c_j`**. Options, in increasing
strength/assumption: empirical CDF on the observed part + (i) a parametric tail
(skew-normal / sinh-arcsinh / generalized-gamma) fit by censored MLE, or (ii)
copula-borrowed extrapolation, or (iii) a monotone semiparametric transform
(e.g. via `icenReg`'s interval-censored AFT to Gaussianise). Sensitivity of the
mean to this choice is itself an experiment (G5).

---

## 3. Candidate solution families (ranked by value / effort / soundness)

**S1 — Gaussian-copula / normal-scores imputation.  [PROTOTYPED 2026-08-03]**
Gaussianise each margin (censored-MLE parametric tail or semiparametric
transform), run the *existing* tobit+Gibbs in latent space, back-transform.
Fixes skew (margins) and turns correlation into tail information (copula)
simultaneously; reuses `survreg` + `rnorm_trunc` almost verbatim. New code is a
marginal-transform layer + a wrapper `imp_model = "copula"`. *Effort: low–med.*

  Implemented in `R/impute-copula.R`, reachable as `gsimp_impute(imp_model =
  "copula")`. Margin = 3-parameter sinh-arcsinh (μ, σ, skew) fit per analyte by
  interval-censored MLE; `eps→0` degrades gracefully to tobit. **Prototype smoke
  (n=150, p=6, ND=35%, tobit vs copula):**

  | DGP | tobit bias | **copula bias** | tobit RMSE | copula RMSE |
  |---|---|---|---|---|
  | sinh-arcsinh skew 0.75 (matched) | −0.211 | **−0.017** | 0.213 | 0.044 |
  | Gamma margins, skewness ~1 (**mismatched**) | −0.129 | **+0.004** | 0.130 | 0.029 |

  Coverage stays 0.98–1.00. Crucially it **passes the cross-family check**
  (Gamma DGP the sinh-arcsinh fit does *not* match), so the gain is not an inverse
  crime. Correlation effect (matched, skew 0.75, ND 35%): copula bias −0.027 at
  rho=0 → −0.009 at rho=0.8 (mild bias help; RMSE ~flat — the copula already near-
  solves the marginal mean at rho=0, so little headroom; G3 will map heavier ND).
  **Cost:** several× tobit (per-column MLE + tobit in latent space). *Still needs
  the full G3/G4/G5 evaluation and wide-p / heteroscedastic stress before default
  consideration; currently opt-in.*

**S2 — Skew-aware conditional model (drop-in `imp_model`).**
Replace the Gaussian conditional with a skewed interval-censored *location-scale*
regression (e.g. censored sinh-arcsinh / skew-normal, or a GAMLSS-style
μ-σ-ν fit with an interval-censored likelihood). Keeps Gibbs; changes only
`resolve_imp_model()`. Models conditional skew *and* heteroscedastic scale
directly, without a margin/copula split. *Effort: med (custom likelihood); watch
stability, speed, wide-p.*

**S3 — Parametric censored multivariate joint (skew-normal / skew-t) via EM.
[benchmark / engine test]**
Fit `X ~ MST/MSN` with left-censoring by ECM; get `E[X_j]` directly and MI draws
from the posterior predictive. Non-Gibbs — its role is to (a) provide a strong
reference and (b) test whether abandoning Gibbs buys anything. *Effort: high;
verify package support (skew-normal `sn`; censored-multivariate-skew EM exists in
the Lachos/Matos/Garay literature — availability TBD, may need custom); wide-p
scaling is the main worry.*

**S4 — Estimand/reporting + diagnostics (partial-identification hedge). [ship-now]**
Independent of S1–S3 and deployable immediately:
(a) a **skew×ND diagnostic** that predicts mean bias via §1.2 and warns;
(b) a **first-order bias correction** `θ̂ + π̂·Δ̂_tail` from an estimated tail;
(c) **Manski-style bounds** on the mean from the identified region (report an
interval, not a point). Complements the median-first guidance already shipped.

> **S4(a) DONE 2026-08-03 — skew-aware pre-flight.** `preflight_from_data()` now
> estimates each analyte's sinh-arcsinh skew by **interval-censored MLE**
> (`fit_shash_margin`, not the sample skewness of imputed data — a symmetric tobit
> fill would attenuate it) and feeds it into the simulation, so the reported
> reliability reflects the skew-induced mean bias instead of assuming log-normal.
> `preflight_reliability(skew=)` accepts a per-analyte vector; the print method
> flags material right-skew (δ ≳ 0.3) and recommends `imp_model="copula"`. Verified:
> on δ=0.75 data the tobit pre-flight now recovers δ≈0.69 and flags 0/3 reliable
> (was falsely 3/3); copula rescues 2/3. Folded into the guidance vignette.

**S0 — status quo + median-first framing.** Already done (vignette/RESEARCH_PLAN).
Baseline to beat.

*Relation to the E3 "non-linear censored conditional" extension:* S1 linearises
dependence in latent space, so it should absorb much of E3's motivation for the
*mean*; a non-linear conditional (S2/GAM) matters only if dependence is non-linear
*within* the latent copula. Test jointly in G3.

---

## 4. Test program (reuse the `validation/` harness)

All grids reuse `simulate_truth*` / `censor_three_tier` / `run_validation` with
new `imp_model` options, scored on bias, **median bias**, coverage, RMSE, CI
width, runtime.

- **G1 — Mechanism confirmation.** skew × ND; log the below-MDL extrapolation
  error and σ̂ inflation; check `bias ≈ π·Δ_tail`. Validates §1.2 (theory → data).
- **G2 — Skew sweep. [DONE 2026-08-03]** `effect_of_skew.R` (n=150, p=6, rho=0.5,
  tobit, 150 reps, M=20; 9 min parallel). skew ∈ {0,.25,.5,.75,1} × ND ∈
  {15,25,35,45}%. **Result confirms the mechanism cleanly:** mean bias is
  downward and scales as **`bias ≈ −0.85·δ·ND`** (δ = sinh-arcsinh skew; fit coef
  −0.85, R-squared essentially 1), so the tolerable skew shrinks as ND rises. The
  **crossing skew** (|bias| = 0.10), in *observable log-scale sample skewness* γ
  (δ→γ: .25→.41, .5→.75, .75→1.00, 1→1.17):

  | ND | max tolerable γ | coverage at crossing |
  |----|----|----|
  | 15% | ~1.11 | 1.00 |
  | 25% | ~0.80 | 1.00 |
  | 35% | ~0.57 | ~0.98 |
  | 45% | ~0.44 | dips below nominal |

  Practical rule: mean dependable while **γ·ND ≲ 0.18**. **Median abs bias ≈ 0
  (machine zero) in every cell** — the robust estimand throughout. Folded into the
  vignette (threshold table + rule) `imputation-guidance.Rmd`. This is the
  status-quo (S0) baseline the S1 copula must beat: S1 should push these
  crossings far higher (ideally "safe over full sweep").
- **G3 — Correlated-skew (question 2). [DONE 2026-08-03]**
  `effect_of_corr_skew.R` (parallel, 7 min; n=150 p=6 skew 0.75, rho ∈ {0,.3,.6,.9}
  × ND ∈ {25,35,45}%, tobit vs copula-hybrid, 60 reps, M=20). Two estimands:
  (a) marginal mean; (b) joint = between-analyte correlation recovery (`cor_mae`).

  **The predicted "flip" did NOT happen — an honest correction.** For the
  **marginal mean, correlation does NOT help the copula**: bias stays low and
  coverage nominal at every rho, but **RMSE and CI width mildly RISE with rho**
  (copula RMSE at ND45: 0.084 → 0.094 → 0.104 → 0.107 for rho 0 → .9; CI width
  0.90 → 0.98). Why the §2 prediction was wrong: the marginal mean's sampling
  variance is a *marginal* quantity that conditional-variance shrinkage does not
  reduce; and under *concordant* censoring at high rho, correlated analytes are
  censored together, so there is *less independent* information — the copula
  honestly reports it as slightly wider intervals. This **confirms & extends E4**
  (correlation doesn't help the marginal mean) from tobit to the copula.

  **Where correlation DOES help: the joint estimand.** `cor_mae` falls as rho
  rises for both methods (ND45 copula: 0.034 → 0.013), and the copula recovers the
  dependence at least as well as tobit. So dependence pays off for a downstream
  *joint* analysis, not the marginal mean — exactly E4's caveat.

  **Takeaway:** the copula's rescue comes from the per-margin censored likelihood
  + Gaussianising transform, NOT from correlation (mirrors the E1 mechanism). The
  method is robust across the whole rho range; "correlation helps" is true only
  for joint-structure recovery.
- **G4/G5 — Method bake-off + anti-inverse-crime. [DONE 2026-08-03]**
  `bakeoff_methods.R` (parallel, ~8 min; n=150 p=6, 4 families × 2 levels × 2 dep
  × 3 ND, 40 reps, M=12; methods naive/ridge/tobit/copula on *paired* data).
  Families break each copula assumption in turn: `sas` (matched), `gamma`+`lnorm`
  (margin broken), **`nonlin` (Gaussian-copula/dependence broken — the honest
  test)**. Realized skew/correlation recorded per cell.

  **Verdict — split bias from coverage (the combined pass-count misleads):**

  | method | mean \|bias\| | bias-pass | coverage-pass | note |
  |---|---|---|---|---|
  | naive  | 0.189 | 8%  | – | always high |
  | ridge  | 0.058 | 77% | 73% | bias-cancellation artifact at low ND; **blows up at ND 45% (0.10) & nonlin (0.11), coverage 0.68–0.75** |
  | tobit  | 0.110 | 54% | 83% | biased on skew (sas 0.178); conservative intervals |
  | **copula** | **0.015** | **100%** | 71% | **best bias in every family & ND, incl. nonlin (0.022) and ND 45% (0.023)** |

  - **copula solves the bias problem decisively** — 100% of cells |bias|≤0.10,
    across all families and censoring, INCLUDING the anti-inverse-crime `nonlin`
    (dependence it does not model). This is the Phase-3 goal met on the primary
    estimand.
  - **copula's one open gap (in the M=12, no-margin-draw prototype) was
    calibration at heavy censoring (ND≈45%):** coverage 0.835 (margin-broken) /
    0.65–0.75 (nonlin) while *bias there was ~0.02* — intervals too NARROW because
    the margin `F_j` was fit once and treated as FIXED, so MI variance omitted
    margin/tail-estimation uncertainty. **FIXED below.**
  - ridge's higher *combined* pass-count is an artifact (bias cancellation at
    low–mid ND); it is not a real competitor — catastrophic at ND 45% / nonlin.

**Margin-uncertainty fix + re-run.  [DONE 2026-08-03]**
`fit_shash_margin` now returns the parameter covariance `V` (inverse observed
information, `optimHess`); `draw_margin()` samples `(μ,logσ,ε)* ~ N(θ̂,V)` per
imputation so Rubin pooling carries margin/tail uncertainty (scales up with
censoring automatically); plug-in fallback when `V` isn't PD. Bake-off re-run at
M=20 with the **MI point estimate (qbar)** as the bias metric (the correct metric
once margins are drawn). Result — **copula is now best on BOTH axes:**

  | method | mean \|bias\| | bias-pass | coverage-pass | both |
  |---|---|---|---|---|
  | ridge  | 0.058 | 77% | 73% | 73% |
  | tobit  | 0.110 | 54% | 85% | 54% |
  | **copula** | **0.029** | 96% | **96%** | **92%** |

  - **Coverage recovered:** copula ND=45% coverage 0.835 → **0.941**; nonlin
    0.75 → **0.89–0.95**. copula now meets both targets in 92% of cells (best).
  - **Bias tradeoff (honest):** the margin draw adds a small upward Jensen bias
    from the nonlinear back-transform, worst at the extreme corner **sas (matched)
    × strong skew 0.75 × ND 45%**, where \|bias\| ≈ 0.12–0.135 (just over 0.10).
    Everywhere else copula bias stays ≤0.06. This corner is already median-first
    territory (ND≥45%, high skew).
  - **nonlin** (dependence broken) coverage now ~0.87–0.92 — the small residual is
    a *dependence*-model limit (linear latent conditional), addressable only by the
    E3 non-linear conditional, not the margin.

**Jensen-bias curb.  [DONE 2026-08-03]**
Added a `margin_draw` control (`gsimp_impute`/`copula_impute`): `FALSE` plugs in
the MLE margin (near-unbiased point), `TRUE` draws (variance). Recommended MI =
**plug-in point + drawn-margin variance** — keeps the coverage fix without the
Jensen bias the all-drawn point picked up at heavy censoring × strong skew. Wired
into the bake-off; full re-run at M=20:

  | method | mean \|bias\| | bias-pass | coverage-pass | both |
  |---|---|---|---|---|
  | ridge  | 0.058 | 77% | 73% | 73% |
  | tobit  | 0.110 | 54% | 85% | 54% |
  | **copula (final)** | **0.015** | **100%** | 92% | **92%** |

  - Best of both: bias-pass back to **100%** (corner `sas strong × ND45` 0.113 →
    0.062) AND coverage kept (ND45 = 0.935). copula overall |bias| 0.015.
  - **Only remaining sub-nominal coverage: the four `nonlin × ND45` cells
    (0.86–0.90)** — every margin-broken family (sas/gamma/lnorm) is now nominal at
    ND45. This residual is the *dependence*-model limit, addressable ONLY by E3
    (non-linear latent conditional), not the margin.

**Wide-p / heteroscedastic stress.  [DONE 2026-08-03]**
`stress_widep_hetero.R` (parallel, ~20 min; n=150, skew 0.75, ND 35%, tobit vs
copula-hybrid, 20 reps, M=15). **copula passes all 6 scenarios; tobit fails bias
on all (skew).**

  | scenario | copula bias | copula coverage |
  |---|---|---|
  | wide p=5 / 15 / 50 | −0.024 / −0.016 / −0.032 | 1.00 |
  | **wide p=150** (p/n=1, PCA regime) | **−0.027** | 0.999 |
  | hetero, marginal-scale (sd 0.4–2.0) | −0.011 | 1.00 |
  | hetero, conditional (shared volatility) | **+0.085** | 0.988 |

  - **Wide p: robust.** copula holds bias ≈ −0.02..−0.03 and coverage ≈ 1.00 all
    the way to p = n = 150 — the many per-analyte margin fits + latent-space PCA do
    not degrade it (mild drift only, within tolerance).
  - **Marginal-scale heterogeneity: a near-freebie** (bias −0.011) — the per-margin
    standardisation absorbs differing analyte scales, as expected.
  - **Conditional heteroscedasticity** (breaks the homoscedastic-latent
    assumption, like `nonlin`) is copula's hardest case here: bias +0.085 (its
    largest, still < 0.10) with coverage 0.988. It passes at ND 35%; expect it to
    be the soft spot at heavier censoring (same conditional-model limit as `nonlin`
    → E3).

  **Status:** S1 (copula) validated and complete as the skew fix — lowest bias
  across all families incl. anti-inverse-crime `nonlin`, near-nominal calibration,
  and now robust to wide-p and heteroscedasticity. Ship as a supported **opt-in
  `imp_model="copula"` for skewed data** (recommend the plug-in/drawn hybrid for
  MI); median-first still fine as a simple default.

**E3 — non-linear latent conditional.  [ATTEMPTED — NEGATIVE 2026-08-03]**
Implemented `spline_tobit` (natural-spline basis on PCA-reduced predictors +
interval-censored `survreg`) and `imp_model="copula_nl"` (copula margins + that
non-linear latent conditional) in `R/impute-nonlinear.R`. Aim: close the residual
under-coverage on `nonlin`/conditional-hetero at heavy ND.
- **First cut DIVERGED:** a flexible basis produces extreme conditional means that,
  through the unbounded non-detect lower tail, feed back and blow up the Gibbs
  sweep (bias −9..−45). Required clamping the conditional mean to the observed
  range (non-linear within range, no runaway) — now stable (bias 0.02–0.07).
- **Even stabilised, it does NOT help — slightly worse.** copula vs copula_nl
  coverage: nonlin 0.7/ND45 **0.871 → 0.825**, nonlin 0.4/ND45 0.917 → 0.908,
  nonlin 0.7/ND35 0.954 → 0.942; conditional-hetero +0.008 (noise). It fails on the
  very cell it targeted.
- **Why:** the residual is a conditional-*shape* problem (homoscedastic-Gaussian
  latent conditional under non-linear dependence), not a conditional-*mean* one; a
  spline bends the mean but keeps the Gaussian shape, and the added estimation
  variance costs more coverage than the bent mean buys (worse at heavier ND).
- **Verdict:** the simple spline conditional is NOT the fix. The `nonlin`×heavy-ND
  residual (~0.87 at the worst corner; ≥0.95 by ND 35%) stands as a documented
  limitation. A real fix would need a non-Gaussian / heteroscedastic conditional
  (or joint model) — a larger undertaking, future work. **`copula` (linear latent)
  remains the recommendation;** `copula_nl` kept only as an experimental/validation
  artifact (do not advertise), pending a decision to keep or revert.
- **G6 — Stress.** heteroscedastic and mixed-skew margins; wide-p (PCA regime);
  small n. Where fixes may break.

---

## 5. Decision criteria & recommended sequence

A method **solves** the skewed-mean problem if, across skew ≤ 0.75 and ND < 50%:
`|mean bias| ≤ 0.10` **and** `coverage ∈ [0.93, 0.97]` (not merely ≥ 0.90 — we
also want to fix the anti-conservative large-n case), **without** materially
worsening RMSE, the symmetric-case results, or runtime — **and** it must hold
under **G5** (family mismatch), not just when DGP = model.

**Sequence:** G1 (confirm) → G2 (map, cheap win for docs) → prototype **S1** →
G3+G4 on S1 vs S0 → if S1 clears G5, ship it as an opt-in `imp_model="copula"`
with S4(a) diagnostic as the safety net; keep S2 as the fallback if the copula's
margin-tail estimation proves fragile, and run S3 once as the engine sanity check.
Escalate to S2/S3 only if S1 fails G5.

**Expected outcome (hypothesis):** keep Gibbs; adopt S1 (copula) as the
skew-robust model; correlation becomes genuinely helpful under skew; the mean
target is recovered to ND < 50% for moderate skew, with S4 diagnostics guarding
the residual high-skew / high-ND corner.

## 6. Open theory to nail down before coding
- Exact MI variance under the copula back-transform (delta method vs full
  posterior-predictive draws through `F_j^{-1}`).
- Bias/variance of the censored marginal-tail estimator, and the copula's
  identification gain as a function of rho (quantify §2's variance-shrinkage).
- Manski bound width for the mean vs ND — is a reported interval useful or too wide?
- Compatibility/convergence of FCS in latent space when margins are *estimated*
  (plug-in `F̂_j` feedback across sweeps).
