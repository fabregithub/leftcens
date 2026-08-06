# =============================================================================
# E8 - how many imputations (m)? stability, FMI, and the effect of sample size
#   (Phase 3 follow-up; see PHASE3_SKEW_PLAN.md and gsimp_mi())
# =============================================================================
#
# E6 studied m for tobit only, on symmetric data, tracking coverage/width. This
# extends it to what matters for guidance and the manuscript:
#   * the COPULA (hybrid MI) vs tobit,
#   * the fraction of missing information (FMI) that governs the required m,
#   * the m = 1 (single-imputation) baseline -- to show it under-covers, and
#   * the effect of SAMPLE SIZE n.
#
# Hypothesis (from field experience: at n > 50k, m ~ 8-10 suffices). Required m
# is set by FMI, not n directly. FMI has an irreducible predictive part (~ the
# censoring fraction, n-invariant) plus an imputation-MODEL-parameter part that
# shrinks as ~1/n. For a PROPER model that draws its parameters (copula, via
# margin_draw), FMI therefore falls with n toward a floor, so required m drops.
# tobit does not draw its conditional parameters, so it should be less
# n-dependent. This driver measures that.
#
# Efficiency: generate M_MAX imputations once per replication, then pool at each
# m in the grid using prefixes (so cost is ~ M_MAX, not sum(m grid)).
#
# Run AFTER the copula capstone (don't oversubscribe): default WORKERS leaves the
# machine some headroom.
#
# Defaults are the manuscript run: REPS=100, NS={75,150,300,600,1200,2400,4800},
# M_MAX=100, both models -- ~4.5 h on 22 workers (copula at the top n dominates).
# Launch it detached with validation/run_m_stability.sh (survives closing the
# terminal). Quick look:
#   REPS=20 NS=75,300,1200 MODELS=copula Rscript validation/effect_of_m_stability.R
# =============================================================================

suppressPackageStartupMessages(library(parallel))
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
source("validation/mc_validation.R", local = TRUE)

p        <- as.integer(Sys.getenv("P", "6"))
rho      <- as.numeric(Sys.getenv("RHO", "0.5"))
nd       <- as.numeric(Sys.getenv("ND", "0.35"))
ns       <- as.integer(strsplit(Sys.getenv("NS", "75,150,300,600,1200,2400,4800"), ",")[[1]])
skews    <- as.numeric(strsplit(Sys.getenv("SKEWS", "0,0.75"), ",")[[1]])
models   <- trimws(strsplit(Sys.getenv("MODELS", "tobit,copula"), ",")[[1]])
m_grid   <- as.integer(strsplit(Sys.getenv("M_GRID", "1,2,3,5,10,20,30,50,75,100"), ",")[[1]])
M_MAX    <- as.integer(Sys.getenv("M_MAX", "100"))
reps     <- as.integer(Sys.getenv("REPS", "100"))
iters    <- as.integer(Sys.getenv("ITERS", "20"))
adapt_tol   <- as.numeric(Sys.getenv("ADAPT_TOL", "0.01"))
adapt_step  <- as.integer(Sys.getenv("ADAPT_STEP", "5"))
workers  <- as.integer(Sys.getenv("WORKERS", as.character(max(1L, detectCores() - 4L))))

# Rubin pooling of per-analyte means from prefix summaries (with FMI).
pool <- function(Pm, Pv, Dm, conf = 0.95) {
  m <- nrow(Dm); qbar <- colMeans(Pm); ubar <- colMeans(Pv)
  if (m < 2L) {                                   # naive single-imputation CI:
    half <- stats::qnorm((1 + conf) / 2) * sqrt(ubar)   # within-var only (ignores
    return(list(qbar = qbar, half = half,               # imputation uncertainty ->
                fmi = NA_real_, width = mean(2 * half))) # under-covers, by design)
  }
  B <- apply(Dm, 2, stats::var); Tv <- ubar + (1 + 1 / m) * B
  r <- (1 + 1 / m) * B / pmax(ubar, .Machine$double.eps)
  df <- (m - 1) * (1 + 1 / r)^2; df[!is.finite(df)] <- m - 1
  lambda <- (r + 2 / (df + 3)) / (r + 1)
  half <- stats::qt((1 + conf) / 2, df) * sqrt(Tv)
  list(qbar = qbar, half = half, fmi = mean(lambda), width = mean(2 * half))
}

run_one <- function(task) {
  n <- task$n; skew <- task$skew; mo <- task$model; r <- task$rep
  set.seed(51000L + task$ti * 1000L + r)
  z <- simulate_truth(n, p, rho, skew = skew)
  b <- build_bounds(censor_three_tier(z, nd, 0))
  truth <- colMeans(z); is_cop <- identical(mo, "copula")

  Pm <- Pv <- matrix(NA_real_, M_MAX, p); Dm <- matrix(NA_real_, M_MAX, p)
  for (mi in seq_len(M_MAX)) {
    set.seed(51000L + task$ti * 1000L + r + 1000L + mi)
    fp <- gsimp_impute(b, iters_all = iters, imp_model = mo, margin_draw = FALSE)
    Pm[mi, ] <- colMeans(fp); Pv[mi, ] <- apply(fp, 2, stats::var) / n
    if (is_cop) {
      set.seed(51000L + task$ti * 1000L + r + 2000L + mi)
      fd <- gsimp_impute(b, iters_all = iters, imp_model = mo, margin_draw = TRUE)
      Dm[mi, ] <- colMeans(fd)
    } else Dm[mi, ] <- Pm[mi, ]
  }

  # pool at each m in the grid (prefixes)
  rows <- lapply(m_grid[m_grid <= M_MAX], function(m) {
    pp <- pool(Pm[seq_len(m), , drop = FALSE], Pv[seq_len(m), , drop = FALSE],
               Dm[seq_len(m), , drop = FALSE])
    cov <- if (all(is.na(pp$half))) NA_real_
           else mean(truth >= pp$qbar - pp$half & truth <= pp$qbar + pp$half)
    data.frame(n = n, skew = skew, model = mo, rep = r, m = m,
               coverage = cov, ci_width = pp$width, fmi = pp$fmi)
  })

  # adaptive m: walk up in steps, stop when mean CI width stabilises
  seqm <- seq(adapt_step, M_MAX, by = adapt_step)
  w_prev <- Inf; chosen <- max(seqm)
  for (m in seqm) {
    w <- pool(Pm[seq_len(m), , drop = FALSE], Pv[seq_len(m), , drop = FALSE],
              Dm[seq_len(m), , drop = FALSE])$width
    if (is.finite(w) && is.finite(w_prev) && abs(w - w_prev) / max(w_prev, 1e-8) < adapt_tol) {
      chosen <- m; break
    }
    w_prev <- w
  }
  attr(rows, "adaptive_m") <- data.frame(n = n, skew = skew, model = mo, rep = r,
                                         adaptive_m = chosen)
  rows
}

tasks <- list(); ti <- 0L
for (n in ns) for (sk in skews) for (mo in models) { ti <- ti + 1L
  for (r in seq_len(reps)) tasks[[length(tasks) + 1L]] <- list(ti = ti, n = n, skew = sk, model = mo, rep = r)
}
message(sprintf("E8 m-stability: p=%d rho=%.1f ND=%.0f%% | n={%s} skew={%s} models={%s}",
                p, rho, 100 * nd, paste(ns, collapse = ","),
                paste(skews, collapse = ","), paste(models, collapse = ",")))
message(sprintf("  m_grid={%s} M_MAX=%d reps=%d | %d tasks, %d workers",
                paste(m_grid, collapse = ","), M_MAX, reps, length(tasks), workers))

t0 <- Sys.time()
out <- mclapply(tasks, run_one, mc.cores = workers, mc.preschedule = TRUE)
res <- do.call(rbind, lapply(out, function(o) do.call(rbind, o)))
adm <- do.call(rbind, lapply(out, function(o) attr(o, "adaptive_m")))
message(sprintf("  done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# aggregate
agg <- do.call(rbind, lapply(
  split(res, interaction(res$n, res$skew, res$model, res$m, drop = TRUE)),
  function(d) data.frame(n = d$n[1], skew = d$skew[1], model = d$model[1], m = d$m[1],
    coverage = mean(d$coverage, na.rm = TRUE), ci_width = mean(d$ci_width, na.rm = TRUE),
    fmi = mean(d$fmi, na.rm = TRUE), row.names = NULL)))
agg <- agg[order(agg$model, agg$skew, agg$n, agg$m), ]
adm_agg <- do.call(rbind, lapply(
  split(adm, interaction(adm$n, adm$skew, adm$model, drop = TRUE)),
  function(d) data.frame(n = d$n[1], skew = d$skew[1], model = d$model[1],
    adaptive_m_median = stats::median(d$adaptive_m),
    adaptive_m_mean = mean(d$adaptive_m), row.names = NULL)))

dir.create("validation/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(per_rep = res, summary = agg, adaptive = adm_agg,
             config = list(p = p, rho = rho, nd = nd, ns = ns, skews = skews,
                           models = models, m_grid = m_grid, reps = reps)),
        "validation/results/effect_of_m_stability.rds")
utils::write.csv(agg, "validation/results/m_stability_summary.csv", row.names = FALSE)
utils::write.csv(adm_agg, "validation/results/m_stability_adaptive.csv", row.names = FALSE)

message("\n== FMI at m=20, by model x skew x n (does FMI fall with n?) ==")
sub <- agg[agg$m == 20, c("model", "skew", "n", "fmi", "coverage")]
print(sub[order(sub$model, sub$skew, sub$n), ], row.names = FALSE, digits = 3)
message("\n== Adaptive m chosen (CI-width tol=", adapt_tol, "), by model x skew x n ==")
print(adm_agg[order(adm_agg$model, adm_agg$skew, adm_agg$n), ], row.names = FALSE, digits = 3)
message("\n== Single-imputation (m=1) vs m=50 coverage -- why MI is needed ==")
w <- reshape(agg[agg$m %in% c(1,50), c("model","skew","n","m","coverage")],
             idvar = c("model","skew","n"), timevar = "m", direction = "wide")
print(w[order(w$model, w$skew, w$n), ], row.names = FALSE, digits = 3)
