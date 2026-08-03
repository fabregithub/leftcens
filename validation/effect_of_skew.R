# =============================================================================
# G2 - the skew sweep  (Phase 3; see validation/PHASE3_SKEW_PLAN.md)
# =============================================================================
#
# D1 found that the tobit imputation's log-scale MEAN is biased downward under
# right-skew, breaching the |bias| <= 0.10 target by ~35% ND, while the MEDIAN
# stays exact. D1 only had two skew points (0 and 0.75). This driver maps the
# boundary finely: mean/median bias, RMSE and MI coverage over a
# skew x non-detect-fraction grid, and interpolates the CROSSING SKEW at which
# |mean bias| first exceeds the tolerance, for each ND level. That turns the
# vignette's qualitative "skew hurts" caveat into a concrete threshold.
#
# Mechanism under test (PHASE3_SKEW_PLAN.md 1.2): bias ~ pi * Delta_tail, so it
# should grow with BOTH skew (deepens Delta_tail) and ND (raises pi). The median,
# being at/above the limit for ND < 50%, should stay ~0 throughout.
#
# Parallel over (skew, nd, rep) tasks (bit-identical to serial; BLAS pinned to
# one thread per worker for the OpenBLAS-pthreads build).
#
# Usage:
#   ./                     # sensible defaults (see below)
#   REPS=200 M=20 WORKERS=22 Rscript validation/effect_of_skew.R
#   SKEWS=0,0.25,0.5,0.75,1 NDS=0.15,0.25,0.35,0.45 Rscript validation/effect_of_skew.R
# =============================================================================

suppressPackageStartupMessages(library(parallel))
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
source("validation/mc_validation.R", local = TRUE)

n        <- as.integer(Sys.getenv("N", "150"))
p        <- as.integer(Sys.getenv("P", "6"))
rho      <- as.numeric(Sys.getenv("RHO", "0.5"))
skews    <- as.numeric(strsplit(Sys.getenv("SKEWS", "0,0.25,0.5,0.75,1"), ",")[[1]])
nds      <- as.numeric(strsplit(Sys.getenv("NDS", "0.15,0.25,0.35,0.45"), ",")[[1]])
reps     <- as.integer(Sys.getenv("REPS", "150"))
iters    <- as.integer(Sys.getenv("ITERS", "20"))
M        <- as.integer(Sys.getenv("M", "20"))
model    <- Sys.getenv("IMP_MODEL", "tobit")
BIAS_TOL <- as.numeric(Sys.getenv("BIAS_TOL", "0.10"))
workers  <- as.integer(Sys.getenv("WORKERS",
                                  as.character(max(1L, detectCores() - 2L))))
base_seed <- 71000L

message(sprintf("G2 skew sweep: n=%d p=%d rho=%.1f model=%s | reps=%d M=%d iters=%d workers=%d",
                n, p, rho, model, reps, M, iters, workers))
message(sprintf("  skews = {%s}  x  ND = {%s}",
                paste(skews, collapse = ","), paste(100 * nds, collapse = ",")))

# one (skew, nd, rep) task -> per-rep metrics
run_one <- function(task) {
  sk <- task$skew; nd <- task$nd; r <- task$rep
  seed <- base_seed + match(sk, skews) * 100000L + match(nd, nds) * 1000L + r
  set.seed(seed)
  z <- simulate_truth(n, p, rho, skew = sk)
  b <- build_bounds(censor_three_tier(z, nd, 0))          # pure left-censoring
  set.seed(seed + 1L)
  f <- gsimp_impute(b, iters_all = iters, imp_model = model)
  mean_bias <- mean(colMeans(f) - colMeans(z))
  med_bias  <- mean(vapply(seq_len(p), function(j)
    stats::median(f[, j]) - stats::median(z[, j]), numeric(1)))
  cov <- mi_coverage(b, z, M = M, iters_all = iters, imp_model = model,
                     seed = seed + 1000L)
  data.frame(skew = sk, nd = nd, rep = r, mean_bias = mean_bias,
             med_bias = med_bias, coverage = cov$coverage, ci_width = cov$ci_width)
}

tasks <- list()
for (sk in skews) for (nd in nds) for (r in seq_len(reps))
  tasks[[length(tasks) + 1]] <- list(skew = sk, nd = nd, rep = r)

t0 <- Sys.time()
res_list <- mclapply(tasks, run_one, mc.cores = workers, mc.preschedule = TRUE)
bad <- vapply(res_list, function(x) !is.data.frame(x), logical(1))
if (any(bad)) message(sprintf("  WARNING: %d task(s) failed", sum(bad)))
per_rep <- do.call(rbind, res_list[!bad])
message(sprintf("  %d tasks in %.1f min", length(tasks),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# aggregate per (skew, nd)
agg <- do.call(rbind, lapply(
  split(per_rep, interaction(per_rep$skew, per_rep$nd, drop = TRUE)),
  function(d) {
    mcse <- function(x) stats::sd(x) / sqrt(length(x))
    data.frame(
      skew = d$skew[1], nd = d$nd[1], n_rep = nrow(d),
      mean_bias = mean(d$mean_bias), mean_bias_mcse = mcse(d$mean_bias),
      abs_mean_bias = abs(mean(d$mean_bias)),
      rmse = sqrt(mean(d$mean_bias^2)),
      med_absbias = mean(abs(d$med_bias)),
      coverage = mean(d$coverage), coverage_mcse = mcse(d$coverage),
      meets_target = abs(mean(d$mean_bias)) <= BIAS_TOL & mean(d$coverage) >= 0.90,
      row.names = NULL)
  }))
agg <- agg[order(agg$nd, agg$skew), ]

# crossing skew per ND: linear interpolation of |mean bias| to the tolerance
crossing <- do.call(rbind, lapply(split(agg, agg$nd), function(d) {
  d <- d[order(d$skew), ]
  y <- d$abs_mean_bias - BIAS_TOL
  cs <- NA_real_
  if (all(y <= 0)) {
    cs <- Inf                                   # never crosses in the swept range
  } else if (y[1] > 0) {
    cs <- 0                                     # already over at skew 0
  } else {
    k <- which(y[-1] > 0 & y[-length(y)] <= 0)[1]   # first up-crossing
    if (!is.na(k)) cs <- d$skew[k] + (d$skew[k + 1] - d$skew[k]) *
        (-y[k]) / (y[k + 1] - y[k])
  }
  data.frame(nd = d$nd[1], crossing_skew = cs)
}))

dir.create("validation/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(config = list(n = n, p = p, rho = rho, model = model, reps = reps,
                           M = M, iters = iters, bias_tol = BIAS_TOL),
             per_rep = per_rep, summary = agg, crossing = crossing),
        "validation/results/effect_of_skew.rds")
utils::write.csv(agg, "validation/results/effect_of_skew.csv", row.names = FALSE)

message("\n== G2 skew sweep (model = ", model, ", n = ", n, ", tol = ", BIAS_TOL, ") ==")
print(agg[, c("skew", "nd", "mean_bias", "abs_mean_bias", "rmse",
              "med_absbias", "coverage", "meets_target")],
      row.names = FALSE, digits = 3)

message("\n== Crossing skew: |mean bias| first exceeds ", BIAS_TOL, " ==")
cr <- crossing
cr$note <- ifelse(is.infinite(cr$crossing_skew), "safe over full sweep",
           ifelse(cr$crossing_skew == 0, "biased even at skew 0",
                  sprintf("crosses near skew %.2f", cr$crossing_skew)))
print(cr, row.names = FALSE, digits = 3)
message("\nRead: at each ND, right-skew up to ~crossing_skew keeps the mean within tolerance.")
message("Median abs bias should stay ~0 throughout (the robust estimand for ND < 50%).")
