# =============================================================================
# Heterogeneous detection-limit experiment for gsimp_impute()
# =============================================================================
#
# The detection-rate sweep (detection_rate_standard.R) found that under a
# UNIFORM detection limit, inter-analyte correlation did NOT extend the reliable
# range. The mechanistic explanation: with exchangeable correlation and a shared
# limit, censoring is itself correlated across analytes -- when one analyte's
# value is low (and censored), the correlated analytes tend to be low (and
# censored) too, so the conditional model has no observed neighbour to borrow
# from exactly for the rows that need it.
#
# This experiment tests the prediction that follows: with HETEROGENEOUS limits
# -- some analytes well observed ("anchors"), others heavily censored
# ("targets") -- correlation SHOULD rescue the heavily-censored analytes, because
# the anchors stay observed and inform them.
#
# Design: p = n_anchor + n_target analytes. Anchors censored lightly
# (`nd_anchor`), targets censored heavily (`nd_target`, in the regime where the
# uniform case failed). Sweep correlation `rho`. Report bias and MI coverage of
# the TARGET analytes as a function of rho. If they improve with rho, correlation
# helps once limits are heterogeneous.
#
# Usage:
#   Rscript validation/heterogeneous_limits.R
#   N_REP=200 M=25 Rscript validation/heterogeneous_limits.R
# =============================================================================

suppressPackageStartupMessages(library(leftcens))
source("validation/mc_validation.R", local = TRUE)

n         <- as.integer(Sys.getenv("N", "200"))
n_anchor  <- as.integer(Sys.getenv("N_ANCHOR", "4"))
n_target  <- as.integer(Sys.getenv("N_TARGET", "4"))
nd_anchor <- as.numeric(Sys.getenv("ND_ANCHOR", "0.05"))
nd_target <- as.numeric(Sys.getenv("ND_TARGET", "0.40"))
rhos      <- as.numeric(strsplit(Sys.getenv("RHOS", "0,0.3,0.6,0.8"), ",")[[1]])
n_rep     <- as.integer(Sys.getenv("N_REP", "60"))
M         <- as.integer(Sys.getenv("M", "15"))
iters_all <- as.integer(Sys.getenv("ITERS_ALL", "30"))
BIAS_TOL     <- as.numeric(Sys.getenv("BIAS_TOL", "0.10"))
COVERAGE_TOL <- as.numeric(Sys.getenv("COVERAGE_TOL", "0.90"))
models       <- trimws(strsplit(Sys.getenv("MODELS", "ridge,tobit"), ",")[[1]])

p <- n_anchor + n_target
nd_fracs <- c(rep(nd_anchor, n_anchor), rep(nd_target, n_target))
anchor_idx <- seq_len(n_anchor)
target_idx <- n_anchor + seq_len(n_target)

# Per-analyte MI coverage of the true mean (returns a logical vector length p).
mi_covered_per_analyte <- function(bnds, z_true, M, iters_all, seed, imp_model) {
  n <- nrow(z_true); p <- ncol(z_true)
  means <- matrix(NA_real_, M, p)
  vars <- matrix(NA_real_, M, p)
  for (mi in seq_len(M)) {
    set.seed(seed + mi)
    fm <- gsimp_impute(bnds, iters_all = iters_all, imp_model = imp_model)
    means[mi, ] <- colMeans(fm)
    vars[mi, ] <- apply(fm, 2, stats::var) / n
  }
  qbar <- colMeans(means)
  ubar <- colMeans(vars)
  b <- apply(means, 2, stats::var)
  Tvar <- ubar + (1 + 1 / M) * b
  df <- (M - 1) * (1 + ubar / ((1 + 1 / M) * pmax(b, .Machine$double.eps)))^2
  df[!is.finite(df)] <- M - 1
  half <- stats::qt(0.975, df) * sqrt(Tvar)
  true_mean <- colMeans(z_true)
  (true_mean >= qbar - half) & (true_mean <= qbar + half)
}

message(sprintf(
  "Heterogeneous-limit experiment: %d anchor analytes @ %.0f%% ND, %d target @ %.0f%% ND",
  n_anchor, 100 * nd_anchor, n_target, 100 * nd_target))
message(sprintf("Models = %s | sweeping rho = %s | %d reps, M=%d, iters=%d",
                paste(models, collapse = ", "),
                paste(rhos, collapse = ", "), n_rep, M, iters_all))

rows <- list()
for (model in models) {
  for (rho in rhos) {
    tb <- ab <- numeric(n_rep)          # target / anchor mean bias
    tc <- ac <- numeric(n_rep)          # target / anchor coverage
    for (r in seq_len(n_rep)) {
      # Seed depends only on (rho, r), so every model sees identical data (paired)
      seed <- 900000L + as.integer(rho * 1000) * 1000L + r
      set.seed(seed)
      z <- simulate_truth(n, p, rho)
      bnds <- build_bounds(censor_three_tier(z, nd_fracs, 0))

      set.seed(seed + 1L)
      f <- gsimp_impute(bnds, iters_all = iters_all, imp_model = model)
      bias_per <- colMeans(f) - colMeans(z)
      tb[r] <- mean(bias_per[target_idx])
      ab[r] <- mean(bias_per[anchor_idx])

      cov_per <- mi_covered_per_analyte(bnds, z, M, iters_all, seed + 1000L, model)
      tc[r] <- mean(cov_per[target_idx])
      ac[r] <- mean(cov_per[anchor_idx])
    }
    rows[[length(rows) + 1]] <- data.frame(
      model = model, rho = rho,
      target_bias = mean(tb), target_bias_mcse = stats::sd(tb) / sqrt(n_rep),
      target_coverage = mean(tc),
      anchor_bias = mean(ab), anchor_coverage = mean(ac),
      target_reliable = abs(mean(tb)) <= BIAS_TOL & mean(tc) >= COVERAGE_TOL
    )
    message(sprintf("  %-6s rho=%.1f: target bias=%+.3f coverage=%.2f",
                    model, rho, mean(tb), mean(tc)))
  }
}
res <- do.call(rbind, rows)

dir.create("validation/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(config = list(n = n, n_anchor = n_anchor, n_target = n_target,
                           nd_anchor = nd_anchor, nd_target = nd_target,
                           models = models, n_rep = n_rep, M = M,
                           iters_all = iters_all),
             results = res),
        "validation/results/heterogeneous_limits.rds")
utils::write.csv(res, "validation/results/heterogeneous_limits.csv",
                 row.names = FALSE)

message(sprintf(
  "\n== Target-analyte reliability vs correlation, by model (targets @ %.0f%% ND, anchors @ %.0f%% ND) ==",
  100 * nd_target, 100 * nd_anchor))
print(res[order(res$model, res$rho),
          c("model", "rho", "target_bias", "target_coverage",
            "anchor_coverage", "target_reliable")],
      row.names = FALSE, digits = 3)
message("\nInterpretation: does target_coverage rise with rho? For the observed-only")
message("model (ridge) it did not (selection bias). If the censored model (tobit)")
message("lets well-observed anchors rescue heavily-censored targets, its target")
message("bias/coverage should hold up where ridge's does not.")
