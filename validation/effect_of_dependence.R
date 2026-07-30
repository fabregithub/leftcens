# =============================================================================
# E3 - parametric (linear) vs non-parametric (non-linear) dependence
# =============================================================================
#
# tobit's conditional model is LINEAR. Does non-linear inter-analyte dependence
# -- strong dependence with near-zero linear correlation -- break the ND < 50%
# rescue, or does tobit's per-analyte censored likelihood keep bias in check
# (leaving only efficiency on the table)? Compares two data-generating regimes at
# matched marginals and a near-target censoring level:
#   * linear     : simulate_truth(rho)        -- Gaussian/linear correlation
#   * nonlinear  : simulate_truth_nl(strength) -- shared latent via orthogonal
#                  non-linear bases (low linear corr, strong dependence)
#
# Point-recovery metrics (mean bias, RMSE, median bias) for ridge vs tobit.
# Expectation from E1/E2: tobit's *bias* rescue is per-analyte, so it should
# survive non-linear dependence; the interesting effect is on RMSE / efficiency.
#
# Usage:
#   REPS=40 ITERS=20 Rscript validation/effect_of_dependence.R
# =============================================================================

suppressPackageStartupMessages(library(leftcens))
source("validation/mc_validation.R", local = TRUE)

n        <- as.integer(Sys.getenv("N", "150"))
p        <- as.integer(Sys.getenv("P", "6"))
nd       <- as.numeric(Sys.getenv("ND", "0.40"))
rho      <- as.numeric(Sys.getenv("RHO", "0.6"))        # linear regime
strength <- as.numeric(Sys.getenv("STRENGTH", "0.6"))   # non-linear regime
reps     <- as.integer(Sys.getenv("REPS", "30"))
iters    <- as.integer(Sys.getenv("ITERS", "20"))
models   <- trimws(strsplit(Sys.getenv("MODELS", "ridge,tobit"), ",")[[1]])
BIAS_TOL <- as.numeric(Sys.getenv("BIAS_TOL", "0.10"))

gen <- list(
  linear    = function(seed) { set.seed(seed); simulate_truth(n, p, rho) },
  nonlinear = function(seed) { set.seed(seed); simulate_truth_nl(n, p, strength) }
)

message(sprintf("E3 dependence: n=%d p=%d ND=%.0f%% | linear rho=%.1f vs nonlinear strength=%.1f",
                n, p, 100 * nd, rho, strength))
message(sprintf("models=%s reps=%d iters=%d | target |mean bias| <= %.2f",
                paste(models, collapse = ","), reps, iters, BIAS_TOL))

rows <- list()
for (dep in names(gen)) {
  for (model in models) {
    mb <- medb <- lincorr <- numeric(reps)
    for (r in seq_len(reps)) {
      z <- gen[[dep]](5000L + r)          # same data per (dep, r) across models
      off <- upper.tri(diag(p))
      lincorr[r] <- mean(abs(cor(z)[off]))
      b <- build_bounds(censor_three_tier(z, nd, 0))
      f <- gsimp_impute(b, iters_all = iters, imp_model = model)
      mb[r] <- mean(colMeans(f) - colMeans(z))
      medb[r] <- mean(vapply(seq_len(p),
                             function(j) stats::median(f[, j]) -
                               stats::median(z[, j]), numeric(1)))
    }
    rows[[length(rows) + 1]] <- data.frame(
      dependence = dep, model = model,
      mean_lin_corr = mean(lincorr),
      mean_bias = mean(mb), rmse = sqrt(mean(mb^2)),
      median_bias = mean(medb),
      meets_target = abs(mean(mb)) <= BIAS_TOL
    )
    message(sprintf("  %-9s %-6s lin_corr=%.2f mean_bias=%+.3f rmse=%.3f",
                    dep, model, mean(lincorr), mean(mb), sqrt(mean(mb^2))))
  }
}
res <- do.call(rbind, rows)

dir.create("validation/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(res, "validation/results/effect_of_dependence.rds")
utils::write.csv(res, "validation/results/effect_of_dependence.csv",
                 row.names = FALSE)

message("\n== E3: linear vs non-linear dependence (ND=", 100 * nd, "%) ==")
print(res[, c("dependence", "model", "mean_lin_corr", "mean_bias", "rmse",
              "median_bias", "meets_target")], row.names = FALSE, digits = 3)
message("\nIf tobit meets_target under BOTH regimes, the bias-rescue survives")
message("non-linear dependence (per-analyte likelihood). Compare rmse across")
message("regimes to see the efficiency cost of unmodelled non-linear dependence.")
