# =============================================================================
# E2 - effect of the number of variables (p) imputed together
# =============================================================================
#
# Does the ND < 50% rescue survive as the analyte count grows, including the
# p >= n regime where tobit switches to a PCA-reduced fit? Fixes a
# near-target censoring level and sweeps p across and beyond n, for ridge vs
# tobit. Point-recovery metrics only (mean bias, RMSE of the mean estimate,
# median bias) -- fast, so the full p sweep is feasible; MI coverage at large p
# is a targeted follow-up.
#
# Measured against the north-star target (ND < 50%): does |mean bias| stay
# <= 0.10 as p grows at ~40% ND?
#
# Usage (env-configurable):
#   REPS=30 ITERS=20 Rscript validation/effect_of_p.R
#   N=300 ND=0.45 PS=6,25,100,300 Rscript validation/effect_of_p.R
# =============================================================================

suppressPackageStartupMessages(library(leftcens))
source("validation/mc_validation.R", local = TRUE)

n     <- as.integer(Sys.getenv("N", "100"))
rho   <- as.numeric(Sys.getenv("RHO", "0.5"))
nd    <- as.numeric(Sys.getenv("ND", "0.40"))     # near the <50% target edge
reps  <- as.integer(Sys.getenv("REPS", "12"))
iters <- as.integer(Sys.getenv("ITERS", "12"))
ps    <- as.integer(strsplit(Sys.getenv("PS", "3,6,12,25,50,100,150,200"),
                             ",")[[1]])
models <- trimws(strsplit(Sys.getenv("MODELS", "ridge,tobit"), ",")[[1]])
BIAS_TOL <- as.numeric(Sys.getenv("BIAS_TOL", "0.10"))

dir.create("validation/results", showWarnings = FALSE, recursive = TRUE)
accum_file <- "validation/results/effect_of_p.rds"
csv_file <- "validation/results/effect_of_p.csv"

message(sprintf("E2 effect of p: n=%d rho=%.1f ND=%.0f%% | p in {%s} | models=%s",
                n, rho, 100 * nd, paste(ps, collapse = ","),
                paste(models, collapse = ",")))
message(sprintf("reps=%d iters=%d | target: |mean bias| <= %.2f at ND=%.0f%%",
                reps, iters, BIAS_TOL, 100 * nd))

rows <- list()
for (model in models) {
  for (pp in ps) {
    mb <- medb <- numeric(reps)
    t0 <- Sys.time()
    for (r in seq_len(reps)) {
      set.seed(6000L + pp * 13L + r)   # same data per (p, r) across models
      z <- simulate_truth(n, pp, rho)
      b <- build_bounds(censor_three_tier(z, nd, 0))
      f <- gsimp_impute(b, iters_all = iters, imp_model = model)
      mb[r] <- mean(colMeans(f) - colMeans(z))
      medb[r] <- mean(vapply(seq_len(pp),
                             function(j) stats::median(f[, j]) -
                               stats::median(z[, j]), numeric(1)))
    }
    row <- data.frame(
      model = model, p = pp, n = n, rho = rho, nd = nd,
      mean_bias = mean(mb),
      rmse = sqrt(mean(mb^2)),
      median_bias = mean(medb),
      meets_target = abs(mean(mb)) <= BIAS_TOL,
      secs = as.numeric(difftime(Sys.time(), t0, units = "secs"))
    )
    rows[[length(rows) + 1]] <- row
    # Incremental save so partial progress survives a long background run.
    saveRDS(do.call(rbind, rows), accum_file)
    utils::write.csv(do.call(rbind, rows), csv_file, row.names = FALSE)
    message(sprintf("  %-6s p=%-4d mean_bias=%+.3f rmse=%.3f median_bias=%+.3f  (%.0fs)",
                    model, pp, row$mean_bias, row$rmse, row$median_bias, row$secs))
  }
}

res <- do.call(rbind, rows)
message("\n== E2: recovery vs number of variables (ND=", 100 * nd, "%) ==")
print(res[, c("model", "p", "mean_bias", "rmse", "median_bias", "meets_target")],
      row.names = FALSE, digits = 3)
message("\nRescue holds where meets_target = TRUE (|mean bias| <= ", BIAS_TOL, ").")
