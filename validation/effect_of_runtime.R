# =============================================================================
# E7 - runtime / scalability of gsimp_impute (tobit vs ridge)
# =============================================================================
# Times one gsimp_impute() call across models and problem sizes, to feed the
# compute-cost line of the guidance. Timing only -- no Monte Carlo needed.
#
#   Rscript validation/effect_of_runtime.R
# =============================================================================

suppressPackageStartupMessages(library(leftcens))
source("validation/mc_validation.R", local = TRUE)

iters <- as.integer(Sys.getenv("ITERS", "10"))
reps  <- as.integer(Sys.getenv("REPS", "3"))     # timing repeats (averaged)
grid <- data.frame(
  n = c(100L, 300L, 100L, 300L, 100L),
  p = c(6L,   6L,   20L,  50L,  100L)
)

message(sprintf("E7 runtime: iters_all=%d, avg of %d calls, ND=35%%", iters, reps))
rows <- list()
for (i in seq_len(nrow(grid))) {
  n <- grid$n[i]; p <- grid$p[i]
  set.seed(1); z <- simulate_truth(n, p, 0.5)
  b <- build_bounds(censor_three_tier(z, 0.35, 0))
  tm <- sapply(c("ridge", "tobit"), function(m) {
    system.time(for (k in seq_len(reps))
      gsimp_impute(b, iters_all = iters, imp_model = m))[["elapsed"]] / reps
  })
  rows[[i]] <- data.frame(n = n, p = p,
                          ridge_s = tm[["ridge"]], tobit_s = tm[["tobit"]],
                          ratio = tm[["tobit"]] / tm[["ridge"]])
  message(sprintf("  n=%-3d p=%-3d ridge=%.3fs tobit=%.3fs (x%.0f)",
                  n, p, tm[["ridge"]], tm[["tobit"]], tm[["tobit"]] / tm[["ridge"]]))
}
res <- do.call(rbind, rows)
dir.create("validation/results", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(res, "validation/results/effect_of_runtime.csv", row.names = FALSE)
message("\n== E7 runtime (seconds per gsimp_impute call) ==")
print(res, row.names = FALSE, digits = 3)
