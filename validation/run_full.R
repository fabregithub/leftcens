# =============================================================================
# Full-grid validation runner with cross-call accumulation.
# =============================================================================
#
# The full mc_validation grid at publication replication takes many hours, which
# exceeds a single foreground slot. This driver runs the full grid in batches:
# each MODE=batch call adds `N_REP` replications (with a batch-specific seed) and
# appends to validation/results/full_accum.rds. MODE=summarise pools all batches
# and writes results/latest.rds for report.Rmd.
#
#   MODE=batch BATCH=1 N_REP=10 M=10 Rscript validation/run_full.R
#   MODE=batch BATCH=2 N_REP=10 M=10 Rscript validation/run_full.R
#   ...
#   MODE=summarise Rscript validation/run_full.R
#
# For the definitive run on capable hardware, prefer a single high-replication
# pass:  CONFIG=full N_REP=500 M=50 Rscript validation/mc_validation.R
# =============================================================================

suppressPackageStartupMessages(library(leftcens))
source("validation/mc_validation.R", local = TRUE)

mode <- Sys.getenv("MODE", "batch")
config <- Sys.getenv("CONFIG", "full")
grid <- make_grid(config)
dir.create("validation/results", showWarnings = FALSE, recursive = TRUE)
accum_file <- "validation/results/full_accum.rds"

if (mode == "batch") {
  batch <- as.integer(Sys.getenv("BATCH", "1"))
  n_rep <- as.integer(Sys.getenv("N_REP", "10"))
  M <- as.integer(Sys.getenv("M", "10"))
  iters_all <- as.integer(Sys.getenv("ITERS_ALL", "30"))
  base_seed <- 20240101L + batch * 7919L        # distinct seeds per batch

  t0 <- Sys.time()
  message(sprintf("Batch %d: %s grid (%d scenarios) x %d reps, M=%d, iters=%d",
                  batch, config, nrow(grid), n_rep, M, iters_all))
  res <- run_validation(grid, n_rep = n_rep, M = M, iters_all = iters_all,
                        imp_model = "ridge", base_seed = base_seed,
                        verbose = FALSE)
  res$batch <- batch

  if (file.exists(accum_file)) res <- rbind(readRDS(accum_file), res)
  saveRDS(res, accum_file)
  message(sprintf("Batch %d done in %.0f s. Cumulative reps/scenario: %.0f",
                  batch, as.numeric(difftime(Sys.time(), t0, units = "secs")),
                  nrow(res) / nrow(grid)))
} else if (mode == "summarise") {
  res <- readRDS(accum_file)
  summ <- summarise_validation(res)
  reps_per <- nrow(res) / nrow(grid)
  saveRDS(list(config = config, n_rep = reps_per, M = NA_integer_,
               imp_model = "ridge", grid = grid, results = res,
               summary = summ),
          "validation/results/latest.rds")
  utils::write.csv(summ, "validation/results/full_summary.csv",
                   row.names = FALSE)
  message(sprintf("Pooled %d batches -> %.0f reps/scenario over %d scenarios",
                  length(unique(res$batch)), reps_per, nrow(grid)))
  message("\n== Full-grid summary (bias, RMSE, MI coverage) ==")
  print(summ[order(summ$censored_frac),
             c("n", "p", "rho", "nd_frac", "skew", "censored_frac",
               "gs_bias", "gs_rmse", "naive_bias", "mi_coverage",
               "mi_coverage_mcse")],
        row.names = FALSE, digits = 3)
} else {
  stop("MODE must be 'batch' or 'summarise'.")
}
