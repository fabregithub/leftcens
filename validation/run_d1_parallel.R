# =============================================================================
# D1 -- definitive high-replication validation run, PARALLEL driver
# =============================================================================
#
# Runs the same study as mc_validation.R's main block, but distributes the
# (scenario x replication) tasks across many cores. Results are bit-identical
# to the serial run: every task sets its own deterministic seed
# (base_seed + scenario*1e5 + rep), so worker assignment never affects output.
#
# The workload is embarrassingly parallel -- each replication simulates,
# censors, imputes and scores an independent dataset; nothing is shared. We
# flatten all tasks, shuffle them (so heavy and cheap conditions are spread
# evenly), split into chunks, and run each chunk with mclapply. Chunking buys
# us progress/ETA logging and periodic checkpoints without hurting throughput.
#
# IMPORTANT (macOS + OpenBLAS pthreads build): each worker must run BLAS
# single-threaded, else 22 forks x N BLAS threads oversubscribes the machine
# (and forking a process with live OpenBLAS threads can deadlock). The launcher
# and this script both pin *_NUM_THREADS=1.
#
# Configure via env (all optional):
#   CONFIG=full N_REP=500 M=50 IMP_MODEL=tobit ITERS_ALL=30 \
#   WORKERS=22 N_CHUNKS=100 BASE_SEED=20240101 \
#   Rscript validation/run_d1_parallel.R
#
# Outputs (validation/results/):
#   d1_<model>_<stamp>.rds   full run (config, grid, per-rep results, summary)
#   d1_latest.rds            copy of the most recent run
#   d1_<model>_<stamp>.csv   per-replication results
#   d1_<model>_<stamp>_summary.csv
#   d1_<model>_<stamp>.checkpoint.rds   partial results, rewritten each chunk
# =============================================================================

suppressPackageStartupMessages(library(parallel))

# Pin BLAS to one thread per worker BEFORE any heavy linear algebra. Set here
# too (not only in the launcher) so a bare `Rscript run_d1_parallel.R` is safe.
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")

# Reuse the engine (simulate_truth, censor_three_tier, build_bounds via
# leftcens, naive_fill, metrics_point, mi_coverage, make_grid). Sourcing does
# NOT trigger mc_validation.R's main block -- it is guarded by sys.nframe()==0.
source("validation/mc_validation.R")

# ---- config ----------------------------------------------------------------
config    <- Sys.getenv("CONFIG", "full")
n_rep     <- as.integer(Sys.getenv("N_REP", "500"))
M         <- as.integer(Sys.getenv("M", "50"))
iters_all <- as.integer(Sys.getenv("ITERS_ALL", "30"))
imp_model <- Sys.getenv("IMP_MODEL", "tobit")
base_seed <- as.integer(Sys.getenv("BASE_SEED", "20240101"))
workers   <- as.integer(Sys.getenv("WORKERS",
                                   as.character(max(1L, parallel::detectCores() - 2L))))
n_chunks  <- as.integer(Sys.getenv("N_CHUNKS", "100"))

grid   <- make_grid(config)
n_scen <- nrow(grid)
n_task <- n_scen * n_rep

logmsg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                                        sprintf(...)))

logmsg("D1 parallel run starting")
logmsg("config=%s  scenarios=%d  n_rep=%d  M=%d  iters_all=%d  model=%s",
       config, n_scen, n_rep, M, iters_all, imp_model)
logmsg("tasks=%d  workers=%d  chunks=%d  base_seed=%d",
       n_task, workers, n_chunks, base_seed)

# ---- one task (identical body to run_validation's inner loop) ---------------
run_one <- function(task) {
  s <- task[[1]]; r <- task[[2]]
  cond <- grid[s, ]
  seed <- base_seed + s * 100000L + r
  set.seed(seed)
  z <- simulate_truth(cond$n, cond$p, cond$rho, skew = cond$skew)
  bnds <- build_bounds(censor_three_tier(z, cond$nd_frac, cond$dnq_frac))

  set.seed(seed + 1L)
  f_gs <- gsimp_impute(bnds, iters_all = iters_all, imp_model = imp_model)
  f_nv <- naive_fill(bnds, "limit")

  m_gs <- metrics_point(f_gs, z)
  m_nv <- metrics_point(f_nv, z)
  cov <- mi_coverage(bnds, z, M = M, iters_all = iters_all,
                     imp_model = imp_model, seed = seed + 1000L)

  data.frame(
    scenario = s, rep = r,
    n = cond$n, p = cond$p, rho = cond$rho,
    nd_frac = cond$nd_frac, dnq_frac = cond$dnq_frac, skew = cond$skew,
    censored_frac = mean(bnds$to_impute),
    gs_mean_bias = m_gs$mean_bias,
    gs_mean_absbias = m_gs$mean_absbias,
    gs_q50_absbias = m_gs$q50_absbias,
    gs_cor_mae = m_gs$cor_mae,
    naive_mean_bias = m_nv$mean_bias,
    naive_mean_absbias = m_nv$mean_absbias,
    naive_cor_mae = m_nv$cor_mae,
    mi_coverage = cov$coverage,
    mi_ci_width = cov$ci_width
  )
}

# ---- build + shuffle the flat task list ------------------------------------
tasks <- expand.grid(rep = seq_len(n_rep), scenario = seq_len(n_scen))
set.seed(20240101L)                         # reproducible shuffle (order only)
tasks <- tasks[sample.int(nrow(tasks)), , drop = FALSE]
task_list <- Map(function(s, r) c(s, r), tasks$scenario, tasks$rep)

chunk_id  <- cut(seq_len(n_task), breaks = min(n_chunks, n_task), labels = FALSE)
chunks    <- split(task_list, chunk_id)

dir.create("validation/results", showWarnings = FALSE, recursive = TRUE)
stamp   <- format(Sys.time(), "%Y%m%d-%H%M%S")
ckpt    <- sprintf("validation/results/d1_%s_%s.checkpoint.rds", imp_model, stamp)

# ---- run chunks ------------------------------------------------------------
t0 <- Sys.time()
acc  <- vector("list", length(chunks))
done <- 0L
for (ci in seq_along(chunks)) {
  res_ci <- mclapply(chunks[[ci]], run_one,
                     mc.cores = workers, mc.preschedule = TRUE)

  errs <- vapply(res_ci, function(x) inherits(x, "try-error") ||
                   inherits(x, "simpleError"), logical(1))
  if (any(errs)) {
    logmsg("WARNING: %d task(s) failed in chunk %d; first error: %s",
           sum(errs), ci, conditionMessage(attr(res_ci[[which(errs)[1]]],
                                                 "condition")))
    res_ci <- res_ci[!errs]
  }
  acc[[ci]] <- do.call(rbind, res_ci)

  done <- done + length(chunks[[ci]])
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  rate <- done / elapsed
  eta_h <- (n_task - done) / rate / 3600
  logmsg("chunk %d/%d | %d/%d tasks (%.1f%%) | %.1f tasks/s | elapsed %.2f h | ETA %.2f h",
         ci, length(chunks), done, n_task, 100 * done / n_task,
         rate, elapsed / 3600, eta_h)

  saveRDS(do.call(rbind, acc[seq_len(ci)]), ckpt)   # crash-safe partial results
}

res  <- do.call(rbind, acc)
res  <- res[order(res$scenario, res$rep), ]         # tidy, deterministic order
summ <- summarise_validation(res)

# ---- persist (same shape as mc_validation.R) -------------------------------
payload <- list(config = config, n_rep = n_rep, M = M, imp_model = imp_model,
                iters_all = iters_all, workers = workers,
                grid = grid, results = res, summary = summ,
                elapsed_hours = as.numeric(difftime(Sys.time(), t0, units = "hours")))

saveRDS(payload, sprintf("validation/results/d1_%s_%s.rds", imp_model, stamp))
saveRDS(payload, sprintf("validation/results/d1_%s_latest.rds", imp_model))
utils::write.csv(res,  sprintf("validation/results/d1_%s_%s.csv", imp_model, stamp),
                 row.names = FALSE)
utils::write.csv(summ, sprintf("validation/results/d1_%s_%s_summary.csv", imp_model, stamp),
                 row.names = FALSE)
if (file.exists(ckpt)) file.remove(ckpt)            # completed -> drop checkpoint

logmsg("DONE in %.2f h. Wrote d1_%s_%s.rds (+ latest, csv, summary).",
       payload$elapsed_hours, imp_model, stamp)

message("\n== D1 summary (bias, RMSE, MI coverage by scenario) ==")
print(summ[, c("n", "p", "rho", "nd_frac", "skew", "censored_frac",
               "gs_bias", "gs_rmse", "naive_bias", "naive_rmse",
               "mi_coverage")], digits = 3)
