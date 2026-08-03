# =============================================================================
# Copula capstone -- definitive high-replication validation of imp_model="copula"
# =============================================================================
#
# The tobit/ridge capstone (D1, run_d1_parallel.R) established the skew problem.
# This is the matching high-replication confirmation for the FIX: the same full
# 72-scenario grid, but imputed with the Gaussian-copula model using the
# RECOMMENDED multiple-imputation workflow --
#   * point estimate (qbar): M plug-in-margin imputations (margin_draw = FALSE),
#   * between-imputation variance: M drawn-margin imputations (margin_draw = TRUE),
# so Rubin's-rules coverage is calibrated without the Jensen point bias the
# all-drawn version carries at heavy censoring x strong skew.
#
# Output shape matches D1 (same columns, summarise_validation) so the copula
# capstone is directly comparable to the tobit/ridge one.
#
# Defaults: N_REP=300, M=30 (~12 h on 22 cores). Copula is ~5-10x tobit and the
# hybrid runs 2*M imputations per rep, so exact D1 parity (N_REP=500 M=50, ~30 h)
# is opt-in via env. Bit-identical to serial; BLAS pinned to one thread/worker.
#
#   CONFIG=full N_REP=300 M=30 WORKERS=22 Rscript validation/run_copula_capstone.R
# =============================================================================

suppressPackageStartupMessages(library(parallel))
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
source("validation/mc_validation.R", local = TRUE)   # engine (main is guarded)

config    <- Sys.getenv("CONFIG", "full")
n_rep     <- as.integer(Sys.getenv("N_REP", "300"))
M         <- as.integer(Sys.getenv("M", "30"))
iters_all <- as.integer(Sys.getenv("ITERS_ALL", "30"))
base_seed <- as.integer(Sys.getenv("BASE_SEED", "20240101"))
workers   <- as.integer(Sys.getenv("WORKERS",
                                   as.character(max(1L, parallel::detectCores() - 2L))))
n_chunks  <- as.integer(Sys.getenv("N_CHUNKS", "100"))

grid   <- make_grid(config)
n_scen <- nrow(grid)
n_task <- n_scen * n_rep

logmsg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                                        sprintf(...)))
logmsg("Copula capstone starting (hybrid MI: plug-in point + drawn variance)")
logmsg("config=%s  scenarios=%d  n_rep=%d  M=%d  iters=%d", config, n_scen, n_rep, M, iters_all)
logmsg("tasks=%d  workers=%d  chunks=%d  (2*M=%d copula imputations per rep)",
       n_task, workers, n_chunks, 2L * M)

# ---- one task: hybrid multiple imputation for the copula ---------------------
run_one <- function(task) {
  s <- task[[1]]; r <- task[[2]]
  cond <- grid[s, ]
  seed <- base_seed + s * 100000L + r
  set.seed(seed)
  z <- simulate_truth(cond$n, cond$p, cond$rho, skew = cond$skew)
  bnds <- build_bounds(censor_three_tier(z, cond$nd_frac, cond$dnq_frac))

  n_ <- nrow(z); p_ <- ncol(z)
  truth <- colMeans(z); truth_med <- apply(z, 2, stats::median)
  ctrue <- if (p_ > 1) stats::cor(z) else NULL
  off   <- if (p_ > 1) upper.tri(ctrue) else NULL

  Pmeans <- Pvars <- Pmeds <- matrix(NA_real_, M, p_)
  Dmeans <- matrix(NA_real_, M, p_)
  Pcor <- numeric(M)
  for (mi in seq_len(M)) {
    set.seed(seed + 1000L + mi)                       # plug-in point set
    fp <- gsimp_impute(bnds, iters_all = iters_all, imp_model = "copula",
                       margin_draw = FALSE)
    Pmeans[mi, ] <- colMeans(fp); Pvars[mi, ] <- apply(fp, 2, stats::var) / n_
    Pmeds[mi, ] <- apply(fp, 2, stats::median)
    Pcor[mi] <- if (p_ > 1) mean(abs((stats::cor(fp) - ctrue)[off])) else NA_real_
    set.seed(seed + 2000L + mi)                       # drawn variance set
    fd <- gsimp_impute(bnds, iters_all = iters_all, imp_model = "copula",
                       margin_draw = TRUE)
    Dmeans[mi, ] <- colMeans(fd)
  }
  qbar <- colMeans(Pmeans); ubar <- colMeans(Pvars); B <- apply(Dmeans, 2, stats::var)
  Tv <- ubar + (1 + 1 / M) * B
  df <- (M - 1) * (1 + ubar / ((1 + 1 / M) * pmax(B, .Machine$double.eps)))^2
  df[!is.finite(df)] <- M - 1
  half <- stats::qt(0.975, df) * sqrt(Tv)

  f_nv <- naive_fill(bnds, "limit"); m_nv <- metrics_point(f_nv, z)

  data.frame(
    scenario = s, rep = r,
    n = cond$n, p = cond$p, rho = cond$rho,
    nd_frac = cond$nd_frac, dnq_frac = cond$dnq_frac, skew = cond$skew,
    censored_frac = mean(bnds$to_impute),
    gs_mean_bias = mean(qbar - truth),
    gs_mean_absbias = mean(abs(qbar - truth)),
    gs_q50_absbias = mean(abs(colMeans(Pmeds) - truth_med)),
    gs_cor_mae = mean(Pcor),
    naive_mean_bias = m_nv$mean_bias,
    naive_mean_absbias = m_nv$mean_absbias,
    naive_cor_mae = m_nv$cor_mae,
    mi_coverage = mean(truth >= qbar - half & truth <= qbar + half),
    mi_ci_width = mean(2 * half)
  )
}

# ---- shuffle + chunk (even load, progress/ETA, checkpoints) ------------------
tasks <- expand.grid(rep = seq_len(n_rep), scenario = seq_len(n_scen))
set.seed(20240101L)
tasks <- tasks[sample.int(nrow(tasks)), , drop = FALSE]
task_list <- Map(function(s, r) c(s, r), tasks$scenario, tasks$rep)
chunks <- split(task_list, cut(seq_len(n_task), breaks = min(n_chunks, n_task),
                               labels = FALSE))

dir.create("validation/results", showWarnings = FALSE, recursive = TRUE)
stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
ckpt  <- sprintf("validation/results/copula_capstone_%s.checkpoint.rds", stamp)

t0 <- Sys.time(); acc <- vector("list", length(chunks)); done <- 0L
for (ci in seq_along(chunks)) {
  res_ci <- mclapply(chunks[[ci]], run_one, mc.cores = workers, mc.preschedule = TRUE)
  errs <- vapply(res_ci, function(x) !is.data.frame(x), logical(1))
  if (any(errs)) {
    logmsg("WARNING: %d task(s) failed in chunk %d", sum(errs), ci)
    res_ci <- res_ci[!errs]
  }
  acc[[ci]] <- do.call(rbind, res_ci)
  done <- done + length(chunks[[ci]])
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs")); rate <- done / el
  logmsg("chunk %d/%d | %d/%d (%.1f%%) | %.2f tasks/s | elapsed %.2f h | ETA %.2f h",
         ci, length(chunks), done, n_task, 100 * done / n_task, rate,
         el / 3600, (n_task - done) / rate / 3600)
  saveRDS(do.call(rbind, acc[seq_len(ci)]), ckpt)
}

res  <- do.call(rbind, acc)
res  <- res[order(res$scenario, res$rep), ]
summ <- summarise_validation(res)
payload <- list(config = config, n_rep = n_rep, M = M, imp_model = "copula",
                iters_all = iters_all, workers = workers, mi = "hybrid",
                grid = grid, results = res, summary = summ,
                elapsed_hours = as.numeric(difftime(Sys.time(), t0, units = "hours")))
saveRDS(payload, sprintf("validation/results/copula_capstone_%s.rds", stamp))
saveRDS(payload, "validation/results/copula_capstone_latest.rds")
utils::write.csv(res,  sprintf("validation/results/copula_capstone_%s.csv", stamp), row.names = FALSE)
utils::write.csv(summ, sprintf("validation/results/copula_capstone_%s_summary.csv", stamp), row.names = FALSE)
if (file.exists(ckpt)) file.remove(ckpt)

logmsg("DONE in %.2f h. Wrote copula_capstone_%s.rds (+ latest, csv, summary).",
       payload$elapsed_hours, stamp)
message("\n== Copula capstone summary (bias, RMSE, MI coverage by scenario) ==")
print(summ[, c("n", "p", "rho", "nd_frac", "skew", "censored_frac",
               "gs_bias", "gs_rmse", "naive_bias", "mi_coverage")], digits = 3)
