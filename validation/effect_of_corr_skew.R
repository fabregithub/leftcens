# =============================================================================
# G3 - correlation under skew: does dependence now HELP?  (PHASE3_SKEW_PLAN.md)
# =============================================================================
#
# Earlier work (E1/E4) found correlation does NOT help tobit recover the marginal
# mean -- but that was a property of the misspecified linear-Gaussian model. The
# Gaussian copula predicts the opposite: the latent conditional variance
# 1 - Sigma_{j,-j} Sigma^{-1} Sigma_{-j,j} SHRINKS as correlation grows, so an
# observed correlated analyte pins a censored latent value into a narrower band
# -> less tail extrapolation -> lower bias, RMSE and interval width. That effect
# should be largest at heavy censoring, where extrapolation matters most.
#
# This driver sweeps rho x ND at strong skew (0.75), comparing tobit vs the
# copula hybrid (plug-in point + drawn-margin variance), on two estimands:
#   (a) the MARGINAL mean  -- does correlation cut bias / RMSE / CI width?
#   (b) a JOINT functional -- recovery of the between-analyte correlation
#       (cor_mae), where dependence should matter even for tobit.
#
# Parallel; bit-identical to serial; BLAS pinned per worker.
#   REPS=100 M=25 WORKERS=22 Rscript validation/effect_of_corr_skew.R
# =============================================================================

suppressPackageStartupMessages(library(parallel))
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
source("validation/mc_validation.R", local = TRUE)

n        <- as.integer(Sys.getenv("N", "150"))
p        <- as.integer(Sys.getenv("P", "6"))
skew     <- as.numeric(Sys.getenv("SKEW", "0.75"))
rhos     <- as.numeric(strsplit(Sys.getenv("RHOS", "0,0.3,0.6,0.9"), ",")[[1]])
nds      <- as.numeric(strsplit(Sys.getenv("NDS", "0.25,0.35,0.45"), ",")[[1]])
methods  <- trimws(strsplit(Sys.getenv("METHODS", "tobit,copula"), ",")[[1]])
reps     <- as.integer(Sys.getenv("REPS", "60"))
iters    <- as.integer(Sys.getenv("ITERS", "20"))
M        <- as.integer(Sys.getenv("M", "20"))
workers  <- as.integer(Sys.getenv("WORKERS", as.character(max(1L, detectCores() - 2L))))

mi_metrics <- function(b, z, mo, seed) {
  p <- ncol(z); truth <- colMeans(z); ctrue <- stats::cor(z); off <- upper.tri(ctrue)
  is_cop <- identical(mo, "copula")
  means <- vars <- matrix(NA_real_, M, p); dmeans <- matrix(NA_real_, M, p)
  cormae <- numeric(M)
  for (mi in seq_len(M)) {
    set.seed(seed + mi)
    fm <- gsimp_impute(b, iters_all = iters, imp_model = mo, margin_draw = FALSE)
    means[mi, ] <- colMeans(fm); vars[mi, ] <- apply(fm, 2, stats::var) / n
    cormae[mi] <- mean(abs((stats::cor(fm) - ctrue)[off]))
    if (is_cop) {
      set.seed(seed + 500L + mi)
      fd <- gsimp_impute(b, iters_all = iters, imp_model = mo, margin_draw = TRUE)
      dmeans[mi, ] <- colMeans(fd)
    }
  }
  qbar <- colMeans(means); ubar <- colMeans(vars)
  bv <- if (is_cop) apply(dmeans, 2, stats::var) else apply(means, 2, stats::var)
  Tv <- ubar + (1 + 1 / M) * bv
  df <- (M - 1) * (1 + ubar / ((1 + 1 / M) * pmax(bv, .Machine$double.eps)))^2
  df[!is.finite(df)] <- M - 1
  half <- stats::qt(0.975, df) * sqrt(Tv)
  list(bias = mean(qbar - truth), coverage = mean(truth >= qbar - half & truth <= qbar + half),
       ci_width = mean(2 * half), cor_mae = mean(cormae))
}

grid <- expand.grid(rho = rhos, nd = nds, KEEP.OUT.ATTRS = FALSE)
run_one <- function(task) {
  g <- grid[task$gi, ]
  set.seed(64000L + task$gi * 1000L + task$rep)
  z <- simulate_truth(n, p, g$rho, skew = skew)
  b <- build_bounds(censor_three_tier(z, g$nd, 0))
  do.call(rbind, lapply(methods, function(mo) {
    m <- tryCatch(mi_metrics(b, z, mo, seed = 64000L + task$gi * 1000L + task$rep + 5L),
                  error = function(e) NULL)
    if (is.null(m)) return(NULL)
    data.frame(rho = g$rho, nd = g$nd, method = mo, rep = task$rep,
               bias = m$bias, coverage = m$coverage, ci_width = m$ci_width,
               cor_mae = m$cor_mae)
  }))
}

tasks <- list()
for (gi in seq_len(nrow(grid))) for (r in seq_len(reps))
  tasks[[length(tasks) + 1]] <- list(gi = gi, rep = r)

message(sprintf("G3 corr-under-skew: n=%d p=%d skew=%.2f | rho={%s} x ND={%s}",
                n, p, skew, paste(rhos, collapse = ","), paste(100 * nds, collapse = ",")))
message(sprintf("  methods={%s} reps=%d M=%d | %d tasks, %d workers",
                paste(methods, collapse = ","), reps, M, length(tasks), workers))

t0 <- Sys.time()
res <- do.call(rbind, mclapply(tasks, run_one, mc.cores = workers, mc.preschedule = TRUE))
message(sprintf("  done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

agg <- do.call(rbind, lapply(
  split(res, interaction(res$rho, res$nd, res$method, drop = TRUE)),
  function(d) data.frame(
    rho = d$rho[1], nd = d$nd[1], method = d$method[1], n_rep = nrow(d),
    mean_bias = mean(d$bias), abs_mean_bias = abs(mean(d$bias)),
    rmse = sqrt(mean(d$bias^2)), coverage = mean(d$coverage),
    ci_width = mean(d$ci_width), cor_mae = mean(d$cor_mae), row.names = NULL)))
agg <- agg[order(agg$nd, match(agg$method, methods), agg$rho), ]

dir.create("validation/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(grid = grid, per_rep = res, summary = agg),
        "validation/results/effect_of_corr_skew.rds")
utils::write.csv(agg, "validation/results/effect_of_corr_skew.csv", row.names = FALSE)

message("\n== G3: correlation under skew=", skew, " (rmse & CI width = efficiency) ==")
print(agg[, c("nd", "method", "rho", "mean_bias", "rmse", "ci_width",
              "coverage", "cor_mae")], row.names = FALSE, digits = 3)

message("\nRead: for COPULA, does rmse / ci_width / |bias| fall as rho rises (esp. at ND 45%)?")
message("For the JOINT estimand, does cor_mae fall with rho (copula < tobit)?")
