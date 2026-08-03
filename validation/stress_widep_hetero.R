# =============================================================================
# Wide-p and heteroscedastic stress test  (Phase 3; PHASE3_SKEW_PLAN.md)
# =============================================================================
#
# Copula (imp_model = "copula", validated in G4/G5) is robust to skew + broken
# margins. This driver stresses the two axes G4/G5 did NOT vary, both under
# strong right-skew, comparing tobit vs the copula hybrid (plug-in point +
# drawn-margin variance):
#
#   WIDE p : p in {5,15,50,150} at n=150 -- p reaches n, so the latent tobit
#            reduces predictors by PCA and the copula fits many per-analyte
#            margins. Does bias/coverage hold as p >> n?
#   HETERO : (a) marginal-scale heterogeneity -- analyte sd spans 0.4..2.0;
#            (b) conditional heteroscedasticity -- a shared volatility factor
#            scales every analyte per row, breaking the homoscedastic-latent
#            assumption both models rely on.
#
# The copula standardises each analyte through its own margin, so it SHOULD shrug
# off marginal-scale heterogeneity; conditional heteroscedasticity is the real
# stress (like `nonlin`, a conditional-model break).
#
# Parallel; bit-identical to serial; BLAS pinned per worker.
#   REPS=40 M=20 WORKERS=22 Rscript validation/stress_widep_hetero.R
# =============================================================================

suppressPackageStartupMessages(library(parallel))
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
source("validation/mc_validation.R", local = TRUE)

n        <- as.integer(Sys.getenv("N", "150"))
skew     <- as.numeric(Sys.getenv("SKEW", "0.75"))
rho      <- as.numeric(Sys.getenv("RHO", "0.5"))
nd       <- as.numeric(Sys.getenv("ND", "0.35"))
methods  <- trimws(strsplit(Sys.getenv("METHODS", "tobit,copula"), ",")[[1]])
reps     <- as.integer(Sys.getenv("REPS", "30"))
iters    <- as.integer(Sys.getenv("ITERS", "20"))
M        <- as.integer(Sys.getenv("M", "15"))
BIAS_TOL <- as.numeric(Sys.getenv("BIAS_TOL", "0.10"))
COV_TOL  <- as.numeric(Sys.getenv("COV_TOL", "0.90"))
workers  <- as.integer(Sys.getenv("WORKERS", as.character(max(1L, detectCores() - 2L))))

# scenarios: kind + p (sd handled inside the generator)
scen <- rbind(
  data.frame(name = sprintf("widep_p%d", c(5,15,50,150)), kind = "homo",
             p = c(5L,15L,50L,150L)),
  data.frame(name = "hetero_marg", kind = "hetero_marg", p = 8L),
  data.frame(name = "hetero_cond", kind = "hetero_cond", p = 8L)
)

gen_data <- function(kind, p) {
  R <- matrix(rho, p, p); diag(R) <- 1
  z <- matrix(stats::rnorm(n * p), n, p) %*% chol(R)
  if (kind == "hetero_cond")                       # shared per-row volatility
    z <- z * exp(0.5 * stats::rnorm(n))            # (recycles down columns: per row)
  if (skew != 0) z <- sinh(asinh(z) + skew)
  sd_vec <- switch(kind,
    homo        = rep(1, p),
    hetero_marg = seq(0.4, 2.0, length.out = p),
    hetero_cond = rep(1, p))
  mu <- seq(0, 1.5, length.out = p)
  sweep(sweep(z, 2, sd_vec, "*"), 2, mu, "+")
}

# hybrid multiple imputation (plug-in point + drawn-margin variance for copula)
mi_metrics <- function(b, z, mo, seed) {
  p <- ncol(z); truth <- colMeans(z); truth_med <- apply(z, 2, stats::median)
  is_cop <- identical(mo, "copula")
  means <- vars <- meds <- matrix(NA_real_, M, p); dmeans <- matrix(NA_real_, M, p)
  for (mi in seq_len(M)) {
    set.seed(seed + mi)
    fm <- gsimp_impute(b, iters_all = iters, imp_model = mo, margin_draw = FALSE)
    means[mi, ] <- colMeans(fm); vars[mi, ] <- apply(fm, 2, stats::var) / n
    meds[mi, ] <- apply(fm, 2, stats::median)
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
  list(bias = mean(qbar - truth), med_bias = mean(colMeans(meds) - truth_med),
       coverage = mean(truth >= qbar - half & truth <= qbar + half))
}

run_one <- function(task) {
  s <- scen[task$si, ]
  set.seed(90000L + task$si * 1000L + task$rep)
  z <- gen_data(s$kind, s$p)
  b <- build_bounds(censor_three_tier(z, nd, 0))
  do.call(rbind, lapply(methods, function(mo) {
    m <- tryCatch(mi_metrics(b, z, mo, seed = 90000L + task$si * 1000L + task$rep + 7L),
                  error = function(e) NULL)
    if (is.null(m)) return(NULL)
    data.frame(name = s$name, kind = s$kind, p = s$p, method = mo,
               rep = task$rep, bias = m$bias, med_bias = m$med_bias,
               coverage = m$coverage)
  }))
}

tasks <- list()
for (si in seq_len(nrow(scen))) for (r in seq_len(reps))
  tasks[[length(tasks) + 1]] <- list(si = si, rep = r)

message(sprintf("Stress test: n=%d skew=%.2f rho=%.1f ND=%.0f%% | methods={%s} reps=%d M=%d",
                n, skew, rho, 100 * nd, paste(methods, collapse = ","), reps, M))
message(sprintf("  scenarios: %s | %d tasks, %d workers",
                paste(scen$name, collapse = ", "), length(tasks), workers))

t0 <- Sys.time()
res <- do.call(rbind, mclapply(tasks, run_one, mc.cores = workers, mc.preschedule = TRUE))
message(sprintf("  done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

agg <- do.call(rbind, lapply(
  split(res, interaction(res$name, res$method, drop = TRUE)),
  function(d) data.frame(
    name = d$name[1], kind = d$kind[1], p = d$p[1], method = d$method[1],
    n_rep = nrow(d), mean_bias = mean(d$bias), abs_mean_bias = abs(mean(d$bias)),
    rmse = sqrt(mean(d$bias^2)), med_absbias = mean(abs(d$med_bias)),
    coverage = mean(d$coverage),
    meets_target = abs(mean(d$bias)) <= BIAS_TOL & mean(d$coverage) >= COV_TOL,
    row.names = NULL)))
agg <- agg[order(match(agg$name, scen$name), match(agg$method, methods)), ]

dir.create("validation/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(scen = scen, per_rep = res, summary = agg),
        "validation/results/stress_widep_hetero.rds")
utils::write.csv(agg, "validation/results/stress_widep_hetero.csv", row.names = FALSE)

message("\n== Wide-p / heteroscedastic stress (skew=", skew, ", ND=", 100*nd, "%) ==")
print(agg[, c("name", "p", "method", "mean_bias", "rmse", "med_absbias",
              "coverage", "meets_target")], row.names = FALSE, digits = 3)
