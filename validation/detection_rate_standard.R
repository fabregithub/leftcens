# =============================================================================
# Detection-rate reliability standard for gsimp_impute()
# =============================================================================
#
# Sweeps the non-detect fraction (pure left-censoring, dnq_frac = 0, so the
# "detection rate" = 1 - nd_frac is unambiguous) crossed with correlation
# strength, and identifies the lowest detection rate at which imputation is
# still "reliable" under two criteria:
#
#   * |mean bias| <= BIAS_TOL on the log scale, AND
#   * MI 95% coverage >= COVERAGE_TOL
#
# Reliability depends on more than censoring alone, so the standard is reported
# CONDITIONAL on correlation strength.
#
# Usage:
#   Rscript validation/detection_rate_standard.R
#   N_REP=120 M=20 Rscript validation/detection_rate_standard.R
# =============================================================================

suppressPackageStartupMessages(library(leftcens))
source("validation/mc_validation.R", local = TRUE)  # defines the engine

BIAS_TOL     <- as.numeric(Sys.getenv("BIAS_TOL", "0.10"))
COVERAGE_TOL <- as.numeric(Sys.getenv("COVERAGE_TOL", "0.90"))
n_rep        <- as.integer(Sys.getenv("N_REP", "80"))
M            <- as.integer(Sys.getenv("M", "15"))
iters_all    <- as.integer(Sys.getenv("ITERS_ALL", "30"))

# Detection-rate sweep x correlation. dnq_frac = 0 -> pure left-censoring.
# Grid resolution is configurable so the study can be run at publication
# resolution or a lighter foreground pass.
ndf_by <- as.numeric(Sys.getenv("NDF_BY", "0.05"))
rhos <- as.numeric(strsplit(Sys.getenv("RHOS", "0,0.3,0.6"), ",")[[1]])
grid <- expand.grid(
  n = 200L, p = 6L,
  rho = rhos,
  nd_frac = seq(0.10, 0.70, by = ndf_by),
  dnq_frac = 0,
  skew = 0,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)

message(sprintf("Detection-rate sweep: %d scenarios x %d reps (M=%d, iters=%d)",
                nrow(grid), n_rep, M, iters_all))
message(sprintf("Reliability: |bias| <= %.2f AND coverage >= %.2f",
                BIAS_TOL, COVERAGE_TOL))

res <- run_validation(grid, n_rep = n_rep, M = M, iters_all = iters_all,
                      imp_model = "ridge", base_seed = 424242L)
summ <- summarise_validation(res)

summ$detection_rate <- 1 - summ$nd_frac
summ$reliable <- abs(summ$gs_bias) <= BIAS_TOL & summ$mi_coverage >= COVERAGE_TOL
summ <- summ[order(summ$rho, -summ$detection_rate), ]

# For each correlation level, the lowest detection rate that is still reliable,
# treating reliability as holding down to the first breach as detection falls.
threshold <- do.call(rbind, lapply(split(summ, summ$rho), function(d) {
  d <- d[order(-d$detection_rate), ]           # high detection -> low
  ok <- d$reliable
  # lowest detection rate before the first FALSE (contiguous reliable region)
  last_ok <- if (!ok[1]) NA_real_ else d$detection_rate[which(cumprod(ok) == 1)]
  data.frame(
    rho = d$rho[1],
    min_reliable_detection_rate = if (all(is.na(last_ok))) NA else min(last_ok),
    max_reliable_nd_frac = if (all(is.na(last_ok))) NA else 1 - min(last_ok)
  )
}))

dir.create("validation/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(grid = grid, results = res, summary = summ, threshold = threshold,
             bias_tol = BIAS_TOL, coverage_tol = COVERAGE_TOL),
        "validation/results/detection_rate_standard.rds")
utils::write.csv(summ, "validation/results/detection_rate_summary.csv",
                 row.names = FALSE)

message("\n== Reliability by detection rate and correlation ==")
print(summ[, c("rho", "detection_rate", "nd_frac", "censored_frac",
               "gs_bias", "gs_rmse", "mi_coverage", "reliable")],
      row.names = FALSE, digits = 3)

message("\n== Rough standard: lowest reliable detection rate by correlation ==")
print(threshold, row.names = FALSE, digits = 3)
