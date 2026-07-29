# =============================================================================
# ridge vs tobit conditional model: does the censored model extend the
# reliable detection-rate range?
# =============================================================================
#
# The detection-rate standard was derived with the observed-only ridge model,
# which is biased under heavy censoring (coverage collapses). The tobit model
# (imp_model = "tobit") fits an interval-censored Gaussian regression on both
# observed and censored rows, removing that selection bias. This script runs the
# SAME simulated datasets through both models and compares bias and MI coverage
# across the detection-rate sweep, at two sample sizes (coverage failure was
# worst at large n), to see how far tobit pushes the reliable threshold.
#
# Reuses run_validation()/summarise_validation() from mc_validation.R, calling
# each model with the same base_seed so the comparison is paired (identical
# data and RNG streams, only the model differs).
#
# Usage:
#   Rscript validation/tobit_vs_ridge.R
#   N_REP=100 M=25 Rscript validation/tobit_vs_ridge.R
# =============================================================================

suppressPackageStartupMessages(library(leftcens))
source("validation/mc_validation.R", local = TRUE)

n_rep     <- as.integer(Sys.getenv("N_REP", "20"))
M         <- as.integer(Sys.getenv("M", "10"))
iters_all <- as.integer(Sys.getenv("ITERS_ALL", "20"))
BIAS_TOL     <- as.numeric(Sys.getenv("BIAS_TOL", "0.10"))
COVERAGE_TOL <- as.numeric(Sys.getenv("COVERAGE_TOL", "0.90"))
base_seed <- 51515L

# Detection sweep x sample size; pure left-censoring, fixed p / rho / no skew
# (the symmetric, worst-case regime where ridge failed).
grid <- expand.grid(
  n = c(100L, 300L),
  p = 6L,
  rho = 0.5,
  nd_frac = seq(0.15, 0.55, by = 0.10),
  dnq_frac = 0,
  skew = 0,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)

message(sprintf("ridge vs tobit: %d scenarios x %d reps (M=%d, iters=%d)",
                nrow(grid), n_rep, M, iters_all))

run_model <- function(model) {
  s <- summarise_validation(
    run_validation(grid, n_rep = n_rep, M = M, iters_all = iters_all,
                   imp_model = model, base_seed = base_seed, verbose = FALSE)
  )
  s$detection_rate <- 1 - s$nd_frac
  s[order(s$n, -s$detection_rate), ]
}

message("-- ridge --"); sr <- run_model("ridge")
message("-- tobit --"); st <- run_model("tobit")

key <- c("n", "detection_rate", "nd_frac")
cmp <- merge(
  sr[, c(key, "censored_frac", "gs_bias", "mi_coverage")],
  st[, c(key, "gs_bias", "mi_coverage")],
  by = key, suffixes = c("_ridge", "_tobit")
)
cmp <- cmp[order(cmp$n, -cmp$detection_rate), ]
cmp$ridge_reliable <- abs(cmp$gs_bias_ridge) <= BIAS_TOL &
  cmp$mi_coverage_ridge >= COVERAGE_TOL
cmp$tobit_reliable <- abs(cmp$gs_bias_tobit) <= BIAS_TOL &
  cmp$mi_coverage_tobit >= COVERAGE_TOL

# Lowest reliable detection rate per model, per sample size.
threshold <- do.call(rbind, lapply(split(cmp, cmp$n), function(d) {
  low <- function(ok) {
    d2 <- d[order(-d$detection_rate), ]
    keep <- cumprod(ok[order(-d$detection_rate)]) == 1
    if (!any(keep)) NA_real_ else min(d2$detection_rate[keep])
  }
  data.frame(n = d$n[1],
             ridge_min_reliable_detection = low(d$ridge_reliable),
             tobit_min_reliable_detection = low(d$tobit_reliable))
}))

dir.create("validation/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(grid = grid, n_rep = n_rep, M = M, iters_all = iters_all,
             comparison = cmp, threshold = threshold),
        "validation/results/tobit_vs_ridge.rds")
utils::write.csv(cmp, "validation/results/tobit_vs_ridge.csv", row.names = FALSE)

message("\n== ridge vs tobit: bias and MI coverage by detection rate and n ==")
print(cmp[, c("n", "detection_rate", "censored_frac",
              "gs_bias_ridge", "mi_coverage_ridge",
              "gs_bias_tobit", "mi_coverage_tobit",
              "ridge_reliable", "tobit_reliable")],
      row.names = FALSE, digits = 3)

message("\n== Lowest reliable detection rate (|bias|<=", BIAS_TOL,
        " & coverage>=", COVERAGE_TOL, ") ==")
print(threshold, row.names = FALSE, digits = 3)
