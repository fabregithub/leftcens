# =============================================================================
# E5 - the three-tier DNQ (detected-not-quantified) band
# =============================================================================
#
# Almost all prior experiments used pure left-censoring (dnq_frac = 0). This one
# exercises the package's actual three-tier model: at a fixed non-detect
# fraction near the target edge (ND = 40%), widen the DNQ band -- the
# interval-censored (MDL, LCMRL) tier -- and see whether recovery holds.
#
# Two things to watch:
#   * DNQ cells are interval-censored (both bounds finite), so *more* informative
#     than non-detects. Recovery should hold at least as well as pure ND.
#   * As DNQ widens, total censored (nd + dnq) crosses 50%, which pushes the
#     MEDIAN into the interval-censored DNQ band even though ND stays < 50%. Does
#     the median stay recoverable when it is no longer directly observed?
#
# The target is defined on the NON-DETECT fraction (ND < 50%); DNQ values are
# "detected" (above MDL). This asks whether the ND<50% rescue is robust to, or
# even helped by, DNQ mass.
#
# Usage:
#   REPS=40 ITERS=20 Rscript validation/effect_of_dnq.R
# =============================================================================

suppressPackageStartupMessages(library(leftcens))
source("validation/mc_validation.R", local = TRUE)

n       <- as.integer(Sys.getenv("N", "150"))
p       <- as.integer(Sys.getenv("P", "6"))
rho     <- as.numeric(Sys.getenv("RHO", "0.5"))
nd      <- as.numeric(Sys.getenv("ND", "0.40"))            # fixed non-detect frac
dnqs    <- as.numeric(strsplit(Sys.getenv("DNQS", "0,0.1,0.2,0.3"), ",")[[1]])
reps    <- as.integer(Sys.getenv("REPS", "30"))
iters   <- as.integer(Sys.getenv("ITERS", "20"))
models  <- trimws(strsplit(Sys.getenv("MODELS", "ridge,tobit"), ",")[[1]])
BIAS_TOL <- as.numeric(Sys.getenv("BIAS_TOL", "0.10"))

message(sprintf("E5 DNQ band: n=%d p=%d rho=%.1f | ND fixed at %.0f%%, DNQ in {%s}",
                n, p, rho, 100 * nd, paste(100 * dnqs, collapse = ",")))
message(sprintf("models=%s reps=%d iters=%d", paste(models, collapse = ","),
                reps, iters))

rows <- list()
for (dnq in dnqs) {
  for (model in models) {
    mb <- medb <- numeric(reps)
    for (r in seq_len(reps)) {
      set.seed(4200L + as.integer(dnq * 100) * 100L + r)  # paired across models
      z <- simulate_truth(n, p, rho)
      b <- build_bounds(censor_three_tier(z, nd, dnq))
      f <- gsimp_impute(b, iters_all = iters, imp_model = model)
      mb[r] <- mean(colMeans(f) - colMeans(z))
      medb[r] <- mean(vapply(seq_len(p),
                             function(j) stats::median(f[, j]) -
                               stats::median(z[, j]), numeric(1)))
    }
    rows[[length(rows) + 1]] <- data.frame(
      model = model, nd = nd, dnq = dnq,
      total_censored = nd + dnq,
      median_in_dnq = (nd < 0.5) & (nd + dnq > 0.5),  # median falls in DNQ band
      mean_bias = mean(mb), rmse = sqrt(mean(mb^2)),
      median_bias = mean(medb),
      meets_target = abs(mean(mb)) <= BIAS_TOL
    )
    message(sprintf("  DNQ=%.0f%% %-6s total_cens=%.0f%% mean_bias=%+.3f median_bias=%+.3f%s",
                    100 * dnq, model, 100 * (nd + dnq), mean(mb), mean(medb),
                    if ((nd < 0.5) && (nd + dnq > 0.5)) "  [median in DNQ band]" else ""))
  }
}
res <- do.call(rbind, rows)

dir.create("validation/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(res, "validation/results/effect_of_dnq.rds")
utils::write.csv(res, "validation/results/effect_of_dnq.csv", row.names = FALSE)

message("\n== E5: three-tier DNQ band (ND fixed at ", 100 * nd, "%) ==")
print(res[, c("model", "dnq", "total_censored", "median_in_dnq",
              "mean_bias", "median_bias", "meets_target")],
      row.names = FALSE, digits = 3)
message("\nWatch median_bias where median_in_dnq = TRUE: does the three-tier")
message("machinery still recover the median once it sits in the DNQ interval?")
