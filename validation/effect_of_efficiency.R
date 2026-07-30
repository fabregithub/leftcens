# =============================================================================
# E4 - efficiency: does correlation shrink tobit's intervals (not its bias)?
# =============================================================================
# E1/E3 showed correlation does not change tobit's bias/coverage (the rescue is
# per-analyte). This checks the predicted flip side: informative correlation
# should improve *efficiency* -- narrower MI intervals -- even when bias is
# already ~0. tobit only; small grid to keep it cheap.
#
#   REPS=30 Rscript validation/effect_of_efficiency.R
# =============================================================================

suppressPackageStartupMessages(library(leftcens))
source("validation/mc_validation.R", local = TRUE)

n     <- as.integer(Sys.getenv("N", "150"))
p     <- as.integer(Sys.getenv("P", "6"))
reps  <- as.integer(Sys.getenv("REPS", "15"))
M     <- as.integer(Sys.getenv("M", "10"))
iters <- as.integer(Sys.getenv("ITERS", "10"))
nds   <- as.numeric(strsplit(Sys.getenv("NDS", "0.25,0.40"), ",")[[1]])
rhos  <- as.numeric(strsplit(Sys.getenv("RHOS", "0,0.5,0.8"), ",")[[1]])

message(sprintf("E4 efficiency (tobit): n=%d p=%d | ND in {%s} x rho in {%s} | reps=%d M=%d",
                n, p, paste(nds, collapse = ","), paste(rhos, collapse = ","),
                reps, M))
rows <- list()
for (nd in nds) {
  for (rho in rhos) {
    bias <- cov <- wid <- numeric(reps)
    for (r in seq_len(reps)) {
      set.seed(8800L + as.integer(nd * 100) * 1000L + as.integer(rho * 100) + r)
      z <- simulate_truth(n, p, rho)
      b <- build_bounds(censor_three_tier(z, nd, 0))
      f <- gsimp_impute(b, iters_all = iters, imp_model = "tobit")
      bias[r] <- mean(colMeans(f) - colMeans(z))
      cc <- mi_coverage(b, z, M = M, iters_all = iters, imp_model = "tobit",
                        seed = 9900L + r)
      cov[r] <- cc$coverage
      wid[r] <- cc$ci_width
    }
    rows[[length(rows) + 1]] <- data.frame(
      nd = nd, rho = rho, mean_bias = mean(bias),
      coverage = mean(cov), ci_width = mean(wid))
    message(sprintf("  ND=%.0f%% rho=%.1f bias=%+.3f coverage=%.2f ci_width=%.3f",
                    100 * nd, rho, mean(bias), mean(cov), mean(wid)))
  }
}
res <- do.call(rbind, rows)
dir.create("validation/results", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(res, "validation/results/effect_of_efficiency.csv", row.names = FALSE)
message("\n== E4: tobit efficiency vs correlation ==")
print(res, row.names = FALSE, digits = 3)
message("\nIf ci_width falls as rho rises (at ~flat bias/coverage), correlation")
message("improves efficiency without affecting validity -- as E1/E3 predicted.")
