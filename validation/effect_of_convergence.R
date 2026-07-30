# =============================================================================
# E6 - convergence: how many Gibbs sweeps (iters_all) and imputations (M)?
# =============================================================================
#
# Produces the practical tuning recommendations for the guidance (D2), using the
# default tobit model at a near-target condition (ND = 40%).
#
#   Part A - iters_all: point-recovery (mean bias, RMSE) vs number of Gibbs
#            sweeps. Find where it plateaus -> recommended iters_all.
#   Part B - M: MI coverage and mean 95% interval width vs number of
#            imputations, at a converged iters_all. Find where they stabilise
#            -> recommended M.
#
# Usage:
#   REPS=40 Rscript validation/effect_of_convergence.R
# =============================================================================

suppressPackageStartupMessages(library(leftcens))
source("validation/mc_validation.R", local = TRUE)

n      <- as.integer(Sys.getenv("N", "150"))
p      <- as.integer(Sys.getenv("P", "6"))
rho    <- as.numeric(Sys.getenv("RHO", "0.5"))
nd     <- as.numeric(Sys.getenv("ND", "0.40"))
reps   <- as.integer(Sys.getenv("REPS", "30"))
model  <- Sys.getenv("MODEL", "tobit")
iters_list <- as.integer(strsplit(Sys.getenv("ITERS_LIST", "1,2,3,5,8,12,20,30"),
                                  ",")[[1]])
m_list <- as.integer(strsplit(Sys.getenv("M_LIST", "3,5,10,20,40"), ",")[[1]])
iters_fixed <- as.integer(Sys.getenv("ITERS_FIXED", "20"))

message(sprintf("E6 convergence: model=%s n=%d p=%d ND=%.0f%% reps=%d",
                model, n, p, 100 * nd, reps))

# ---- Part A: sweeps (iters_all) -------------------------------------------
message("\n-- Part A: point recovery vs iters_all --")
rowsA <- list()
for (it in iters_list) {
  mb <- numeric(reps)
  for (r in seq_len(reps)) {
    set.seed(3100L + r)                 # same data across iters levels (paired)
    z <- simulate_truth(n, p, rho)
    b <- build_bounds(censor_three_tier(z, nd, 0))
    f <- gsimp_impute(b, iters_all = it, imp_model = model)
    mb[r] <- mean(colMeans(f) - colMeans(z))
  }
  rowsA[[length(rowsA) + 1]] <- data.frame(
    iters_all = it, mean_bias = mean(mb), rmse = sqrt(mean(mb^2)))
  message(sprintf("  iters_all=%-3d mean_bias=%+.3f rmse=%.3f",
                  it, mean(mb), sqrt(mean(mb^2))))
}
resA <- do.call(rbind, rowsA)

# ---- Part B: imputations (M) ----------------------------------------------
message(sprintf("\n-- Part B: MI coverage / width vs M (iters_all=%d) --", iters_fixed))
rowsB <- list()
for (M in m_list) {
  cov <- wid <- numeric(reps)
  for (r in seq_len(reps)) {
    set.seed(3200L + r)                 # same data across M levels (paired)
    z <- simulate_truth(n, p, rho)
    b <- build_bounds(censor_three_tier(z, nd, 0))
    cc <- mi_coverage(b, z, M = M, iters_all = iters_fixed,
                      imp_model = model, seed = 7000L + r)
    cov[r] <- cc$coverage
    wid[r] <- cc$ci_width
  }
  rowsB[[length(rowsB) + 1]] <- data.frame(
    M = M, coverage = mean(cov), ci_width = mean(wid))
  message(sprintf("  M=%-3d coverage=%.2f mean_ci_width=%.3f",
                  M, mean(cov), mean(wid)))
}
resB <- do.call(rbind, rowsB)

dir.create("validation/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(iters = resA, M = resB, config = list(n = n, p = p, rho = rho,
             nd = nd, reps = reps, model = model, iters_fixed = iters_fixed)),
        "validation/results/effect_of_convergence.rds")
utils::write.csv(resA, "validation/results/convergence_iters.csv", row.names = FALSE)
utils::write.csv(resB, "validation/results/convergence_M.csv", row.names = FALSE)

message("\n== E6 Part A: point recovery vs iters_all ==")
print(resA, row.names = FALSE, digits = 3)
message("== E6 Part B: MI coverage / width vs M ==")
print(resB, row.names = FALSE, digits = 3)
message("\nRecommended iters_all = smallest where bias/RMSE have plateaued;")
message("recommended M = smallest where coverage and interval width have settled.")
