## Generate the bundled example datasets:
##   groundwater  -- three-tier .left/.right layout
##   surfacewater -- two-tier value + qualifier-flag layout
## Reproducible; run with:  source("data-raw/sample_data.R")

set.seed(20240101)

metals <- c("As", "Cd", "Pb", "Cu", "Zn")
mdl   <- c(As = 0.5,  Cd = 0.05, Pb = 0.1, Cu = 1.0, Zn = 2.0)   # detection limit (ug/L)
lcmrl <- c(As = 2.0,  Cd = 0.20, Pb = 0.5, Cu = 5.0, Zn = 10.0)  # quantitation limit (ug/L)
mu    <- c(As = 1.5,  Cd = 0.15, Pb = 0.4, Cu = 4.0, Zn = 12.0)  # median true conc (ug/L)
sdlog <- c(As = 1.1,  Cd = 1.0,  Pb = 1.2, Cu = 0.9, Zn = 0.8)

sim_conc <- function(n) {
  sapply(metals, function(m) stats::rlnorm(n, meanlog = log(mu[[m]]),
                                           sdlog = sdlog[[m]]))
}

## ---- groundwater: three-tier left/right ----------------------------------
n_gw <- 28
conc_gw <- sim_conc(n_gw)
groundwater <- data.frame(sample_id = sprintf("GW-%03d", seq_len(n_gw)),
                          stringsAsFactors = FALSE)
for (m in metals) {
  v <- conc_gw[, m]
  left  <- ifelse(v < mdl[[m]], 0,        ifelse(v < lcmrl[[m]], mdl[[m]],   v))
  right <- ifelse(v < mdl[[m]], mdl[[m]], ifelse(v < lcmrl[[m]], lcmrl[[m]], v))
  groundwater[[paste0(m, ".left")]]  <- round(left, 3)
  groundwater[[paste0(m, ".right")]] <- round(right, 3)
}

## ---- surfacewater: two-tier flagged (value + "<" qualifier) ---------------
n_sw <- 30
conc_sw <- sim_conc(n_sw)
surfacewater <- data.frame(sample_id = sprintf("SW-%03d", seq_len(n_sw)),
                           stringsAsFactors = FALSE)
for (m in metals) {
  v <- conc_sw[, m]
  nd <- v < mdl[[m]]                       # non-detect: report the limit with "<"
  surfacewater[[m]] <- ifelse(nd, mdl[[m]], round(v, 3))
  surfacewater[[paste0(m, ".flag")]] <- ifelse(nd, "<", "")
}

usethis::use_data(groundwater, surfacewater, overwrite = TRUE)
