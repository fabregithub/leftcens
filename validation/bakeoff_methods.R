# =============================================================================
# G4 / G5 - method bake-off across skew and dependence families
#   (Phase 3; see validation/PHASE3_SKEW_PLAN.md)
# =============================================================================
#
# G4: compare imputation methods {naive, ridge, tobit, copula} head-to-head on
# IDENTICAL datasets (paired), scoring mean bias, median bias, RMSE, MI coverage.
#
# G5 (the honest part): the copula model assumes a Gaussian copula with
# sinh-arcsinh margins. To avoid an inverse crime we generate from families that
# break each assumption in turn, and record the REALIZED marginal skewness /
# linear correlation so methods are compared on what the data actually are:
#
#   family   margin            dependence         breaks...
#   ------   ------            ----------         --------
#   sas      sinh-arcsinh      Gaussian copula    nothing (matched -> ceiling)
#   gamma    Gamma             Gaussian copula    the MARGIN (parametric)
#   lnorm    log-normal        Gaussian copula    the MARGIN (heavy right tail)
#   nonlin   non-linear latent (shared basis)     the COPULA / dependence shape
#
# A method only "wins" if it holds the target (|bias| <= 0.10 AND coverage
# >= 0.90) across ALL families -- especially `nonlin`, which no method here is
# matched to.
#
# Parallel over (family, level, dep, nd, rep) tasks; each task builds one dataset
# and runs every method on it. Bit-identical to serial; BLAS pinned per worker.
#
# Usage (defaults ~30-40 min on 22 cores; scale up on the workstation):
#   ./validation/... nohup style, or:
#   REPS=150 M=25 WORKERS=22 Rscript validation/bakeoff_methods.R
#   FAMILIES=sas,nonlin METHODS=tobit,copula Rscript validation/bakeoff_methods.R
# =============================================================================

suppressPackageStartupMessages(library(parallel))
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
source("validation/mc_validation.R", local = TRUE)

n        <- as.integer(Sys.getenv("N", "150"))
p        <- as.integer(Sys.getenv("P", "6"))
families <- trimws(strsplit(Sys.getenv("FAMILIES", "sas,gamma,lnorm,nonlin"), ",")[[1]])
levels_  <- trimws(strsplit(Sys.getenv("LEVELS", "mod,strong"), ",")[[1]])
nds      <- as.numeric(strsplit(Sys.getenv("NDS", "0.25,0.35,0.45"), ",")[[1]])
methods  <- trimws(strsplit(Sys.getenv("METHODS", "naive,ridge,tobit,copula"), ",")[[1]])
reps     <- as.integer(Sys.getenv("REPS", "40"))
iters    <- as.integer(Sys.getenv("ITERS", "20"))
M        <- as.integer(Sys.getenv("M", "12"))
BIAS_TOL <- as.numeric(Sys.getenv("BIAS_TOL", "0.10"))
COV_TOL  <- as.numeric(Sys.getenv("COV_TOL", "0.90"))
workers  <- as.integer(Sys.getenv("WORKERS", as.character(max(1L, detectCores() - 2L))))
base_seed <- 82000L
mu <- seq(0, 1.5, length.out = p)

# dependence knob values per family: Gaussian-copula rho, or non-linear strength.
dep_of <- function(family) if (family == "nonlin") c(0.4, 0.7) else c(0, 0.6)

# ---- data-generating families (all return an n x p working/log-scale matrix) --
gen_data <- function(family, level, dep, n, p) {
  cholR <- function(rho) { R <- matrix(rho, p, p); diag(R) <- 1; chol(R) }
  switch(family,
    sas = {                                   # sinh-arcsinh margin, Gaussian copula
      skew <- c(mod = 0.5, strong = 0.75)[[level]]
      simulate_truth(n, p, dep, skew = skew)
    },
    gamma = {                                 # Gamma margin, Gaussian copula
      shape <- c(mod = 8, strong = 4)[[level]]
      z <- matrix(stats::rnorm(n * p), n, p) %*% cholR(dep)
      X <- matrix(0, n, p)
      for (j in seq_len(p))
        X[, j] <- mu[j] + (stats::qgamma(stats::pnorm(z[, j]), shape) - shape) /
          sqrt(shape)
      X
    },
    lnorm = {                                 # log-normal margin, Gaussian copula
      s <- c(mod = 0.25, strong = 0.35)[[level]]
      z <- matrix(stats::rnorm(n * p), n, p) %*% cholR(dep)
      m <- exp(s^2 / 2); v <- (exp(s^2) - 1) * exp(s^2)
      X <- matrix(0, n, p)
      for (j in seq_len(p))
        X[, j] <- mu[j] + (exp(s * z[, j]) - m) / sqrt(v)
      X
    },
    nonlin = {                                # non-linear dependence (breaks copula)
      simulate_truth_nl(n, p, strength = dep)
    },
    stop("unknown family: ", family)
  )
}

marg_skew <- function(z) mean(vapply(seq_len(ncol(z)), function(j) {
  x <- z[, j]; m <- mean(x); s <- stats::sd(x)
  if (s <= 0) return(0); mean((x - m)^3) / s^3 }, numeric(1)))
mean_abscor <- function(z) { if (ncol(z) < 2) return(NA_real_)
  C <- stats::cor(z); mean(abs(C[upper.tri(C)])) }

# ---- one dataset, all methods (paired) -------------------------------------
run_one <- function(task) {
  fi <- match(task$family, families); li <- match(task$level, levels_)
  di <- round(task$dep * 100); ni <- round(task$nd * 100)
  seed <- base_seed + fi * 1000000L + li * 100000L + di * 1000L + ni * 10L + task$rep
  set.seed(seed)
  z <- gen_data(task$family, task$level, task$dep, n, p)
  b <- build_bounds(censor_three_tier(z, task$nd, 0))
  sk <- marg_skew(z); cr <- mean_abscor(z)

  truth <- colMeans(z); truth_med <- apply(z, 2, stats::median)
  do.call(rbind, lapply(methods, function(mo) {
    out <- tryCatch({
      if (mo == "naive") {
        f <- naive_fill(b, "limit")
        list(mb = mean(colMeans(f) - truth),
             medb = mean(apply(f, 2, stats::median) - truth_med),
             cov = NA_real_, ciw = NA_real_)
      } else {
        # Proper multiple imputation. For copula use the recommended hybrid: the
        # POINT (qbar, ubar) from plug-in-margin imputations (near-unbiased, no
        # Jensen bias) and the between-imputation variance from DRAWN-margin
        # imputations (keeps the coverage fix). Other models ignore margin_draw,
        # so their point and variance come from the one set.
        is_cop <- identical(mo, "copula")
        means <- matrix(NA_real_, M, p); vars <- matrix(NA_real_, M, p)
        meds <- matrix(NA_real_, M, p); dmeans <- matrix(NA_real_, M, p)
        for (mi in seq_len(M)) {
          set.seed(seed + 1000L + mi)
          fm <- gsimp_impute(b, iters_all = iters, imp_model = mo, margin_draw = FALSE)
          means[mi, ] <- colMeans(fm); vars[mi, ] <- apply(fm, 2, stats::var) / n
          meds[mi, ] <- apply(fm, 2, stats::median)
          if (is_cop) {
            set.seed(seed + 2000L + mi)
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
        list(mb = mean(qbar - truth),
             medb = mean(colMeans(meds) - truth_med),
             cov = mean(truth >= qbar - half & truth <= qbar + half),
             ciw = mean(2 * half))
      }
    }, error = function(e) NULL)
    if (is.null(out)) return(NULL)
    data.frame(family = task$family, level = task$level, dep = task$dep,
               nd = task$nd, rep = task$rep, method = mo,
               skew_real = sk, cor_real = cr,
               mean_bias = out$mb, med_bias = out$medb,
               coverage = out$cov, ci_width = out$ciw)
  }))
}

# ---- task list --------------------------------------------------------------
tasks <- list()
for (fam in families) for (lev in levels_) for (dp in dep_of(fam))
  for (nd in nds) for (r in seq_len(reps))
    tasks[[length(tasks) + 1]] <- list(family = fam, level = lev, dep = dp,
                                       nd = nd, rep = r)

message(sprintf("G4/G5 bake-off: n=%d p=%d | families={%s} levels={%s} nd={%s}",
                n, p, paste(families, collapse = ","),
                paste(levels_, collapse = ","), paste(100 * nds, collapse = ",")))
message(sprintf("  methods={%s} reps=%d M=%d iters=%d | %d tasks, %d workers",
                paste(methods, collapse = ","), reps, M, iters,
                length(tasks), workers))

t0 <- Sys.time()
res_list <- mclapply(tasks, run_one, mc.cores = workers, mc.preschedule = TRUE)
per_rep <- do.call(rbind, res_list[!vapply(res_list, is.null, logical(1))])
message(sprintf("  done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- aggregate --------------------------------------------------------------
agg <- do.call(rbind, lapply(
  split(per_rep, interaction(per_rep$family, per_rep$level, per_rep$dep,
                             per_rep$nd, per_rep$method, drop = TRUE)),
  function(d) {
    mcse <- function(x) stats::sd(x) / sqrt(length(x))
    data.frame(
      family = d$family[1], level = d$level[1], dep = d$dep[1], nd = d$nd[1],
      method = d$method[1], n_rep = nrow(d),
      skew_real = mean(d$skew_real), cor_real = mean(d$cor_real),
      mean_bias = mean(d$mean_bias), abs_mean_bias = abs(mean(d$mean_bias)),
      bias_mcse = mcse(d$mean_bias), rmse = sqrt(mean(d$mean_bias^2)),
      med_absbias = mean(abs(d$med_bias)),
      coverage = mean(d$coverage), ci_width = mean(d$ci_width),
      meets_target = abs(mean(d$mean_bias)) <= BIAS_TOL &
        (is.na(mean(d$coverage)) | mean(d$coverage) >= COV_TOL),
      row.names = NULL)
  }))
agg <- agg[order(agg$family, agg$level, agg$dep, agg$nd,
                 match(agg$method, methods)), ]

dir.create("validation/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(config = list(n = n, p = p, families = families, levels = levels_,
                           nds = nds, methods = methods, reps = reps, M = M,
                           iters = iters, bias_tol = BIAS_TOL, cov_tol = COV_TOL),
             per_rep = per_rep, summary = agg),
        "validation/results/bakeoff_methods.rds")
utils::write.csv(agg, "validation/results/bakeoff_methods.csv", row.names = FALSE)

# ---- report ----------------------------------------------------------------
message("\n== G4/G5 bake-off summary (bias / coverage by family x ND x method) ==")
print(agg[, c("family", "level", "dep", "nd", "method", "skew_real",
              "mean_bias", "rmse", "coverage", "meets_target")],
      row.names = FALSE, digits = 3)

# per-method pass rate, split by whether the family is matched (sas) or not.
message("\n== Target pass-rate by method (|bias|<=", BIAS_TOL,
        " AND coverage>=", COV_TOL, ") ==")
for (mo in methods) {
  d <- agg[agg$method == mo, ]
  matched <- d[d$family == "sas", ]; broken <- d[d$family != "sas", ]
  message(sprintf("  %-7s overall %2d/%2d | matched(sas) %d/%d | broken-family %2d/%2d",
                  mo, sum(d$meets_target), nrow(d),
                  sum(matched$meets_target), nrow(matched),
                  sum(broken$meets_target), nrow(broken)))
}
message("\nDecision: copula must beat tobit on the broken families too (esp. nonlin),")
message("not only on the matched sinh-arcsinh case, to justify becoming a default.")
