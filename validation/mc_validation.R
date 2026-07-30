# =============================================================================
# Monte Carlo validation study for leftcens::gsimp_impute()
# =============================================================================
#
# This is a standalone reproducible study (NOT part of the package test suite or
# build). It quantifies how well the Gibbs-sampler imputation recovers known
# ground truth across a grid of data-generating conditions, and whether its
# multiple-imputation uncertainty is calibrated. Intended as the basis for a
# methods/validation write-up.
#
# What it measures, per condition, averaged over many replications:
#   * point bias and RMSE of recovered analyte means (log scale)
#   * bias of recovered quantiles (10th / 50th / 90th)
#   * recovery of the between-analyte correlation structure
#   * multiple-imputation (Rubin's rules) 95% CI COVERAGE of the true mean
#     -- the key calibration metric
# gsimp is compared against naive substitution at the detection limit.
#
# Conditions varied: sample size n, number of analytes p, correlation rho,
# censoring fraction, and marginal skewness (log-normal vs sinh-arcsinh-skewed),
# so the method is stressed on the axes that matter (weak correlation, heavy
# censoring, and model misspecification).
#
# Usage:
#   # from the package root, with leftcens installed:
#   Rscript validation/mc_validation.R                 # default (quick) config
#   N_REP=500 M=50 CONFIG=full Rscript validation/mc_validation.R
#
#   # or source the functions and drive them yourself:
#   source("validation/mc_validation.R")
#   res <- run_validation(make_grid("full"), n_rep = 500, M = 50)
#
# Results are written to validation/results/ as .rds and .csv.
# =============================================================================

suppressPackageStartupMessages(library(leftcens))

# ---- data-generating process -----------------------------------------------

#' Simulate ground-truth log-concentrations.
#'
#' Multivariate normal on the log scale with exchangeable correlation `rho`.
#' A non-zero `skew` applies a sinh-arcsinh transform to each margin, so the
#' data departs from log-normality (a model-misspecification stress test) while
#' preserving monotone ordering.
simulate_truth <- function(n, p, rho, mu = NULL, sd = 1, skew = 0) {
  if (is.null(mu)) mu <- seq(0, 1.5, length.out = p)
  mu <- rep_len(mu, p)
  R <- matrix(rho, p, p)
  diag(R) <- 1
  z <- matrix(stats::rnorm(n * p), n, p) %*% chol(R)
  if (skew != 0) z <- sinh(asinh(z) + skew)   # sinh-arcsinh skew
  sweep(z * sd, 2, mu, "+")                    # n x p, log scale
}

#' Simulate ground truth with NON-PARAMETRIC (non-linear) dependence.
#'
#' All analytes are driven by a shared standard-normal latent through different
#' *non-linear* basis functions, plus idiosyncratic noise. The default bases are
#' the (orthogonal) Hermite polynomials H2, H3 and the standardised sine and
#' absolute-value maps, none of which is the identity -- so the analytes are
#' strongly dependent (common latent) yet have near-zero pairwise *linear*
#' correlation. This is the stress test for the linear conditional model in
#' `imp_model = "tobit"`: a linear regression cannot see the dependence, whereas
#' a non-parametric model could.
#'
#' Margins are standardised to mean `mu[j]`, sd `sd` (unit by default), matching
#' [simulate_truth()] so censoring at a given quantile yields comparable
#' non-detect fractions.
#'
#' @param n,p,mu,sd As in [simulate_truth()].
#' @param strength Weight (in `[0, 1]`) of the shared latent signal vs.
#'   idiosyncratic noise; analogous to correlation strength. Higher = stronger
#'   (non-linear) dependence.
#' @param forms Optional character vector of basis names cycled across analytes;
#'   any of `"lin"`, `"quad"`, `"cub"`, `"sine"`, `"abs"`. Default cycles the
#'   four non-linear bases (excludes `"lin"`), maximising non-linearity.
#' @return An `n x p` log-scale matrix.
simulate_truth_nl <- function(n, p, strength = 0.6, mu = NULL, sd = 1,
                              forms = NULL) {
  if (is.null(mu)) mu <- seq(0, 1.5, length.out = p)
  mu <- rep_len(mu, p)
  strength <- min(max(strength, 0), 1)

  u <- stats::rnorm(n)                       # shared latent
  std <- function(g) (g - mean(g)) / stats::sd(g)   # empirical unit variance
  basis <- list(
    lin  = function(z) z,                    # identity (linear) -- excluded by default
    quad = function(z) (z^2 - 1) / sqrt(2),  # Hermite H2 (unit var)
    cub  = function(z) (z^3 - 3 * z) / sqrt(6),  # Hermite H3 (unit var)
    sine = function(z) std(sin(2 * z)),
    abs  = function(z) std(abs(z))
  )
  if (is.null(forms)) forms <- c("quad", "cub", "sine", "abs")
  forms <- rep(forms, length.out = p)

  X <- matrix(0, n, p)
  for (j in seq_len(p)) {
    g <- basis[[forms[j]]](u)                # unit-variance non-linear signal
    X[, j] <- mu[j] + sd * (sqrt(strength) * g +
                              sqrt(1 - strength) * stats::rnorm(n))
  }
  X
}

#' Censor a log-scale truth matrix into the three-tier interval structure.
#'
#' Per analyte, `nd_frac` of values fall below the MDL (non-detect) and a further
#' `dnq_frac` between MDL and LCMRL (detected-not-quantified). Returns a list of
#' `cens_data` columns. `nd_frac` and `dnq_frac` may be scalars (applied to every
#' analyte) or length-`p` vectors for per-analyte (heterogeneous) detection
#' limits.
censor_three_tier <- function(z_true, nd_frac, dnq_frac) {
  conc <- exp(z_true)
  p <- ncol(conc)
  nd_frac <- rep_len(nd_frac, p)
  dnq_frac <- rep_len(dnq_frac, p)
  lapply(seq_len(p), function(j) {
    cj <- conc[, j]
    mdl <- as.numeric(stats::quantile(cj, nd_frac[j]))
    lcmrl <- as.numeric(stats::quantile(cj, nd_frac[j] + dnq_frac[j]))
    left <- ifelse(cj < mdl, 0, ifelse(cj < lcmrl, mdl, cj))
    right <- ifelse(cj < mdl, mdl, ifelse(cj < lcmrl, lcmrl, cj))
    as_interval_data(left, right, log_transform = TRUE)
  })
}

# ---- competing fills -------------------------------------------------------

#' Naive substitution of the detection limit (or a fraction of it).
naive_fill <- function(bnds, at = c("limit", "half", "sqrt2")) {
  at <- match.arg(at)
  filled <- bnds$data_wide
  m <- bnds$to_impute
  hi <- bnds$hi_mat[m]                          # log(limit)
  filled[m] <- switch(at,
    limit = hi,
    half  = hi - log(2),                        # limit / 2
    sqrt2 = hi - 0.5 * log(2)                   # limit / sqrt(2)
  )
  filled
}

# ---- metrics ---------------------------------------------------------------

#' Point-recovery metrics of a completed matrix against the log-scale truth.
metrics_point <- function(filled, z_true) {
  p <- ncol(z_true)
  probs <- c(0.1, 0.5, 0.9)
  q_bias <- vapply(seq_len(p), function(j) {
    as.numeric(stats::quantile(filled[, j], probs) -
                 stats::quantile(z_true[, j], probs))
  }, numeric(length(probs)))
  off <- upper.tri(diag(p))
  list(
    mean_bias = mean(colMeans(filled) - colMeans(z_true)),
    mean_absbias = mean(abs(colMeans(filled) - colMeans(z_true))),
    q50_absbias = mean(abs(q_bias[2, ])),
    q10_absbias = mean(abs(q_bias[1, ])),
    q90_absbias = mean(abs(q_bias[3, ])),
    cor_mae = if (p > 1) mean(abs((cor(filled) - cor(z_true))[off])) else NA_real_
  )
}

#' Multiple-imputation coverage of the true analyte means (Rubin's rules).
#'
#' Generates `M` gsimp imputations, pools the per-analyte means, and returns the
#' fraction of analytes whose 95% MI confidence interval covers the true mean.
#' A calibrated method returns ~0.95.
mi_coverage <- function(bnds, z_true, M, iters_all, imp_model, seed) {
  n <- nrow(z_true)
  p <- ncol(z_true)
  means <- matrix(NA_real_, M, p)
  vars <- matrix(NA_real_, M, p)
  for (mi in seq_len(M)) {
    set.seed(seed + mi)
    fm <- gsimp_impute(bnds, iters_all = iters_all, imp_model = imp_model)
    means[mi, ] <- colMeans(fm)
    vars[mi, ] <- apply(fm, 2, stats::var) / n   # sampling variance of the mean
  }
  qbar <- colMeans(means)
  ubar <- colMeans(vars)                          # within-imputation variance
  b <- apply(means, 2, stats::var)                # between-imputation variance
  Tvar <- ubar + (1 + 1 / M) * b                  # total variance
  df <- (M - 1) * (1 + ubar / ((1 + 1 / M) * pmax(b, .Machine$double.eps)))^2
  df[!is.finite(df)] <- M - 1
  tcrit <- stats::qt(0.975, df)
  half <- tcrit * sqrt(Tvar)
  true_mean <- colMeans(z_true)
  covered <- (true_mean >= qbar - half) & (true_mean <= qbar + half)
  list(coverage = mean(covered), ci_width = mean(2 * half))
}

# ---- scenario grid ---------------------------------------------------------

#' Build the condition grid. "quick" is for a fast sanity run; "full" is the
#' publication grid. Edit freely -- the runner is agnostic to grid size.
make_grid <- function(config = c("quick", "full")) {
  config <- match.arg(config)
  if (config == "quick") {
    expand.grid(
      n = 200L, p = 5L,
      rho = c(0, 0.6),
      nd_frac = c(0.15, 0.45),
      dnq_frac = 0.10,
      skew = c(0, 0.5),
      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
  } else {
    expand.grid(
      n = c(100L, 300L),
      p = c(5L, 15L),
      rho = c(0, 0.3, 0.6),
      nd_frac = c(0.15, 0.35, 0.55),
      dnq_frac = 0.10,
      skew = c(0, 0.75),
      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
  }
}

# ---- runner ----------------------------------------------------------------

#' Run the study: for every scenario x replication, simulate, censor, impute,
#' and record metrics for gsimp vs naive plus MI coverage. Returns a tidy
#' data frame (one row per scenario x replication).
run_validation <- function(grid, n_rep = 100L, M = 20L, iters_all = 30L,
                            imp_model = "ridge", base_seed = 20240101L,
                            verbose = TRUE) {
  out <- vector("list", nrow(grid) * n_rep)
  k <- 0L
  for (s in seq_len(nrow(grid))) {
    cond <- grid[s, ]
    for (r in seq_len(n_rep)) {
      seed <- base_seed + s * 100000L + r
      set.seed(seed)
      z <- simulate_truth(cond$n, cond$p, cond$rho, skew = cond$skew)
      bnds <- build_bounds(censor_three_tier(z, cond$nd_frac, cond$dnq_frac))

      set.seed(seed + 1L)
      f_gs <- gsimp_impute(bnds, iters_all = iters_all, imp_model = imp_model)
      f_nv <- naive_fill(bnds, "limit")

      m_gs <- metrics_point(f_gs, z)
      m_nv <- metrics_point(f_nv, z)
      cov <- mi_coverage(bnds, z, M = M, iters_all = iters_all,
                         imp_model = imp_model, seed = seed + 1000L)

      k <- k + 1L
      out[[k]] <- data.frame(
        scenario = s, rep = r,
        n = cond$n, p = cond$p, rho = cond$rho,
        nd_frac = cond$nd_frac, dnq_frac = cond$dnq_frac, skew = cond$skew,
        censored_frac = mean(bnds$to_impute),
        gs_mean_bias = m_gs$mean_bias,
        gs_mean_absbias = m_gs$mean_absbias,
        gs_q50_absbias = m_gs$q50_absbias,
        gs_cor_mae = m_gs$cor_mae,
        naive_mean_bias = m_nv$mean_bias,
        naive_mean_absbias = m_nv$mean_absbias,
        naive_cor_mae = m_nv$cor_mae,
        mi_coverage = cov$coverage,
        mi_ci_width = cov$ci_width
      )
    }
    if (verbose) {
      message(sprintf("scenario %d/%d done (%d reps)", s, nrow(grid), n_rep))
    }
  }
  do.call(rbind, out)
}

#' Aggregate replications to per-scenario summaries, with Monte Carlo SEs and
#' RMSE of the mean estimate.
summarise_validation <- function(res) {
  by <- res[, c("scenario", "n", "p", "rho", "nd_frac", "skew")]
  key <- interaction(by, drop = TRUE)
  split_res <- split(res, key)
  do.call(rbind, lapply(split_res, function(d) {
    mcse <- function(x) stats::sd(x) / sqrt(length(x))
    data.frame(
      d[1, c("scenario", "n", "p", "rho", "nd_frac", "skew")],
      n_rep = nrow(d),
      censored_frac = mean(d$censored_frac),
      gs_bias = mean(d$gs_mean_bias),
      gs_bias_mcse = mcse(d$gs_mean_bias),
      gs_rmse = sqrt(mean(d$gs_mean_bias^2)),
      naive_bias = mean(d$naive_mean_bias),
      naive_rmse = sqrt(mean(d$naive_mean_bias^2)),
      gs_cor_mae = mean(d$gs_cor_mae),
      naive_cor_mae = mean(d$naive_cor_mae),
      mi_coverage = mean(d$mi_coverage),
      mi_coverage_mcse = mcse(d$mi_coverage),
      mi_ci_width = mean(d$mi_ci_width),
      row.names = NULL
    )
  }))
}

# ---- main ------------------------------------------------------------------

if (sys.nframe() == 0) {
  config <- Sys.getenv("CONFIG", "quick")
  n_rep <- as.integer(Sys.getenv("N_REP", "50"))
  M <- as.integer(Sys.getenv("M", "20"))
  imp_model <- Sys.getenv("IMP_MODEL", "ridge")

  message(sprintf("Running validation: config=%s, n_rep=%d, M=%d, imp_model=%s",
                  config, n_rep, M, imp_model))

  grid <- make_grid(config)
  res <- run_validation(grid, n_rep = n_rep, M = M, imp_model = imp_model)
  summ <- summarise_validation(res)

  dir.create("validation/results", showWarnings = FALSE, recursive = TRUE)
  stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
  saveRDS(list(config = config, n_rep = n_rep, M = M, imp_model = imp_model,
               grid = grid, results = res, summary = summ),
          sprintf("validation/results/validation_%s.rds", stamp))
  saveRDS(list(config = config, n_rep = n_rep, M = M, imp_model = imp_model,
               grid = grid, results = res, summary = summ),
          "validation/results/latest.rds")
  utils::write.csv(res, "validation/results/validation_results.csv",
                   row.names = FALSE)
  utils::write.csv(summ, "validation/results/validation_summary.csv",
                   row.names = FALSE)

  message("\n== Summary (bias, RMSE, MI coverage by scenario) ==")
  print(summ[, c("n", "p", "rho", "nd_frac", "skew", "censored_frac",
                 "gs_bias", "gs_rmse", "naive_bias", "naive_rmse",
                 "mi_coverage")], digits = 3)
}
