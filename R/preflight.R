# Pre-flight reliability check: simulate a study design like the user's, impute,
# and report the bias / coverage they can expect -- so imputation can be vetted
# before it is trusted on real data. This is a self-contained Monte Carlo tool
# built on the package's own gsimp_impute(); it is the productised form of the
# validation study shipped in the source `validation/` directory.

# ---- internal simulation / evaluation helpers ------------------------------

#' Simulate correlated log-scale truth (exchangeable correlation, optional skew)
#'
#' `skew` is the sinh-arcsinh shape and may be a scalar or a length-`p` vector
#' (per-analyte skewness).
#' @keywords internal
#' @noRd
pf_simulate <- function(n, p, rho, mu, sd, skew) {
  R <- matrix(rho, p, p)
  diag(R) <- 1
  z <- matrix(stats::rnorm(n * p), n, p) %*% chol(R)
  skew <- rep_len(skew, p)
  if (any(skew != 0)) z <- sinh(sweep(asinh(z), 2, skew, "+"))  # per-analyte skew
  z <- sweep(z, 2, rep_len(sd, p), "*")
  sweep(z, 2, rep_len(mu, p), "+")
}

#' Censor a log-scale truth matrix into per-analyte three-tier `cens_data`
#' @keywords internal
#' @noRd
pf_censor <- function(z, nd_frac, dnq_frac) {
  conc <- exp(z)
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

#' Per-analyte multiple-imputation coverage of the true mean (Rubin's rules)
#' @keywords internal
#' @noRd
pf_mi_coverage <- function(bounds, z, M, iters_all, imp_model, seed) {
  n <- nrow(z)
  p <- ncol(z)
  means <- vars <- matrix(NA_real_, M, p)
  for (mi in seq_len(M)) {
    set.seed(seed + mi)
    fm <- gsimp_impute(bounds, iters_all = iters_all, imp_model = imp_model)
    means[mi, ] <- colMeans(fm)
    vars[mi, ] <- apply(fm, 2, stats::var) / n
  }
  qbar <- colMeans(means)
  ubar <- colMeans(vars)
  b <- apply(means, 2, stats::var)
  Tvar <- ubar + (1 + 1 / M) * b
  df <- (M - 1) * (1 + ubar / ((1 + 1 / M) * pmax(b, .Machine$double.eps)))^2
  df[!is.finite(df)] <- M - 1
  half <- stats::qt(0.975, df) * sqrt(Tvar)
  true_mean <- colMeans(z)
  (true_mean >= qbar - half) & (true_mean <= qbar + half)
}

# ---- main pre-flight -------------------------------------------------------

#' Pre-flight reliability check for a study design
#'
#' Simulates many datasets matching a described study design (sample size,
#' per-analyte non-detect fractions, correlation, ...), imputes each with
#' [gsimp_impute()], and reports the bias and multiple-imputation coverage you
#' can expect --- so imputation can be vetted *before* it is trusted on real
#' data. An analyte is flagged `reliable` when its mean bias is small, its 95%
#' interval covers the truth often enough, and its non-detect fraction is below
#' 50% (the point above which the median itself is censored).
#'
#' This is a Monte Carlo procedure and takes some time (it fits `n_rep * (M + 1)`
#' imputations). It assumes an approximately log-normal, exchangeably-correlated
#' data-generating process; use it as a guide, not a guarantee.
#'
#' @param n Number of samples (rows) in the design.
#' @param nd_frac Non-detect fraction, as a single value applied to all analytes
#'   or a length-`p` vector of per-analyte fractions.
#' @param p Number of analytes. Defaults to `length(nd_frac)` when that is a
#'   vector.
#' @param rho Exchangeable correlation between analytes on the log scale, in
#'   `[0, 1)`.
#' @param dnq_frac Detected-not-quantified fraction (the middle tier), scalar or
#'   length-`p`. Default `0` (pure left-censoring).
#' @param skew sinh-arcsinh skewness of the marginals (`0` = log-normal); scalar
#'   or length-`p` for per-analyte skewness. [preflight_from_data()] estimates
#'   this from your data. Right-skew is the main driver of mean-imputation bias, so
#'   a non-zero value here is what makes the reported reliability honest for
#'   skewed data.
#' @param mu,sd Log-scale mean and standard deviation of the analytes; scalar or
#'   length-`p`. Defaults span a modest range of means at unit sd.
#' @param n_rep Number of simulated datasets.
#' @param M Number of imputations per dataset for the coverage estimate.
#' @param iters_all Gibbs sweeps per imputation (see [gsimp_impute()]).
#' @param imp_model Conditional model (see [gsimp_impute()]); default `"copula"`.
#' @param bias_tol,coverage_tol Reliability thresholds: maximum |mean bias| and
#'   minimum coverage.
#' @param seed Optional integer for reproducibility.
#'
#' @return A `preflight` object: a list with `by_analyte` (a data frame of
#'   per-analyte `nd_frac`, `mean_bias`, `median_bias`, `rmse`, `coverage`, and
#'   `reliable`) and the `config` used.
#'
#' @seealso [preflight_from_data()] to derive the design from your own data.
#' @examples
#' \donttest{
#' # three analytes at 30% / 45% / 55% non-detects, n = 120
#' preflight_reliability(n = 120, nd_frac = c(0.30, 0.45, 0.55),
#'                       n_rep = 10, M = 5)
#' }
#' @export
preflight_reliability <- function(n, nd_frac, p = NULL, rho = 0.5,
                                  dnq_frac = 0, skew = 0, mu = NULL, sd = 1,
                                  n_rep = 30, M = 20, iters_all = 10,
                                  imp_model = "copula",
                                  bias_tol = 0.10, coverage_tol = 0.90,
                                  seed = NULL) {
  if (is.null(p)) p <- length(nd_frac)
  if (p < 1) stop("`p` must be at least 1.", call. = FALSE)
  nd_frac <- rep_len(nd_frac, p)
  dnq_frac <- rep_len(dnq_frac, p)
  if (any(nd_frac < 0 | nd_frac >= 1) || any(dnq_frac < 0) ||
      any(nd_frac + dnq_frac >= 1)) {
    stop("`nd_frac` and `nd_frac + dnq_frac` must lie in [0, 1).", call. = FALSE)
  }
  if (rho < 0 || rho >= 1) {
    stop("`rho` must be in [0, 1).", call. = FALSE)
  }
  if (any(c(n, n_rep, M, iters_all) < 1)) {
    stop("`n`, `n_rep`, `M`, and `iters_all` must be positive.", call. = FALSE)
  }
  if (is.null(mu)) mu <- seq(0, 1.5, length.out = p)
  mu <- rep_len(mu, p)
  sd <- rep_len(sd, p)
  skew <- rep_len(skew, p)

  base <- if (is.null(seed)) sample.int(1e6L, 1L) else as.integer(seed)
  bias <- medbias <- cov <- matrix(NA_real_, n_rep, p)

  for (r in seq_len(n_rep)) {
    set.seed(base + r)
    z <- pf_simulate(n, p, rho, mu, sd, skew)
    bounds <- build_bounds(pf_censor(z, nd_frac, dnq_frac))

    set.seed(base + r + 500000L)
    f <- gsimp_impute(bounds, iters_all = iters_all, imp_model = imp_model)
    bias[r, ] <- colMeans(f) - colMeans(z)
    medbias[r, ] <- vapply(seq_len(p),
                           function(j) stats::median(f[, j]) -
                             stats::median(z[, j]), numeric(1))
    cov[r, ] <- pf_mi_coverage(bounds, z, M, iters_all, imp_model,
                               seed = base + r + 1000000L)
  }

  mean_bias <- colMeans(bias)
  coverage <- colMeans(cov)
  by_analyte <- data.frame(
    analyte = seq_len(p),
    nd_frac = nd_frac,
    mean_bias = mean_bias,
    median_bias = colMeans(medbias),
    rmse = sqrt(colMeans(bias^2)),
    coverage = coverage,
    reliable = abs(mean_bias) <= bias_tol & coverage >= coverage_tol &
      nd_frac < 0.5
  )

  structure(
    list(
      by_analyte = by_analyte,
      config = list(n = n, p = p, rho = rho, dnq_frac = dnq_frac, skew = skew,
                    n_rep = n_rep, M = M, iters_all = iters_all,
                    imp_model = imp_model, bias_tol = bias_tol,
                    coverage_tol = coverage_tol, seed = base)
    ),
    class = "preflight"
  )
}

#' Pre-flight reliability check derived from your own data
#'
#' A convenience wrapper around [preflight_reliability()] that reads the study
#' design directly from your `cens_data`: it takes the per-analyte non-detect
#' fractions from the data, fits each analyte's log-scale margin (location,
#' spread, and **skewness**) by interval-censored MLE, and estimates the
#' correlation from a quick imputation --- all used only to parameterise the
#' simulation. It then runs the pre-flight for that design. Estimating the skew
#' matters, and from the censored likelihood rather than the imputed values: a
#' log-normal assumption would make the report overly optimistic for skewed
#' analytes, and a symmetric imputation would hide the skew. Pass an explicit
#' `skew` to override the estimate.
#'
#' Because the moments and correlation are estimated (in part from imputed
#' values), treat the result as an approximate guide to what the imputation can
#' deliver for data like yours.
#'
#' @param x A `cens_data` object or a list of them (one per analyte), as passed
#'   to [build_bounds()].
#' @param imp_model,iters_all Passed to [gsimp_impute()] for both the
#'   parameter-estimation impute and the pre-flight.
#' @param ... Further arguments passed to [preflight_reliability()] (e.g.
#'   `n_rep`, `M`, `seed`).
#'
#' @return A `preflight` object; see [preflight_reliability()].
#' @examples
#' \donttest{
#' a <- as_interval_data(c(0, 2, 5, 8, 3), c(1, 4, 5, 8, 3))
#' b <- as_interval_data(c(3, 0, 7, 2, 6), c(3, 1, 7, 4, 6))
#' preflight_from_data(list(a = a, b = b), n_rep = 10, M = 5)
#' }
#' @export
preflight_from_data <- function(x, imp_model = "copula", iters_all = 10, ...) {
  cols <- as_cens_data_list(x)
  bounds <- build_bounds(cols)
  filled <- gsimp_impute(bounds, iters_all = iters_all, imp_model = imp_model)

  n <- nrow(filled)
  p <- ncol(filled)
  nd_frac <- vapply(cols,
                    function(cd) mean(cd$category == "non_detect", na.rm = TRUE),
                    numeric(1))
  R <- suppressWarnings(stats::cor(filled))
  rho <- if (p > 1) mean(R[upper.tri(R)], na.rm = TRUE) else 0
  rho <- min(max(rho, 0), 0.95)
  if (!is.finite(rho)) rho <- 0

  # Per-analyte margin (location, scale, sinh-arcsinh skew) by interval-censored
  # MLE. Estimating skew from the censored likelihood is robust; the sample
  # skewness of the imputed data would be ATTENUATED by a symmetric (tobit) fill,
  # hiding the very skew that biases the mean. `mu`/`sd` also come from the fit.
  miss <- bounds$to_impute
  fits <- lapply(seq_len(p), function(j)
    fit_shash_margin(filled[!miss[, j], j], bounds$lo_mat[miss[, j], j],
                     bounds$hi_mat[miss[, j], j]))
  mu <- vapply(fits, function(f) f$mu, numeric(1))
  sd <- vapply(fits, function(f) f$sigma, numeric(1))
  skew_hat <- vapply(fits, function(f) f$eps, numeric(1))
  skew_hat[!is.finite(skew_hat)] <- 0

  dots <- list(...)
  if (is.null(dots$skew)) dots$skew <- skew_hat

  do.call(preflight_reliability,
          c(list(n = n, nd_frac = nd_frac, p = p, rho = rho, mu = mu, sd = sd,
                 imp_model = imp_model, iters_all = iters_all), dots))
}

#' @export
print.preflight <- function(x, ...) {
  cfg <- x$config
  cat(sprintf("<preflight> reliability for %d analyte%s, n = %d (%s model)\n",
              cfg$p, if (cfg$p == 1L) "" else "s", cfg$n, cfg$imp_model))
  cat(sprintf("  design: rho=%.2f, dnq=%.2f, skew=%.2f | %d reps, M=%d, iters=%d\n",
              cfg$rho, mean(cfg$dnq_frac), mean(cfg$skew), cfg$n_rep, cfg$M,
              cfg$iters_all))
  cat(sprintf("  reliable if |bias| <= %.2f AND coverage >= %.2f AND ND < 50%%\n\n",
              cfg$bias_tol, cfg$coverage_tol))
  d <- x$by_analyte
  print(d, row.names = FALSE, digits = 3)

  n_rel <- sum(d$reliable)
  cat(sprintf("\nVerdict: %d of %d analyte%s reliably imputable.\n",
              n_rel, nrow(d), if (nrow(d) == 1L) "" else "s"))
  hi_nd <- d$analyte[d$nd_frac >= 0.5]
  if (length(hi_nd) > 0) {
    cat(sprintf("  At/above 50%% ND (median not estimable): %s\n",
                paste(hi_nd, collapse = ", ")))
  }
  # Skew is the main cause of mean-imputation bias; point to the skew-robust model.
  if (mean(abs(cfg$skew)) >= 0.3 && !identical(cfg$imp_model, "copula")) {
    cat(sprintf(paste0("  Right-skew detected (sinh-arcsinh delta ~ %.2f): the mean",
                       " may be biased low.\n  Consider imp_model = \"copula\"",
                       " (skew-robust), or report the median.\n"),
                mean(cfg$skew)))
  }
  invisible(x)
}
