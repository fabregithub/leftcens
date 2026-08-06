# Multiple imputation with pooled inference (Rubin's rules) for gsimp_impute().
#
# gsimp_impute() returns ONE completed dataset -- the multiple-imputation
# primitive. Valid inference needs several: generate m completions, pool the
# per-analyte means with Rubin's rules, and propagate the imputation uncertainty
# into the standard errors. gsimp_mi() wraps that loop, reports the fraction of
# missing information (FMI), can pick m adaptively (increment until the pooled
# inference stabilises), and returns the completed datasets for a downstream
# analysis if wanted. For imp_model = "copula" it uses the recommended hybrid:
# the plug-in margin for the point estimate, drawn margins for the between-
# imputation variance.

#' Pool per-analyte means across imputations (Rubin's rules) with FMI.
#'
#' `point_means`/`point_vars` are `m x p` (the per-imputation analyte means and
#' their within-imputation sampling variances); `draw_means` is `m x p`, the set
#' used for the between-imputation variance (identical to `point_means` unless the
#' copula hybrid supplies a separate drawn set).
#' @keywords internal
#' @noRd
pool_means <- function(point_means, point_vars, draw_means, conf = 0.95) {
  m <- nrow(draw_means)
  qbar <- colMeans(point_means)
  ubar <- colMeans(point_vars)
  if (m < 2L) {                                   # no between-variance from m=1
    p <- length(qbar)
    return(data.frame(estimate = qbar, within_var = ubar,
                      between_var = NA_real_, total_var = ubar,
                      df = Inf, ci_lo = NA_real_, ci_hi = NA_real_,
                      fmi = NA_real_))
  }
  B  <- apply(draw_means, 2, stats::var)
  Tv <- ubar + (1 + 1 / m) * B
  r  <- (1 + 1 / m) * B / pmax(ubar, .Machine$double.eps)  # rel. variance increase
  df <- (m - 1) * (1 + 1 / r)^2
  df[!is.finite(df)] <- m - 1
  lambda <- (r + 2 / (df + 3)) / (r + 1)          # fraction of missing information
  half <- stats::qt((1 + conf) / 2, df) * sqrt(Tv)
  data.frame(estimate = qbar, within_var = ubar, between_var = B, total_var = Tv,
             df = df, ci_lo = qbar - half, ci_hi = qbar + half, fmi = lambda)
}

#' Multiple imputation of censored data with pooled inference
#'
#' Generates several imputations with [gsimp_impute()] and pools the per-analyte
#' means by Rubin's rules, so the imputation uncertainty is carried into the
#' standard errors (a single imputation understates it). Reports the fraction of
#' missing information (FMI) per analyte, which governs how many imputations are
#' needed. With `adaptive = TRUE` it increases `m` in steps until the pooled
#' interval width stabilises --- useful because the required `m` falls with sample
#' size for a proper model (the imputation-model uncertainty shrinks as the fit
#' tightens), so a fixed `m` is wasteful for large `n` and risky for small `n`.
#'
#' For `imp_model = "copula"` the recommended hybrid is used automatically: the
#' point estimate comes from plug-in-margin imputations and the between-imputation
#' variance from drawn-margin imputations, so intervals stay calibrated without a
#' point bias.
#'
#' @param x A `cens_bounds` object from [build_bounds()].
#' @param m Number of imputations (the starting number when `adaptive = TRUE`).
#' @param imp_model,iters_all,initial Passed to [gsimp_impute()].
#' @param conf Confidence level for the pooled interval.
#' @param adaptive If `TRUE`, add imputations in batches of `m_step` until the
#'   mean interval width changes by less than `tol` (relative) or `m_max` is hit.
#' @param m_max,m_step,tol Adaptive-mode controls.
#' @param return_imputations If `TRUE`, keep the completed datasets in the result
#'   (the proper drawn-margin set for the copula) for a downstream analysis.
#' @param n_cores Fork workers for parallel imputation (serial on Windows or when
#'   `1`).
#' @param seed Optional integer for reproducibility.
#'
#' @return A `gsimp_mi` object: `pooled` (per-analyte estimate, variances, df,
#'   CI, and `fmi`), the `m` used, the adaptive `path`, and optionally
#'   `imputations`.
#' @seealso [gsimp_impute()] for a single imputation.
#' @examples
#' \donttest{
#' a <- as_interval_data(c(0, 2, 5, 8, 3), c(1, 4, 5, 8, 3))
#' b <- as_interval_data(c(3, 0, 7, 2, 6), c(3, 1, 7, 4, 6))
#' fit <- gsimp_mi(build_bounds(list(A = a, B = b)), m = 10)
#' fit$pooled
#' }
#' @export
gsimp_mi <- function(x, m = 20L, imp_model = "copula", iters_all = 10L,
                     initial = "bounds", conf = 0.95,
                     adaptive = FALSE, m_max = 100L, m_step = 10L, tol = 0.01,
                     return_imputations = FALSE, n_cores = 1L, seed = NULL) {
  if (!is_cens_bounds(x)) {
    stop("`x` must be a <cens_bounds> object; see `build_bounds()`.", call. = FALSE)
  }
  is_cop <- identical(imp_model, "copula")
  p <- ncol(x$data_wide)
  base <- if (is.null(seed)) sample.int(1e6L, 1L) else as.integer(seed)

  gen <- function(seeds, draw) {
    fn <- function(s) {
      set.seed(s)
      gsimp_impute(x, iters_all = iters_all, initial = initial,
                   imp_model = imp_model, margin_draw = draw)
    }
    if (n_cores > 1L) parallel::mclapply(seeds, fn, mc.cores = n_cores) else lapply(seeds, fn)
  }
  summ <- function(imps) list(
    means = t(vapply(imps, colMeans, numeric(p))),
    vars  = t(vapply(imps, function(f) apply(f, 2, stats::var) / nrow(f), numeric(p)))
  )

  Pm <- Pv <- Dm <- NULL                          # accumulated summaries
  imps_store <- list(); path <- list(); done <- 0L
  repeat {
    k <- if (adaptive && done > 0L) m_step else m
    seeds <- base + done + seq_len(k)
    pt <- gen(seeds, draw = FALSE)                # plug-in (copula) / the set (else)
    sp <- summ(pt); Pm <- rbind(Pm, sp$means); Pv <- rbind(Pv, sp$vars)
    if (is_cop) {
      dr <- gen(seeds + 500000L, draw = TRUE)     # drawn-margin -> variance
      Dm <- rbind(Dm, summ(dr)$means)
      if (return_imputations) imps_store <- c(imps_store, dr)
    } else {
      Dm <- Pm
      if (return_imputations) imps_store <- c(imps_store, pt)
    }
    done <- done + k
    pooled <- pool_means(Pm, Pv, Dm, conf)
    path[[length(path) + 1L]] <- data.frame(
      m = done, mean_ci_width = mean(pooled$ci_hi - pooled$ci_lo),
      mean_fmi = mean(pooled$fmi))
    if (!adaptive || done >= m_max) break
    if (length(path) >= 2L) {
      pw <- path[[length(path) - 1L]]$mean_ci_width
      cw <- path[[length(path)]]$mean_ci_width
      if (is.finite(pw) && is.finite(cw) &&
          abs(cw - pw) / max(pw, 1e-8) < tol) break
    }
  }

  nm <- colnames(x$data_wide)
  pooled <- cbind(analyte = if (is.null(nm)) seq_len(p) else nm, pooled)
  structure(list(
    pooled = pooled, m = done, imp_model = imp_model, adaptive = adaptive,
    path = do.call(rbind, path),
    imputations = if (return_imputations) imps_store else NULL,
    config = list(iters_all = iters_all, conf = conf, n_cores = n_cores, seed = base)
  ), class = "gsimp_mi")
}

#' @export
print.gsimp_mi <- function(x, ...) {
  cat(sprintf("<gsimp_mi> %d imputation%s, %s model%s\n",
              x$m, if (x$m == 1L) "" else "s", x$imp_model,
              if (x$adaptive) " (adaptive m)" else ""))
  cat(sprintf("  pooled per-analyte mean (%.0f%% CI, Rubin's rules); FMI = fraction of missing info\n",
              100 * x$config$conf))
  d <- x$pooled
  print(d[, c("analyte", "estimate", "ci_lo", "ci_hi", "fmi")],
        row.names = FALSE, digits = 3)
  fmax <- max(x$pooled$fmi, na.rm = TRUE)
  if (is.finite(fmax)) {
    cat(sprintf("  max FMI = %.2f; a rough guide is m >= 100*FMI ~ %d imputations.\n",
                fmax, ceiling(100 * fmax)))
  }
  invisible(x)
}
