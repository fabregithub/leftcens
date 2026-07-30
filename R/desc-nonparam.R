#' @importFrom stats quantile
NULL

# ---- internal helpers ------------------------------------------------------

#' Extract a survival step curve from a fitted NPMLE, on the original scale
#'
#' @param object A `cens_np_fit` object.
#' @return A data frame with `conc` (concentration) and `surv` (`S(conc)`).
#' @keywords internal
#' @noRd
np_surv_curve <- function(object) {
  if (object$backend == "survival") {
    data.frame(conc = object$fit$time, surv = object$fit$surv)
  } else {
    sc <- icenReg::getSCurves(object$fit)
    data.frame(
      conc = as.numeric(sc$Tbull_ints[, "upper"]),
      surv = as.numeric(sc$S_curves$baseline)
    )
  }
}

# ---- fitting ---------------------------------------------------------------

#' Non-parametric NPMLE summary of interval-censored data
#'
#' Fits a non-parametric maximum-likelihood estimate (NPMLE) of the distribution
#' underlying a `cens_data` object, using the Turnbull estimator. This is the
#' distribution-free descriptive summary for the three-tier interval-censoring
#' structure: it needs no parametric assumption about the concentration
#' distribution and correctly accounts for non-detects and
#' detected-but-not-quantified observations.
#'
#' The estimate is computed on the **original measurement scale**. Because the
#' Turnbull NPMLE depends only on the ordering and overlap of the intervals, it
#' is invariant to the monotone log transform stored in the `cens_data` object;
#' fitting on the original scale simply means the resulting quantiles are
#' reported directly as concentrations.
#'
#' @param x A `cens_data` object (see [as_interval_data()]).
#' @param method NPMLE algorithm. `"turnbull"` (default) uses
#'   [survival::survfit()] with `type = "interval2"`; `"wang"` uses
#'   [icenReg::ic_np()] (the EMICM algorithm).
#'
#' @return A `cens_np_fit` object: a list with the fitted `fit`, the `method`
#'   and `backend` used, `log_transform`, the number of observations `n`, the
#'   number of rows `dropped` for missing bounds, and the `category` counts.
#'
#' @seealso [quantile.cens_np_fit()] for quantile extraction and
#'   [plot.cens_np_fit()] for the cumulative-distribution plot.
#'
#' @examples
#' x <- as_interval_data(left = c(0, 1, 5, 8), right = c(1, 3, 5, 8))
#' fit <- desc_np(x)
#' quantile(fit)
#' @export
desc_np <- function(x, method = c("turnbull", "wang")) {
  if (!is_cens_data(x)) {
    stop("`x` must be a <cens_data> object; see `as_interval_data()`.",
         call. = FALSE)
  }
  method <- match.arg(method)

  b <- cens_original_bounds(x)
  lo <- b$lo[b$keep]
  hi <- b$hi[b$keep]
  if (length(lo) == 0) {
    stop("No non-missing observations to fit.", call. = FALSE)
  }

  backend <- if (method == "turnbull") "survival" else "icenReg"
  fit <- switch(
    method,
    turnbull = survival::survfit(
      survival::Surv(lo, hi, type = "interval2") ~ 1
    ),
    wang = icenReg::ic_np(cbind(lo, hi))
  )

  # Reporting limits on the original scale, derived from the censored cells: the
  # detection limit is the highest non-detect upper bound (MDL); the
  # quantitation limit is the highest upper bound among all non-quantified cells
  # (the LCMRL where a detected-not-quantified tier is present, else the MDL).
  cat_keep <- x$category[b$keep]
  nd_cell <- !is.na(cat_keep) & cat_keep == "non_detect"
  nq_cell <- !is.na(cat_keep) & cat_keep != "quantified"
  detection_limit <- if (any(nd_cell)) max(hi[nd_cell]) else NA_real_
  quantitation_limit <- if (any(nq_cell)) max(hi[nq_cell]) else NA_real_

  structure(
    list(
      fit = fit,
      method = method,
      backend = backend,
      log_transform = x$log_transform,
      n = length(lo),
      dropped = b$dropped,
      category = table(factor(x$category, levels = cens_categories)),
      detection_limit = detection_limit,
      quantitation_limit = quantitation_limit
    ),
    class = "cens_np_fit"
  )
}

# ---- methods ---------------------------------------------------------------

#' Quantiles of a non-parametric NPMLE fit
#'
#' Extracts quantiles (percentiles) of the estimated concentration distribution
#' from a [desc_np()] fit, on the original measurement scale.
#'
#' A quantile that falls below the quantitation limit is not reliably estimable
#' --- it lies in the censored region where the fit is extrapolating --- so by
#' default such quantiles are returned as `NA` rather than a misleading point
#' value. Set `ql = 0` (or `NA`) to disable this and return the raw NPMLE
#' estimates, or pass a different threshold (e.g. `ql = fit$detection_limit`).
#'
#' @param x A `cens_np_fit` object.
#' @param probs Numeric vector of probabilities in `[0, 1]`.
#' @param ql Reporting threshold on the original scale; quantile estimates below
#'   `ql` are set to `NA`. Defaults to the quantitation limit stored in the fit.
#'   Use `0` or `NA` to keep all estimates.
#' @param ... Unused, for S3 consistency.
#'
#' @return A named numeric vector of quantile estimates, with entries below `ql`
#'   set to `NA`.
#' @examples
#' fit <- desc_np(as_interval_data(c(0, 1, 5, 8), c(1, 3, 5, 8)))
#' quantile(fit, probs = c(0.5, 0.9))
#' quantile(fit, probs = c(0.5, 0.9), ql = 0)   # raw estimates, no masking
#' @export
quantile.cens_np_fit <- function(x, probs = c(.1, .25, .5, .75, .9),
                                 ql = x$quantitation_limit, ...) {
  q <- switch(
    x$backend,
    survival = as.numeric(quantile(x$fit, probs = probs, conf.int = FALSE)),
    icenReg  = as.numeric(icenReg::getFitEsts(x$fit, p = probs))
  )
  # Estimates below the quantitation limit are extrapolations into the censored
  # region; report them as NA (not reliably quantified).
  if (length(ql) == 1L && !is.na(ql) && ql > 0) {
    q[!is.na(q) & q < ql] <- NA_real_
  }
  stats::setNames(q, paste0(format(100 * probs, trim = TRUE), "%"))
}

#' Coerce a non-parametric NPMLE fit to a data frame
#'
#' Returns the quantile summary as a tidy data frame --- one row per requested
#' probability, with columns `probability`, `quantile` (the `"10%"`-style label),
#' and `estimate` --- ready to hand to a table-formatting package such as `gt`,
#' `flextable`, `knitr::kable()`, or `DT`. Estimates below the quantitation limit
#' are `NA` (see [quantile.cens_np_fit()]); the detection and quantitation limits
#' are attached as attributes.
#'
#' @param x A `cens_np_fit` object.
#' @param row.names,optional Passed along for S3 consistency; `row.names` sets
#'   the row names if supplied.
#' @param probs Probabilities to tabulate.
#' @param ql Quantitation-limit threshold (see [quantile.cens_np_fit()]).
#' @param ... Unused.
#'
#' @return A data frame with columns `probability`, `quantile`, `estimate`.
#' @examples
#' fit <- desc_np(as_interval_data(c(0, 1, 5, 8), c(1, 3, 5, 8)))
#' as.data.frame(fit)
#' # e.g. knitr::kable(as.data.frame(fit))  or  gt::gt(as.data.frame(fit))
#' @export
as.data.frame.cens_np_fit <- function(x, row.names = NULL, optional = FALSE,
                                      ..., probs = c(.1, .25, .5, .75, .9),
                                      ql = x$quantitation_limit) {
  q <- quantile(x, probs = probs, ql = ql)
  df <- data.frame(
    probability = probs,
    quantile = names(q),
    estimate = unname(as.numeric(q)),
    stringsAsFactors = FALSE
  )
  if (!is.null(row.names)) rownames(df) <- row.names
  attr(df, "detection_limit") <- x$detection_limit
  attr(df, "quantitation_limit") <- x$quantitation_limit
  df
}

#' Plot the cumulative distribution of a non-parametric NPMLE fit
#'
#' Draws the estimated cumulative distribution function `P(X <= conc)` as a step
#' function on the original measurement scale.
#'
#' @param x A `cens_np_fit` object.
#' @param xlab,ylab Axis labels.
#' @param ... Further arguments passed to [plot()].
#'
#' @return `x`, invisibly.
#' @export
plot.cens_np_fit <- function(x, xlab = "Concentration",
                             ylab = "Cumulative probability", ...) {
  cv <- np_surv_curve(x)
  ord <- order(cv$conc)
  conc <- cv$conc[ord]
  cdf <- 1 - cv$surv[ord]
  plot(
    c(conc[1], conc), c(0, cdf),
    type = "s", ylim = c(0, 1),
    xlab = xlab, ylab = ylab, ...
  )
  invisible(x)
}

#' @export
print.cens_np_fit <- function(x, ...) {
  cat(sprintf(
    "<cens_np_fit> NPMLE via %s (method = \"%s\")\n",
    x$backend, x$method
  ))
  cat(sprintf("  %d observation%s", x$n, if (x$n == 1L) "" else "s"))
  if (x$dropped > 0) cat(sprintf(" (%d dropped for missing bounds)", x$dropped))
  cat("\n")
  for (nm in cens_categories) {
    cat(sprintf("  %-24s %d\n", nm, x$category[[nm]]))
  }
  if (!is.na(x$quantitation_limit)) {
    cat(sprintf("  detection limit ~ %g | quantitation limit ~ %g\n",
                x$detection_limit, x$quantitation_limit))
  }
  cat("\nQuantiles (NA = below quantitation limit):\n")
  print(quantile(x))
  invisible(x)
}
