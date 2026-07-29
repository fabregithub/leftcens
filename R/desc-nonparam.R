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

  structure(
    list(
      fit = fit,
      method = method,
      backend = backend,
      log_transform = x$log_transform,
      n = length(lo),
      dropped = b$dropped,
      category = table(factor(x$category, levels = cens_categories))
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
#' @param x A `cens_np_fit` object.
#' @param probs Numeric vector of probabilities in `[0, 1]`.
#' @param ... Unused, for S3 consistency.
#'
#' @return A named numeric vector of quantile estimates.
#' @examples
#' fit <- desc_np(as_interval_data(c(0, 1, 5, 8), c(1, 3, 5, 8)))
#' quantile(fit, probs = c(0.5, 0.9))
#' @export
quantile.cens_np_fit <- function(x, probs = c(.1, .25, .5, .75, .9), ...) {
  q <- switch(
    x$backend,
    survival = as.numeric(quantile(x$fit, probs = probs, conf.int = FALSE)),
    icenReg  = as.numeric(icenReg::getFitEsts(x$fit, p = probs))
  )
  stats::setNames(q, paste0(format(100 * probs, trim = TRUE), "%"))
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
  cat("\nQuantiles:\n")
  print(quantile(x))
  invisible(x)
}
