# Methods for the multi-analyte descriptive collections returned by desc_np()
# and desc_sp() when given a list of cens_data: `cens_np_fits` and
# `cens_sp_fits` are named lists of the corresponding single-analyte fits. The
# combining methods stack each analyte's tidy table with an `analyte` column.

#' @keywords internal
#' @noRd
bind_by_analyte <- function(fits, fun, ...) {
  parts <- lapply(names(fits), function(a) {
    d <- fun(fits[[a]], ...)
    if (nrow(d) == 0L) return(NULL)
    cbind(analyte = a, d, stringsAsFactors = FALSE)
  })
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (length(parts) == 0L) return(data.frame())
  do.call(rbind, parts)
}

# ---- non-parametric collection ---------------------------------------------

#' @export
print.cens_np_fits <- function(x, ...) {
  cat(sprintf("<cens_np_fits> non-parametric summaries for %d analyte%s\n",
              length(x), if (length(x) == 1L) "" else "s"))
  cat("\nQuantiles (NA = below quantitation limit):\n")
  print(quantile(x), ...)
  invisible(x)
}

#' Quantiles for a multi-analyte non-parametric summary
#'
#' @param x A `cens_np_fits` collection.
#' @param probs Probabilities to report.
#' @param ... Passed to [quantile.cens_np_fit()] (e.g. `ql`).
#' @return A matrix with one row per analyte and one column per probability.
#' @export
quantile.cens_np_fits <- function(x, probs = c(.1, .25, .5, .75, .9), ...) {
  m <- t(vapply(x, function(f) as.numeric(quantile(f, probs = probs, ...)),
                numeric(length(probs))))
  colnames(m) <- paste0(format(100 * probs, trim = TRUE), "%")
  rownames(m) <- names(x)
  m
}

#' Coerce a multi-analyte non-parametric summary to a data frame
#'
#' Stacks each analyte's quantile table (see [as.data.frame.cens_np_fit()]) with
#' a leading `analyte` column.
#'
#' @param x A `cens_np_fits` collection.
#' @param row.names,optional For S3 consistency.
#' @param ... Passed to [as.data.frame.cens_np_fit()] (e.g. `probs`, `ql`).
#' @return A tidy data frame.
#' @export
as.data.frame.cens_np_fits <- function(x, row.names = NULL, optional = FALSE,
                                       ...) {
  bind_by_analyte(x, as.data.frame, ...)
}

#' @export
tidy.cens_np_fits <- function(x, ...) {
  bind_by_analyte(x, generics::tidy, ...)
}

#' @export
glance.cens_np_fits <- function(x, ...) {
  bind_by_analyte(x, generics::glance, ...)
}

#' Plot the estimated distributions of several analytes
#'
#' Overlays each analyte's cumulative distribution on the original scale.
#'
#' @param x A `cens_np_fits` collection.
#' @param xlab,ylab Axis labels.
#' @param col Optional colours, one per analyte.
#' @param ... Further arguments passed to [plot()].
#' @return `x`, invisibly.
#' @export
plot.cens_np_fits <- function(x, xlab = "Concentration",
                              ylab = "Cumulative probability", col = NULL, ...) {
  curves <- lapply(x, np_surv_curve)
  k <- length(curves)
  if (is.null(col)) col <- seq_len(k)
  all_conc <- unlist(lapply(curves, `[[`, "conc"))
  plot(NA, xlim = range(all_conc, finite = TRUE), ylim = c(0, 1),
       xlab = xlab, ylab = ylab, ...)
  for (i in seq_len(k)) {
    cv <- curves[[i]]
    ord <- order(cv$conc)
    conc <- cv$conc[ord]
    graphics::lines(c(conc[1], conc), c(0, 1 - cv$surv[ord]),
                    type = "s", col = col[i])
  }
  graphics::legend("bottomright", legend = names(x), col = col, lty = 1,
                   bty = "n")
  invisible(x)
}

# ---- semi-parametric collection --------------------------------------------

#' @export
print.cens_sp_fits <- function(x, ...) {
  cat(sprintf("<cens_sp_fits> semi-parametric summaries for %d analyte%s\n",
              length(x), if (length(x) == 1L) "" else "s"))
  df <- as.data.frame(x)
  if (nrow(df) == 0L) {
    cat("  (intercept-only fits; no covariates)\n")
  } else {
    cat("\nCoefficients:\n")
    print(df, row.names = FALSE, ...)
  }
  invisible(x)
}

#' Coerce a multi-analyte semi-parametric summary to a data frame
#'
#' Stacks each analyte's coefficient table (see [as.data.frame.cens_sp_fit()])
#' with a leading `analyte` column.
#'
#' @param x A `cens_sp_fits` collection.
#' @param row.names,optional For S3 consistency.
#' @param ... Passed to [as.data.frame.cens_sp_fit()].
#' @return A tidy data frame, one row per analyte-by-term coefficient.
#' @export
as.data.frame.cens_sp_fits <- function(x, row.names = NULL, optional = FALSE,
                                       ...) {
  bind_by_analyte(x, as.data.frame, ...)
}

#' @export
tidy.cens_sp_fits <- function(x, ...) {
  bind_by_analyte(x, generics::tidy, ...)
}

#' @export
glance.cens_sp_fits <- function(x, ...) {
  bind_by_analyte(x, generics::glance, ...)
}
