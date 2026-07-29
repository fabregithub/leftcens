# Shared internal helpers used across leftcens modules.

#' Recover original-scale interval bounds from a `cens_data` object
#'
#' Back-transforms the stored (possibly log-scale) bounds to the original
#' measurement scale, so non-detects become clean `[0, MDL]` intervals rather
#' than `(-Inf, log(MDL)]`. Returns full-length vectors plus a `keep` mask
#' flagging observations with both bounds present, so callers can subset the
#' bounds and any aligned covariates consistently.
#'
#' @param x A `cens_data` object.
#' @return A list with numeric `lo` and `hi` (full length, original scale), a
#'   logical `keep` mask, and the integer count of `dropped` rows.
#' @keywords internal
#' @noRd
cens_original_bounds <- function(x) {
  lo <- x$left
  hi <- x$right
  if (isTRUE(x$log_transform)) {
    lo <- exp(lo)   # exp(-Inf) == 0, the non-detect lower bound
    hi <- exp(hi)
  }
  keep <- !(is.na(lo) | is.na(hi))
  list(lo = lo, hi = hi, keep = keep, dropped = sum(!keep))
}
