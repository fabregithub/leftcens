# Build per-cell imputation bounds from cens_data objects.
#
# The GSimp-lineage Gibbs sampler imputes each censored/missing cell of a wide
# analyte matrix from a truncated normal constrained to that cell's (lo, hi)
# interval. This module turns the `cens_data` interval representation into the
# `lo_mat` / `hi_mat` / `data_wide` matrices that machinery consumes. It is
# self-contained original code: it does not depend on, and is not derived from,
# the (CC BY-NC-SA) GSimp source.

# ---- internal helpers ------------------------------------------------------

#' Normalise the `build_bounds()` input to a list of `cens_data` columns
#'
#' @param x A `cens_data` object or a list of them.
#' @return A named list of `cens_data` objects.
#' @keywords internal
#' @noRd
as_cens_data_list <- function(x) {
  if (is_cens_data(x)) {
    return(list(V1 = x))
  }
  if (is.list(x) && length(x) > 0 && all(vapply(x, is_cens_data, logical(1)))) {
    nms <- names(x)
    if (is.null(nms) || any(nms == "")) {
      names(x) <- paste0("V", seq_along(x))
    }
    return(x)
  }
  stop("`x` must be a <cens_data> object or a non-empty list of them.",
       call. = FALSE)
}

#' Bounds for a single `cens_data` column
#'
#' @param cd A `cens_data` object.
#' @return A list of four length-`n` vectors: `data`, `lo`, `hi`, `to_impute`.
#' @keywords internal
#' @noRd
column_bounds <- function(cd) {
  observed <- !is.na(cd$category) & cd$category == "quantified"

  # Quantified cells are observed (left == right); everything else is imputed.
  data <- ifelse(observed, cd$left, NA_real_)
  lo <- ifelse(observed, NA_real_,
               ifelse(is.na(cd$left), -Inf, cd$left))
  hi <- ifelse(observed, NA_real_,
               ifelse(is.na(cd$right), Inf, cd$right))

  list(data = data, lo = lo, hi = hi, to_impute = !observed)
}

# ---- constructor -----------------------------------------------------------

#' Build per-cell imputation bounds from interval-censored data
#'
#' Converts one or more `cens_data` columns into the per-cell bound matrices
#' consumed by `gsimp_impute()` (the imputation module). Each analyte column
#' becomes a column of three
#' aligned matrices: the observed values (`data_wide`, with censored cells
#' `NA`), and the lower/upper bounds (`lo_mat` / `hi_mat`, `NA` at observed
#' cells). Building all three together guarantees the imputation input and its
#' bounds share exactly the same missing-cell pattern.
#'
#' All matrices are on the stored scale of the input (log scale when the
#' `cens_data` objects were built with `log_transform = TRUE`, which is the
#' scale the Gibbs sampler operates on). On the log scale a non-detect `(0, MDL)`
#' becomes the left-censored bound `(-Inf, log MDL)` directly, with no special
#' casing. Cells whose bounds are both `NA` are treated as ordinary
#' (unconstrained) missing values and imputed from `(-Inf, Inf)`.
#'
#' @param x A single `cens_data` object (one analyte column) or a list of
#'   `cens_data` objects, one per analyte. When a list, all elements must have
#'   the same length and the same `log_transform` flag; list names become the
#'   output column names.
#'
#' @return A `cens_bounds` object: a list with numeric matrices `data_wide`,
#'   `lo_mat`, and `hi_mat` (all `n` observations x `p` analytes), a logical
#'   matrix `to_impute`, and the `log_transform` flag.
#'
#' @seealso [as_interval_data()] for building the `cens_data` inputs.
#'
#' @examples
#' a <- as_interval_data(left = c(0, 2, 5), right = c(1, 4, 5))   # analyte A
#' b <- as_interval_data(left = c(3, 0, 7), right = c(3, 1, 7))   # analyte B
#' bnds <- build_bounds(list(A = a, B = b))
#' bnds$lo_mat
#' @export
build_bounds <- function(x) {
  cols <- as_cens_data_list(x)

  lens <- vapply(cols, function(cd) length(cd$left), integer(1))
  if (length(unique(lens)) != 1L) {
    stop("All `cens_data` columns must have the same length; got: ",
         paste(lens, collapse = ", "), ".", call. = FALSE)
  }
  logs <- vapply(cols, function(cd) isTRUE(cd$log_transform), logical(1))
  if (length(unique(logs)) != 1L) {
    stop("All `cens_data` columns must share the same `log_transform` flag.",
         call. = FALSE)
  }

  n <- lens[[1]]
  parts <- lapply(cols, column_bounds)
  col_nms <- names(cols)

  to_mat <- function(field) {
    m <- vapply(parts, function(p) p[[field]], numeric(n))
    dim(m) <- c(n, length(parts))
    colnames(m) <- col_nms
    m
  }

  structure(
    list(
      data_wide = to_mat("data"),
      lo_mat = to_mat("lo"),
      hi_mat = to_mat("hi"),
      to_impute = to_mat("to_impute") == 1,
      log_transform = logs[[1]]
    ),
    class = "cens_bounds"
  )
}

# ---- methods ---------------------------------------------------------------

#' Test whether an object is a `cens_bounds`
#'
#' @param x Any object.
#' @return `TRUE` if `x` is a `cens_bounds` object, otherwise `FALSE`.
#' @export
is_cens_bounds <- function(x) {
  inherits(x, "cens_bounds")
}

#' @export
print.cens_bounds <- function(x, ...) {
  d <- dim(x$data_wide)
  n_imp <- sum(x$to_impute)
  cat(sprintf(
    "<cens_bounds> %d observation%s x %d analyte%s (log_transform = %s)\n",
    d[1], if (d[1] == 1L) "" else "s",
    d[2], if (d[2] == 1L) "" else "s",
    x$log_transform
  ))
  cat(sprintf("  cells to impute: %d of %d (%.1f%%)\n",
              n_imp, prod(d), 100 * n_imp / prod(d)))
  # How many of the imputable cells are left-censored (lo == -Inf) vs bounded.
  imp_lo <- x$lo_mat[x$to_impute]
  left_cens <- sum(is.infinite(imp_lo) & imp_lo < 0)
  cat(sprintf("    left-censored (lo = -Inf): %d\n", left_cens))
  cat(sprintf("    finite interval:           %d\n", n_imp - left_cens))
  invisible(x)
}
