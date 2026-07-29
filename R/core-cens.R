#' Interval-censored data categories
#'
#' The three-tier interval-censoring structure used throughout `leftcens`. The
#' category of each observation is *derived* from its `(left, right)` bounds, not
#' supplied by the user, so the two representations cannot drift out of sync.
#'
#' \describe{
#'   \item{`non_detect`}{`(0, MDL)` --- below the method detection limit.}
#'   \item{`detected_not_quantified`}{`(MDL, LCMRL)` --- detected but below the
#'     lowest concentration minimum reporting level.}
#'   \item{`quantified`}{`(value, value)` --- a point measurement.}
#' }
#'
#' @keywords internal
#' @noRd
cens_categories <- c("non_detect", "detected_not_quantified", "quantified")

# ---- low-level constructor -------------------------------------------------

#' Low-level `cens_data` constructor
#'
#' Assembles a `cens_data` object from already-validated, already-transformed
#' components. Performs no checking --- use [as_interval_data()] to build one from
#' user input.
#'
#' @param left,right Numeric vectors of interval bounds on the *working* scale
#'   (log scale when `log_transform` is `TRUE`).
#' @param category Character vector of censoring categories (see
#'   `cens_categories`).
#' @param log_transform Single logical: whether `left`/`right` are on the log
#'   scale.
#'
#' @return A `cens_data` object.
#' @keywords internal
#' @noRd
new_cens_data <- function(left, right, category, log_transform) {
  stopifnot(
    is.numeric(left), is.numeric(right),
    length(left) == length(right),
    is.logical(log_transform), length(log_transform) == 1L
  )
  structure(
    list(
      left = left,
      right = right,
      category = factor(category, levels = cens_categories),
      log_transform = log_transform
    ),
    class = "cens_data"
  )
}

# ---- validator -------------------------------------------------------------

#' Validate a `cens_data` object
#'
#' Checks the structural invariants of a `cens_data` object and returns it
#' invisibly, or throws an informative error.
#'
#' @param x A `cens_data` object.
#' @return `x`, invisibly, if valid.
#' @keywords internal
#' @noRd
validate_cens_data <- function(x) {
  if (!inherits(x, "cens_data")) {
    stop("`x` must be a <cens_data> object.", call. = FALSE)
  }
  left <- x$left
  right <- x$right

  n <- length(left)
  if (length(right) != n || length(x$category) != n) {
    stop("`left`, `right`, and `category` must have the same length.",
         call. = FALSE)
  }

  # left <= right elementwise, ignoring NA pairs.
  bad <- which(left > right)
  if (length(bad) > 0) {
    stop(
      "`left` must be <= `right` for every observation. Offending index(es): ",
      paste(utils::head(bad, 5), collapse = ", "),
      if (length(bad) > 5) ", ..." else "",
      ".",
      call. = FALSE
    )
  }

  if (any(!x$category %in% cens_categories & !is.na(x$category))) {
    stop("`category` must be one of: ",
         paste(cens_categories, collapse = ", "), ".", call. = FALSE)
  }

  invisible(x)
}

# ---- classification --------------------------------------------------------

#' Classify interval bounds into censoring categories
#'
#' @param left,right Numeric vectors of bounds on the *original* (untransformed)
#'   scale.
#' @return Character vector of censoring categories (`cens_categories`), or `NA`
#'   where a bound is `NA`.
#' @keywords internal
#' @noRd
classify_interval <- function(left, right) {
  ifelse(
    left == right, "quantified",
    ifelse(left <= 0, "non_detect", "detected_not_quantified")
  )
}

# ---- user-facing constructors ----------------------------------------------

#' Create an interval-censored data object
#'
#' Constructs a `cens_data` object from per-observation lower and upper interval
#' bounds. This is the single shared data representation consumed by every
#' descriptive-statistics and imputation function in `leftcens`, so the interval
#' semantics are defined exactly once.
#'
#' Each observation is classified (see `cens_categories`) from its bounds:
#' `(0, MDL)` is a non-detect, `(MDL, LCMRL)` is detected-but-not-quantified, and
#' `(value, value)` is a quantified point measurement. The classification is
#' derived, never supplied, so it cannot disagree with the bounds.
#'
#' Classic single-limit (Helsel-style) left censoring is a degenerate case of
#' this same representation --- see [as_interval_data_from_single()].
#'
#' @param left,right Numeric vectors of interval lower and upper bounds on the
#'   original measurement scale. Recycled to a common length. `left = 0` denotes
#'   a non-detect lower bound; under `log_transform` it maps to `-Inf`, which is
#'   expected, not an error.
#' @param log_transform Single logical. When `TRUE` (default), bounds are stored
#'   on the natural-log scale, matching how the imputation and NPMLE machinery
#'   operate. `left`/`right` must then be non-negative.
#'
#' @return A `cens_data` object: a list with elements `left`, `right`,
#'   `category` (a factor), and `log_transform`, with class `"cens_data"`.
#'
#' @examples
#' # non-detect, detected-not-quantified, quantified
#' as_interval_data(left = c(0, 1, 5), right = c(1, 3, 5))
#' @export
as_interval_data <- function(left, right, log_transform = TRUE) {
  if (!is.numeric(left) || !is.numeric(right)) {
    stop("`left` and `right` must be numeric vectors.", call. = FALSE)
  }
  if (length(log_transform) != 1L || !is.logical(log_transform) ||
      is.na(log_transform)) {
    stop("`log_transform` must be a single `TRUE` or `FALSE`.", call. = FALSE)
  }

  # Recycle to a common length (length-1 recycling only, like arithmetic).
  n <- max(length(left), length(right))
  if (n == 0L) {
    return(new_cens_data(double(), double(), character(), log_transform))
  }
  if (!length(left) %in% c(1L, n) || !length(right) %in% c(1L, n)) {
    stop("`left` and `right` must have the same length, ",
         "or one of them must have length 1.", call. = FALSE)
  }
  left <- rep_len(left, n)
  right <- rep_len(right, n)

  # Validate ordering on the original scale for clear error messages.
  bad <- which(left > right)
  if (length(bad) > 0) {
    stop(
      "`left` must be <= `right`. Offending index(es): ",
      paste(utils::head(bad, 5), collapse = ", "),
      if (length(bad) > 5) ", ..." else "", ".",
      call. = FALSE
    )
  }

  if (log_transform) {
    neg <- which(left < 0 | right < 0)
    if (length(neg) > 0) {
      stop(
        "With `log_transform = TRUE`, bounds must be non-negative. ",
        "Offending index(es): ",
        paste(utils::head(neg, 5), collapse = ", "),
        if (length(neg) > 5) ", ..." else "", ".",
        call. = FALSE
      )
    }
  }

  category <- classify_interval(left, right)

  if (log_transform) {
    left <- log(left)   # log(0) == -Inf, the expected non-detect lower bound
    right <- log(right)
  }

  validate_cens_data(
    new_cens_data(left, right, category, log_transform)
  )
}

#' Create interval-censored data from single-limit (Helsel-style) input
#'
#' A convenience wrapper around [as_interval_data()] for datasets in the classic
#' single-detection-limit format: one value per observation plus a flag for
#' whether it was censored. This is the two-tier collapse of the interval model
#' (`non_detect` and `quantified` only), so it flows through the same
#' `cens_data` object as fully interval-censored data.
#'
#' @param value Numeric vector of measured values (for detected observations) or
#'   detection limits (for censored observations).
#' @param censored Logical vector: `TRUE` where `value` is a detection limit
#'   (non-detect), `FALSE` where it is a quantified measurement. Recycled against
#'   `value`.
#' @param log_transform Passed to [as_interval_data()].
#'
#' @return A `cens_data` object.
#'
#' @examples
#' # 2.0 is a non-detect at DL = 2.0; 5.1 is a quantified measurement
#' as_interval_data_from_single(value = c(2.0, 5.1), censored = c(TRUE, FALSE))
#' @export
as_interval_data_from_single <- function(value, censored, log_transform = TRUE) {
  if (!is.numeric(value)) {
    stop("`value` must be a numeric vector.", call. = FALSE)
  }
  if (!is.logical(censored)) {
    stop("`censored` must be a logical vector.", call. = FALSE)
  }
  n <- max(length(value), length(censored))
  if (n > 0 && (!length(value) %in% c(1L, n) ||
                !length(censored) %in% c(1L, n))) {
    stop("`value` and `censored` must have the same length, ",
         "or one of them must have length 1.", call. = FALSE)
  }
  value <- rep_len(value, n)
  censored <- rep_len(censored, n)

  # censored -> (0, value); detected -> (value, value)
  left <- ifelse(censored, 0, value)
  right <- value

  as_interval_data(left, right, log_transform = log_transform)
}

# ---- methods ---------------------------------------------------------------

#' @export
print.cens_data <- function(x, ...) {
  n <- length(x$left)
  cat(sprintf(
    "<cens_data> %d observation%s (log_transform = %s)\n",
    n, if (n == 1L) "" else "s", x$log_transform
  ))

  if (n > 0) {
    counts <- table(factor(x$category, levels = cens_categories))
    for (nm in cens_categories) {
      cat(sprintf("  %-24s %d\n", nm, counts[[nm]]))
    }

    # Show bounds back-transformed to the original scale for readability.
    show <- seq_len(min(n, 6L))
    lo <- x$left[show]
    hi <- x$right[show]
    if (isTRUE(x$log_transform)) {
      lo <- exp(lo)
      hi <- exp(hi)
    }
    cat("\n")
    preview <- data.frame(
      left = lo, right = hi,
      category = as.character(x$category[show])
    )
    print(preview, row.names = FALSE)
    if (n > length(show)) cat(sprintf("  ... and %d more\n", n - length(show)))
  }
  invisible(x)
}

#' @export
format.cens_data <- function(x, ...) {
  sprintf("<cens_data[%d]>", length(x$left))
}

#' Test whether an object is a `cens_data`
#'
#' @param x Any object.
#' @return `TRUE` if `x` is a `cens_data` object, otherwise `FALSE`.
#' @export
is_cens_data <- function(x) {
  inherits(x, "cens_data")
}
