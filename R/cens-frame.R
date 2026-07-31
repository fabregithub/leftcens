# Build a named list of cens_data (one per analyte) from a wide data frame,
# so a whole lab export can flow into desc_np()/desc_sp() and build_bounds()
# in one call instead of a hand-written lapply.

#' Build a list of `cens_data` from a wide data frame
#'
#' Turns a wide laboratory-style table --- one row per sample, one or two columns
#' per analyte --- into a named list of [as_interval_data()] objects, one per
#' analyte. That list is exactly what [desc_np()], [desc_sp()], and
#' [build_bounds()] accept, so this is the single on-ramp from a raw data frame to
#' every part of the package.
#'
#' Two common layouts are supported:
#'
#' \describe{
#'   \item{`"left_right"`}{Each analyte has a lower- and upper-bound column,
#'     e.g. `A.left` / `A.right`. This carries the full three-tier structure
#'     (non-detect, detected-not-quantified, quantified).}
#'   \item{`"flagged"`}{Each analyte has a value column and a qualifier-flag
#'     column, e.g. `A` / `A.flag`, where a flag in `nd_flags` marks a non-detect
#'     (the value being the detection/reporting limit) and any other flag a
#'     quantified result. This two-tier layout is the most common lab export.}
#' }
#'
#' @param data A data frame (or object coercible to one).
#' @param analytes Optional character vector of analyte names. When `NULL`
#'   (default) they are inferred from the column names for the chosen `format`.
#' @param format `"left_right"` (columns `<analyte><left_suffix>` /
#'   `<analyte><right_suffix>`) or `"flagged"` (columns `<analyte>` /
#'   `<analyte><flag_suffix>`).
#' @param left_suffix,right_suffix Column-name suffixes for the `"left_right"`
#'   format.
#' @param flag_suffix Column-name suffix for the qualifier flag in the
#'   `"flagged"` format.
#' @param nd_flags Flag values that mark a non-detect in the `"flagged"` format.
#' @param log_transform Passed to [as_interval_data()].
#'
#' @return A named list of `cens_data` objects, with class `"cens_list"`.
#'
#' @seealso [as_interval_data()] for a single analyte; the "Preparing your own
#'   data" vignette.
#' @examples
#' # left/right layout
#' d <- data.frame(ID = 1:3, A.left = c(0, 2, 5), A.right = c(1, 4, 5),
#'                 B.left = c(3, 0, 7), B.right = c(3, 1, 7))
#' as_cens_list(d)
#'
#' # flagged layout
#' f <- data.frame(ID = 1:3, A = c(2.1, 0.5, 3.4), A.flag = c("", "<", ""))
#' as_cens_list(f, format = "flagged")
#' @export
as_cens_list <- function(data, analytes = NULL,
                         format = c("left_right", "flagged"),
                         left_suffix = ".left", right_suffix = ".right",
                         flag_suffix = ".flag", nd_flags = c("<", "ND", "U"),
                         log_transform = TRUE) {
  data <- as.data.frame(data)
  format <- match.arg(format)
  nms <- names(data)

  if (format == "left_right") {
    if (is.null(analytes)) {
      have_left <- nms[endsWith(nms, left_suffix)]
      analytes <- sub(paste0(gsub("([.\\\\])", "\\\\\\1", left_suffix), "$"),
                      "", have_left)
      analytes <- analytes[paste0(analytes, right_suffix) %in% nms]
    }
    if (length(analytes) == 0) {
      stop("No analytes found: expected columns like `<analyte>", left_suffix,
           "` and `<analyte>", right_suffix, "`.", call. = FALSE)
    }
    cl <- lapply(analytes, function(a) {
      lcol <- paste0(a, left_suffix)
      rcol <- paste0(a, right_suffix)
      miss <- setdiff(c(lcol, rcol), nms)
      if (length(miss)) {
        stop("Missing column(s) for analyte '", a, "': ",
             paste(miss, collapse = ", "), ".", call. = FALSE)
      }
      as_interval_data(data[[lcol]], data[[rcol]],
                       log_transform = log_transform)
    })
  } else {
    if (is.null(analytes)) {
      have_flag <- nms[endsWith(nms, flag_suffix)]
      analytes <- sub(paste0(gsub("([.\\\\])", "\\\\\\1", flag_suffix), "$"),
                      "", have_flag)
      analytes <- analytes[analytes %in% nms]
    }
    if (length(analytes) == 0) {
      stop("No analytes found: expected value columns with matching `<analyte>",
           flag_suffix, "` columns.", call. = FALSE)
    }
    cl <- lapply(analytes, function(a) {
      fcol <- paste0(a, flag_suffix)
      miss <- setdiff(c(a, fcol), nms)
      if (length(miss)) {
        stop("Missing column(s) for analyte '", a, "': ",
             paste(miss, collapse = ", "), ".", call. = FALSE)
      }
      value <- data[[a]]
      nd <- as.character(data[[fcol]]) %in% nd_flags
      as_interval_data(left = ifelse(nd, 0, value), right = value,
                       log_transform = log_transform)
    })
  }

  names(cl) <- analytes
  structure(cl, class = "cens_list")
}

#' @export
print.cens_list <- function(x, ...) {
  cat(sprintf("<cens_list> %d analyte%s\n", length(x),
              if (length(x) == 1L) "" else "s"))
  for (a in names(x)) {
    cd <- x[[a]]
    tab <- table(factor(cd$category, levels = cens_categories))
    cat(sprintf("  %-14s n = %-4d  ND %d | DNQ %d | Q %d\n",
                a, length(cd$left), tab[["non_detect"]],
                tab[["detected_not_quantified"]], tab[["quantified"]]))
  }
  invisible(x)
}

#' Test whether an object is a `cens_list`
#'
#' @param x Any object.
#' @return `TRUE` if `x` is a `cens_list`, otherwise `FALSE`.
#' @export
is_cens_list <- function(x) inherits(x, "cens_list")
