# broom-style tidy() / glance() methods for the descriptive fits. The generics
# come from the lightweight 'generics' package (also used by broom), re-exported
# here so `library(leftcens); tidy(fit)` works without loading broom. Methods
# return plain data frames (no tibble dependency).

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

# ---- non-parametric fit ----------------------------------------------------

#' Tidy a non-parametric NPMLE fit
#'
#' Returns the estimated quantiles as a one-row-per-probability data frame.
#'
#' @param x A `cens_np_fit` object.
#' @param probs Probabilities to report.
#' @param ql Quantitation-limit threshold; estimates below it are `NA`
#'   (see [quantile.cens_np_fit()]). Use `0` to keep the raw estimates.
#' @param ... Unused.
#' @return A data frame with `term`, `probability`, and `estimate`.
#' @examples
#' fit <- desc_np(as_interval_data(c(0, 1, 5, 8), c(1, 3, 5, 8)))
#' tidy(fit)
#' @export
tidy.cens_np_fit <- function(x, probs = c(.1, .25, .5, .75, .9),
                             ql = x$quantitation_limit, ...) {
  q <- quantile(x, probs = probs, ql = ql)
  data.frame(
    term = names(q),
    probability = probs,
    estimate = unname(as.numeric(q)),
    stringsAsFactors = FALSE
  )
}

#' Glance at a non-parametric NPMLE fit
#'
#' One-row model-level summary.
#'
#' @param x A `cens_np_fit` object.
#' @param ... Unused.
#' @return A one-row data frame.
#' @examples
#' glance(desc_np(as_interval_data(c(0, 1, 5, 8), c(1, 3, 5, 8))))
#' @export
glance.cens_np_fit <- function(x, ...) {
  data.frame(
    method = x$method,
    backend = x$backend,
    nobs = x$n,
    n_non_detect = as.integer(x$category[["non_detect"]]),
    n_detected_not_quantified =
      as.integer(x$category[["detected_not_quantified"]]),
    n_quantified = as.integer(x$category[["quantified"]]),
    detection_limit = x$detection_limit,
    quantitation_limit = x$quantitation_limit,
    stringsAsFactors = FALSE
  )
}

# ---- semi-parametric fit ---------------------------------------------------

#' Tidy a semi-parametric fit
#'
#' Returns the covariate coefficient table (one row per term).
#'
#' @param x A `cens_sp_fit` object.
#' @param exponentiate If `TRUE`, report `estimate` on the ratio scale (hazard
#'   ratio for `"ph"`, odds ratio for `"po"`).
#' @param ... Unused.
#' @return A data frame with `term`, `estimate`, and (when the fit was built with
#'   `bs_samples > 0`) `std.error`, `statistic`, and `p.value`.
#' @examples
#' x <- as_interval_data(c(0, 1, 5, 8, 0, 2), c(1, 3, 5, 8, 1, 4))
#' grp <- factor(c("a", "a", "a", "b", "b", "b"))
#' tidy(desc_sp(x, covariates = grp, model = "ph"))
#' @export
tidy.cens_sp_fit <- function(x, exponentiate = FALSE, ...) {
  cf <- stats::coef(x$fit)
  if (length(cf) == 0) {
    return(data.frame(term = character(0), estimate = numeric(0),
                      stringsAsFactors = FALSE))
  }
  df <- data.frame(term = names(cf), estimate = as.numeric(cf),
                   stringsAsFactors = FALSE)
  if (length(x$fit$bsMat) > 0 && !is.null(x$fit$var)) {
    se <- as.numeric(sqrt(diag(as.matrix(x$fit$var)))[df$term])
    df$std.error <- se
    df$statistic <- df$estimate / se
    df$p.value <- 2 * stats::pnorm(-abs(df$statistic))
  }
  if (exponentiate) df$estimate <- exp(df$estimate)
  df
}

#' Glance at a semi-parametric fit
#'
#' One-row model-level summary.
#'
#' @param x A `cens_sp_fit` object.
#' @param ... Unused.
#' @return A one-row data frame.
#' @examples
#' x <- as_interval_data(c(0, 1, 5, 8, 0, 2), c(1, 3, 5, 8, 1, 4))
#' grp <- factor(c("a", "a", "a", "b", "b", "b"))
#' glance(desc_sp(x, covariates = grp, model = "ph"))
#' @export
glance.cens_sp_fit <- function(x, ...) {
  data.frame(
    model = x$model,
    nobs = x$n,
    dropped = x$dropped,
    npar = length(stats::coef(x$fit)),
    logLik = tryCatch(as.numeric(x$fit$llk), error = function(e) NA_real_),
    bootstrap_samples = { nb <- nrow(x$fit$bsMat)
                          if (is.null(nb)) 0L else as.integer(nb) },
    stringsAsFactors = FALSE
  )
}
