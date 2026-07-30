#' @importFrom graphics lines legend
NULL

# ---- internal helpers ------------------------------------------------------

#' Coerce a `covariates` argument to a data frame
#'
#' @param covariates `NULL`, an atomic vector/factor, or a data frame.
#' @return A data frame, or `NULL`.
#' @keywords internal
#' @noRd
as_cov_df <- function(covariates) {
  if (is.null(covariates)) return(NULL)
  if (is.data.frame(covariates)) return(covariates)
  if (is.atomic(covariates)) return(data.frame(covariate = covariates))
  stop("`covariates` must be NULL, a vector/factor, or a data frame.",
       call. = FALSE)
}

#' Extract per-group survival step curves from a semi-parametric fit
#'
#' @param object A `cens_sp_fit` object.
#' @param newdata Optional data frame of covariate combinations.
#' @return A named list of data frames, each with `conc` and `surv`.
#' @keywords internal
#' @noRd
sp_surv_curves <- function(object, newdata = NULL) {
  sc <- if (is.null(newdata)) {
    icenReg::getSCurves(object$fit)
  } else {
    icenReg::getSCurves(object$fit, newdata = newdata)
  }
  conc <- as.numeric(sc$Tbull_ints[, "upper"])
  lapply(sc$S_curves, function(s) {
    data.frame(conc = conc, surv = as.numeric(s))
  })
}

# ---- fitting ---------------------------------------------------------------

#' Semi-parametric summary of interval-censored data
#'
#' Fits a semi-parametric proportional-hazards or proportional-odds model to a
#' `cens_data` object with a non-parametric baseline, via [icenReg::ic_sp()].
#' With no covariates this is a semi-parametric estimate of the baseline
#' distribution; with covariates it additionally estimates how each covariate
#' shifts the concentration distribution (as a hazard ratio under `"ph"` or an
#' odds ratio under `"po"`).
#'
#' Like [desc_np()], the model is fit on the original measurement scale, so
#' baseline curves and quantiles are expressed directly as concentrations.
#' Covariates are matched to `x` by position and are subset in step with any
#' observations dropped for missing bounds.
#'
#' @param x A `cens_data` object (see [as_interval_data()]).
#' @param covariates Optional covariates aligned to `x` (one row/element per
#'   observation): a data frame, an atomic vector, a factor, or `NULL` for an
#'   intercept-only baseline model.
#' @param model `"ph"` (proportional hazards, the default) or `"po"`
#'   (proportional odds).
#' @param bs_samples Number of bootstrap samples for coefficient standard errors
#'   and p-values. The default `0` fits point estimates only (fast and
#'   deterministic); set it higher for inference via [summary()].
#'
#' @return A `cens_sp_fit` object: a list with the fitted `fit`, the `model`,
#'   the `covariates` names, `log_transform`, the number of observations `n`,
#'   the count `dropped` for missing bounds, and the `category` counts.
#'
#' @seealso [plot.cens_sp_fit()] for baseline / per-group distribution curves.
#'
#' @examples
#' x <- as_interval_data(left = c(0, 1, 5, 8, 0, 2), right = c(1, 3, 5, 8, 1, 4))
#' grp <- factor(c("a", "a", "a", "b", "b", "b"))
#' fit <- desc_sp(x, covariates = grp, model = "ph")
#' fit
#' @export
desc_sp <- function(x, covariates = NULL, model = c("ph", "po"),
                    bs_samples = 0) {
  if (!is_cens_data(x)) {
    stop("`x` must be a <cens_data> object; see `as_interval_data()`.",
         call. = FALSE)
  }
  model <- match.arg(model)

  b <- cens_original_bounds(x)
  df <- data.frame(.lo = b$lo, .hi = b$hi)

  covdf <- as_cov_df(covariates)
  if (is.null(covdf)) {
    rhs <- "0"
    cov_names <- character(0)
  } else {
    if (nrow(covdf) != nrow(df)) {
      stop("`covariates` must have one row/element per observation in `x` ",
           sprintf("(got %d, expected %d).", nrow(covdf), nrow(df)),
           call. = FALSE)
    }
    df <- cbind(df, covdf)
    cov_names <- names(covdf)
    rhs <- paste(cov_names, collapse = " + ")
  }

  df <- df[b$keep, , drop = FALSE]
  if (nrow(df) == 0) {
    stop("No non-missing observations to fit.", call. = FALSE)
  }

  form <- stats::as.formula(paste("cbind(.lo, .hi) ~", rhs))
  fit <- icenReg::ic_sp(form, data = df, model = model,
                        bs_samples = bs_samples)

  structure(
    list(
      fit = fit,
      model = model,
      covariates = cov_names,
      log_transform = x$log_transform,
      n = nrow(df),
      dropped = b$dropped,
      category = table(factor(x$category, levels = cens_categories))
    ),
    class = "cens_sp_fit"
  )
}

# ---- methods ---------------------------------------------------------------

#' @export
coef.cens_sp_fit <- function(object, ...) {
  stats::coef(object$fit)
}

#' @export
summary.cens_sp_fit <- function(object, ...) {
  summary(object$fit)
}

#' Coerce a semi-parametric fit to a data frame
#'
#' Returns the covariate coefficient table as a tidy data frame --- `term`,
#' `coef` (log-scale), and the exponentiated `hazard_ratio` (PH) or `odds_ratio`
#' (PO), plus `std_error` and `p_value` when the fit was built with
#' `bs_samples > 0` --- ready to hand to a table-formatting package. An
#' intercept-only (no-covariate) fit yields a zero-row frame.
#'
#' @param x A `cens_sp_fit` object.
#' @param row.names,optional Passed along for S3 consistency.
#' @param ... Unused.
#'
#' @return A data frame, one row per covariate coefficient.
#' @examples
#' x <- as_interval_data(c(0, 1, 5, 8, 0, 2), c(1, 3, 5, 8, 1, 4))
#' grp <- factor(c("a", "a", "a", "b", "b", "b"))
#' as.data.frame(desc_sp(x, covariates = grp, model = "ph"))
#' @export
as.data.frame.cens_sp_fit <- function(x, row.names = NULL, optional = FALSE,
                                      ...) {
  cf <- stats::coef(x$fit)
  ratio_lab <- if (x$model == "ph") "hazard_ratio" else "odds_ratio"
  if (length(cf) == 0) {
    return(data.frame(term = character(0), coef = numeric(0),
                      stringsAsFactors = FALSE))
  }
  df <- data.frame(term = names(cf), coef = as.numeric(cf),
                   stringsAsFactors = FALSE)
  df[[ratio_lab]] <- exp(df$coef)
  if (length(x$fit$bsMat) > 0 && !is.null(x$fit$var)) {
    se <- sqrt(diag(as.matrix(x$fit$var)))[df$term]
    df$std_error <- as.numeric(se)
    df$p_value <- 2 * stats::pnorm(-abs(df$coef / df$std_error))
  }
  if (!is.null(row.names)) rownames(df) <- row.names
  df
}

#' Plot baseline / per-group distribution curves of a semi-parametric fit
#'
#' Draws the estimated cumulative distribution function `P(X <= conc)` on the
#' original measurement scale. With no `newdata`, the model baseline is drawn;
#' supply `newdata` (a data frame of covariate combinations) to overlay one
#' curve per group.
#'
#' @param x A `cens_sp_fit` object.
#' @param newdata Optional data frame of covariate combinations, one row per
#'   curve.
#' @param xlab,ylab Axis labels.
#' @param col Optional vector of colours, one per curve.
#' @param ... Further arguments passed to [plot()].
#'
#' @return `x`, invisibly.
#' @export
plot.cens_sp_fit <- function(x, newdata = NULL, xlab = "Concentration",
                             ylab = "Cumulative probability", col = NULL, ...) {
  curves <- sp_surv_curves(x, newdata)
  k <- length(curves)
  labels <- if (!is.null(newdata)) {
    apply(newdata, 1L, function(r) {
      paste(names(newdata), r, sep = "=", collapse = ", ")
    })
  } else {
    names(curves)
  }
  if (is.null(col)) col <- seq_len(k)

  all_conc <- unlist(lapply(curves, `[[`, "conc"))
  plot(NA, xlim = range(all_conc, finite = TRUE), ylim = c(0, 1),
       xlab = xlab, ylab = ylab, ...)
  for (i in seq_len(k)) {
    cv <- curves[[i]]
    ord <- order(cv$conc)
    conc <- cv$conc[ord]
    lines(c(conc[1], conc), c(0, 1 - cv$surv[ord]), type = "s", col = col[i])
  }
  if (!is.null(newdata) || k > 1) {
    legend("bottomright", legend = labels, col = col, lty = 1, bty = "n")
  }
  invisible(x)
}

#' @export
print.cens_sp_fit <- function(x, ...) {
  model_lab <- switch(
    x$model,
    ph = "proportional hazards (Cox PH)",
    po = "proportional odds"
  )
  cat(sprintf("<cens_sp_fit> semi-parametric %s\n", model_lab))
  cat(sprintf("  %d observation%s", x$n, if (x$n == 1L) "" else "s"))
  if (x$dropped > 0) cat(sprintf(" (%d dropped for missing bounds)", x$dropped))
  cat("\n")
  for (nm in cens_categories) {
    cat(sprintf("  %-24s %d\n", nm, x$category[[nm]]))
  }

  cf <- stats::coef(x$fit)
  if (length(cf) == 0) {
    cat("\n  Intercept-only baseline (no covariates).\n")
  } else {
    ratio_lab <- if (x$model == "ph") "hazard_ratio" else "odds_ratio"
    tab <- data.frame(
      coef = as.numeric(cf),
      ratio = exp(as.numeric(cf)),
      row.names = names(cf)
    )
    names(tab)[2] <- ratio_lab
    cat("\nCoefficients:\n")
    print(tab)
    if (length(x$fit$bsMat) == 0) {
      cat("(no standard errors; refit with bs_samples > 0 for inference)\n")
    }
  }
  invisible(x)
}
