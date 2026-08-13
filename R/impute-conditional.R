# =============================================================================
# General-predictor, interval-censored, skew-aware conditional imputation.
# -----------------------------------------------------------------------------
# The single-variable core behind congenial imputation of a censored variable on
# an arbitrary design matrix (e.g. an outcome Y, mixed-type covariates Z, and
# other analytes). Unlike the copula imputer -- which conditions an analyte on the
# OTHER ANALYTES only -- this conditions on a general set of predictors, so it can
# be used inside a fully-conditional-specification (FCS) loop or as the engine of
# a `mice.impute.*` method.
#
# It reuses the sinh-arcsinh margin (R/shash-margin.R): fit the margin, map the
# variable and its censoring bounds to a latent normal scale, run a linear-
# Gaussian interval-censored regression there (where it is correctly specified),
# draw truncated within the latent bounds, and map back. Proper MI draws both the
# margin parameters and the regression coefficients per imputation.
# =============================================================================

#' Impute a censored variable conditional on a general design matrix
#'
#' Multiple imputation of a single left-/interval-censored variable, conditioning
#' on an arbitrary predictor matrix and respecting the per-cell censoring bounds.
#' Skew is handled by a sinh-arcsinh margin so non-detects are drawn from a
#' correctly shaped tail; the conditional mean is modelled linearly on the latent
#' normal scale. This is the general-predictor complement to
#' [gsimp_impute()] (which conditions on other analytes only) and the reusable
#' core for FCS / `mice` integration.
#'
#' @param y Numeric vector; `NA` at the cells to impute (censored or missing),
#'   the observed value elsewhere.
#' @param x Predictors: a matrix or data frame with one row per element of `y`
#'   (an outcome, covariates, other analytes; mixed types allowed). No intercept
#'   column is needed -- one is added.
#' @param lower,upper Per-cell censoring bounds, same length as `y`, used only at
#'   the `NA` cells: the imputed value is drawn within `(lower, upper]`. For a
#'   left-censored non-detect at limit `L` use `lower = -Inf`, `upper = L`; for an
#'   interval use both finite bounds. Values at observed cells are ignored.
#' @param m Number of completed vectors to return.
#' @param margin `"shash"` (default, skew-aware sinh-arcsinh margin) or
#'   `"gaussian"` (a plain tobit conditional; use when skew is known negligible).
#' @param proper If `TRUE` (default), draw the margin and regression parameters
#'   from their posteriors per imputation, so between-imputation variance is
#'   honest. `FALSE` plugs in the MLEs (single-imputation / point use).
#' @return A numeric matrix with `length(y)` rows and `m` columns: each column is
#'   a completed copy of `y` (observed cells carried through unchanged, censored
#'   cells imputed).
#' @seealso [fit_shash_margin()], [draw_margin()], [rnorm_trunc()], [gsimp_impute()].
#' @examples
#' set.seed(1)
#' n <- 300
#' z  <- rnorm(n)
#' x1 <- 0.6 * z + rnorm(n)                     # a covariate-correlated exposure
#' lod <- quantile(x1, 0.3)
#' y  <- x1                                     # observed values
#' y[x1 < lod] <- NA                            # left-censored below the LOD
#' lower <- ifelse(is.na(y), -Inf, NA)
#' upper <- ifelse(is.na(y), lod,  NA)
#' imp <- impute_censored_conditional(y, x = data.frame(z), lower, upper, m = 5)
#' dim(imp)                                     # n x 5 completed exposures
#' @export
impute_censored_conditional <- function(y, x, lower, upper, m = 1L,
                                        margin = c("shash", "gaussian"),
                                        proper = TRUE) {
  margin <- match.arg(margin)
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required for impute_censored_conditional().",
         call. = FALSE)
  }
  n <- length(y)
  xdf <- as.data.frame(x)
  if (nrow(xdf) != n) {
    stop("`x` must have one row per element of `y` (got ",
         nrow(xdf), " rows for length-", n, " `y`).", call. = FALSE)
  }
  if (length(lower) != n || length(upper) != n) {
    stop("`lower` and `upper` must have the same length as `y`.", call. = FALSE)
  }
  m <- as.integer(m)
  if (!is.finite(m) || m < 1L) stop("`m` must be a positive integer.", call. = FALSE)

  ry <- !is.na(y); wy <- !ry
  if (!any(wy)) {                                # nothing to impute
    return(matrix(y, nrow = n, ncol = m))
  }
  if (!any(ry)) {
    stop("`y` has no observed (non-NA) values to fit the imputation model.",
         call. = FALSE)
  }

  # Margin: skew-aware sinh-arcsinh, or an identity transform (-> plain tobit).
  mfit <- if (margin == "shash") {
    fit_shash_margin(y[ry], lower[wy], upper[wy])
  } else {
    list(mu = 0, sigma = 1, eps = 0, par = c(0, 0, 0), V = NULL)
  }

  pred_names <- names(xdf)
  fm <- stats::as.formula(paste(
    "survival::Surv(.t1, .t2, type = 'interval2') ~",
    if (length(pred_names)) paste(sprintf("`%s`", pred_names), collapse = " + ") else "1"))

  out <- matrix(NA_real_, nrow = n, ncol = m)
  for (i in seq_len(m)) {
    mp <- if (proper) draw_margin(mfit) else list(mu = mfit$mu, sigma = mfit$sigma, eps = mfit$eps)

    z_obs <- x_to_z(y,     mp$mu, mp$sigma, mp$eps)   # NA at wy
    z_lo  <- x_to_z(lower, mp$mu, mp$sigma, mp$eps)
    z_hi  <- x_to_z(upper, mp$mu, mp$sigma, mp$eps)

    # interval2 coding: exact where observed; open ends coded NA.
    t1 <- z_obs; t2 <- z_obs
    t1[wy] <- ifelse(is.finite(z_lo[wy]), z_lo[wy], NA_real_)   # -Inf -> left-censored
    t2[wy] <- ifelse(is.finite(z_hi[wy]), z_hi[wy], NA_real_)   #  Inf -> right-censored

    dat <- cbind(data.frame(.t1 = t1, .t2 = t2), xdf)
    sr <- tryCatch(survival::survreg(fm, data = dat, dist = "gaussian"),
                   error = function(e) NULL)
    if (is.null(sr)) { out[, i] <- z_to_x(z_obs, mp$mu, mp$sigma, mp$eps); next }

    XX <- stats::model.matrix(sr)
    beta <- stats::coef(sr); k <- length(beta)
    if (proper) {                                # draw (beta, log scale) jointly
      V <- stats::vcov(sr); par <- c(beta, log(sr$scale))
      drawn <- tryCatch(as.numeric(par + t(chol(V)) %*% stats::rnorm(length(par))),
                        error = function(e) par)
      beta_i <- drawn[seq_len(k)]; s_i <- exp(drawn[k + 1L])
    } else {
      beta_i <- beta; s_i <- sr$scale
    }
    mu_lin <- as.vector(XX %*% beta_i)

    z1 <- z_obs
    z1[wy] <- rnorm_trunc(sum(wy), mu_lin[wy], s_i,
                          lower = z_lo[wy], upper = z_hi[wy])
    out[, i] <- z_to_x(z1, mp$mu, mp$sigma, mp$eps)
  }
  out
}
