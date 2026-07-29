# Gibbs-sampler imputation for interval-censored data.
#
# This is an independent (clean-room) implementation of the GSimp imputation
# method described in Wei et al. (2018), "GSimp: A Gibbs sampler based
# left-censored missing value imputation approach for metabolomics studies"
# (PLoS Comput. Biol. 14(1):e1005973, doi:10.1371/journal.pcbi.1005973). It is
# written from the published algorithm, NOT ported from the GSimp source code
# (which is CC BY-NC-SA licensed); the statistical method itself is not
# copyrightable. This keeps `leftcens` under its MIT license.
#
# The extension over the original method is per-cell truncation bounds
# (lo_mat / hi_mat), so a single analyte can carry non-detect (0, MDL) and
# detected-not-quantified (MDL, LCMRL) cells simultaneously.

# ---- truncated-normal draw -------------------------------------------------

#' Draw from a truncated normal distribution
#'
#' Samples from `N(mean, sd^2)` conditioned to lie in `(lower, upper)`, using
#' inverse-CDF sampling. One-sided-infinite intervals (e.g. the `(-Inf, log MDL)`
#' bound of a log-scale non-detect) are handled on the log-probability scale so
#' draws deep in a tail do not underflow to `+/-Inf`. Arguments are recycled to
#' length `n`.
#'
#' @param n Number of draws.
#' @param mean,sd Mean and standard deviation of the untruncated normal.
#' @param lower,upper Truncation bounds (may be `-Inf` / `Inf`).
#'
#' @return A numeric vector of length `n`.
#' @examples
#' # left-censored: values known only to be below log(2)
#' summary(rnorm_trunc(1000, mean = 1, sd = 1, lower = -Inf, upper = log(2)))
#' @export
rnorm_trunc <- function(n, mean, sd, lower, upper) {
  mean <- rep_len(mean, n)
  sd <- rep_len(sd, n)
  lower <- rep_len(lower, n)
  upper <- rep_len(upper, n)

  out <- numeric(n)
  if (n == 0) return(out)

  # Degenerate cases first.
  degen <- !is.finite(sd) | sd <= 0 | lower == upper
  if (any(degen)) {
    out[degen] <- pmin(pmax(mean[degen], lower[degen]), upper[degen])
  }
  live <- !degen
  if (!any(live)) return(out)

  a <- (lower - mean) / sd
  b <- (upper - mean) / sd
  u <- stats::runif(n)
  # Keep u away from the exact endpoints so log(u) and qnorm() stay finite.
  u <- pmin(pmax(u, 1e-300), 1 - 1e-16)

  left_inf <- is.infinite(a) & a < 0
  right_inf <- is.infinite(b) & b > 0

  z <- numeric(n)

  # Both infinite -> plain standard normal.
  m <- live & left_inf & right_inf
  if (any(m)) z[m] <- stats::qnorm(u[m])

  # Lower = -Inf, upper finite: draw below b via the lower tail on log scale.
  m <- live & left_inf & !right_inf
  if (any(m)) {
    z[m] <- stats::qnorm(log(u[m]) + stats::pnorm(b[m], log.p = TRUE),
                         log.p = TRUE)
  }

  # Lower finite, upper = +Inf: mirror the lower-tail trick.
  m <- live & !left_inf & right_inf
  if (any(m)) {
    z[m] <- -stats::qnorm(log(u[m]) + stats::pnorm(-a[m], log.p = TRUE),
                          log.p = TRUE)
  }

  # Both finite: standard inverse-CDF, with a midpoint guard when the
  # probability mass in the interval underflows to zero.
  m <- live & !left_inf & !right_inf
  if (any(m)) {
    Fa <- stats::pnorm(a[m])
    Fb <- stats::pnorm(b[m])
    uu <- Fa + u[m] * (Fb - Fa)
    zz <- stats::qnorm(uu)
    collapsed <- !is.finite(zz) | (Fb - Fa) <= 0
    zz[collapsed] <- (a[m] + b[m])[collapsed] / 2
    z[m] <- zz
  }

  out[live] <- mean[live] + sd[live] * z[live]
  out
}

# ---- conditional prediction models -----------------------------------------
#
# Each model is a function(yo, Xo, Xm) returning list(mean, sd): the predicted
# conditional mean at the missing rows (Xm) and a residual sd, both estimated
# from the observed rows (yo ~ Xo).

#' Ridge-regression conditional model (default, dependency-free)
#' @keywords internal
#' @noRd
ridge_predict <- function(yo, Xo, Xm, lambda = 1e-3) {
  if (ncol(Xo) == 0) {
    return(list(mean = rep(mean(yo), nrow(Xm)), sd = stats::sd(yo)))
  }
  ctr <- colMeans(Xo)
  scl <- apply(Xo, 2, stats::sd)
  scl[!is.finite(scl) | scl == 0] <- 1
  Xoc <- sweep(sweep(Xo, 2, ctr, "-"), 2, scl, "/")
  Xmc <- sweep(sweep(Xm, 2, ctr, "-"), 2, scl, "/")
  ybar <- mean(yo)
  yc <- yo - ybar
  p <- ncol(Xoc)
  beta <- solve(crossprod(Xoc) + diag(lambda, p), crossprod(Xoc, yc))
  pred <- as.numeric(Xmc %*% beta) + ybar
  resid <- yc - as.numeric(Xoc %*% beta)
  sdr <- sqrt(sum(resid^2) / max(1L, length(yo) - 1L))
  if (!is.finite(sdr) || sdr <= 0) sdr <- stats::sd(yo)
  list(mean = pred, sd = sdr)
}

#' Ordinary-least-squares conditional model
#' @keywords internal
#' @noRd
lm_predict <- function(yo, Xo, Xm) {
  if (ncol(Xo) == 0) {
    return(list(mean = rep(mean(yo), nrow(Xm)), sd = stats::sd(yo)))
  }
  cn <- paste0("x", seq_len(ncol(Xo)))
  dfo <- as.data.frame(Xo)
  names(dfo) <- cn
  dfo$.y <- yo
  fit <- tryCatch(stats::lm(.y ~ ., data = dfo), error = function(e) NULL)
  if (is.null(fit)) return(ridge_predict(yo, Xo, Xm))
  dfm <- as.data.frame(Xm)
  names(dfm) <- cn
  pred <- as.numeric(stats::predict(fit, dfm))
  sdr <- stats::sigma(fit)
  if (!is.finite(sdr) || sdr <= 0) sdr <- stats::sd(yo)
  list(mean = pred, sd = sdr)
}

#' Elastic-net conditional model (matches the original GSimp default)
#' @keywords internal
#' @noRd
glmnet_predict <- function(yo, Xo, Xm) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("imp_model = \"glmnet\" requires the 'glmnet' package.", call. = FALSE)
  }
  if (ncol(Xo) < 2) return(ridge_predict(yo, Xo, Xm))
  cvfit <- glmnet::cv.glmnet(Xo, yo, alpha = 0.5)
  pred <- as.numeric(stats::predict(cvfit, newx = Xm, s = "lambda.min"))
  fitted_o <- as.numeric(stats::predict(cvfit, newx = Xo, s = "lambda.min"))
  sdr <- stats::sd(yo - fitted_o)
  if (!is.finite(sdr) || sdr <= 0) sdr <- stats::sd(yo)
  list(mean = pred, sd = sdr)
}

#' @keywords internal
#' @noRd
resolve_imp_model <- function(imp_model) {
  if (is.function(imp_model)) return(imp_model)
  switch(
    imp_model,
    ridge = ridge_predict,
    lm = lm_predict,
    glmnet = glmnet_predict,
    stop("`imp_model` must be \"ridge\", \"lm\", \"glmnet\", or a function.",
         call. = FALSE)
  )
}

# ---- initialisation --------------------------------------------------------

#' Initialise censored cells to a valid in-bounds starting point
#' @keywords internal
#' @noRd
initialise_cells <- function(data, lo, hi, miss, initial) {
  cur <- data
  p <- ncol(data)

  if (identical(initial, "qrilc")) {
    if (!requireNamespace("imputeLCMD", quietly = TRUE)) {
      stop("initial = \"qrilc\" requires the 'imputeLCMD' package.",
           call. = FALSE)
    }
    filled <- imputeLCMD::impute.QRILC(data)[[1]]
    # Clamp QRILC draws into each cell's interval.
    cur[miss] <- pmin(pmax(filled[miss], lo[miss]), hi[miss])
    return(cur)
  }

  if (!identical(initial, "bounds")) {
    stop("`initial` must be \"bounds\" or \"qrilc\".", call. = FALSE)
  }

  all_obs <- data[!miss]
  gm <- mean(all_obs)
  gs <- stats::sd(all_obs)
  if (!is.finite(gs) || gs <= 0) gs <- 1
  for (j in seq_len(p)) {
    mj <- which(miss[, j])
    if (length(mj) == 0) next
    obs <- data[!miss[, j], j]
    if (length(obs) >= 2 && is.finite(stats::sd(obs)) && stats::sd(obs) > 0) {
      m <- mean(obs); s <- stats::sd(obs)
    } else {
      m <- gm; s <- gs
    }
    cur[mj, j] <- rnorm_trunc(length(mj), m, s, lo[mj, j], hi[mj, j])
  }
  cur
}

# ---- main sampler ----------------------------------------------------------

#' Gibbs-sampler imputation of interval-censored data
#'
#' Imputes the censored cells of a wide analyte matrix by Gibbs sampling: each
#' analyte is repeatedly redrawn from its conditional distribution given the
#' other analytes, with every draw constrained to that cell's censoring interval
#' `(lo, hi)` via a truncated normal ([rnorm_trunc()]). This is a clean-room
#' implementation of the GSimp method (Wei et al. 2018), extended to per-cell
#' bounds so non-detect and detected-not-quantified cells of the same analyte
#' are imputed within their respective bands.
#'
#' @param x A `cens_bounds` object from [build_bounds()], which supplies the
#'   working-scale data matrix and the per-cell `(lo, hi)` bound matrices.
#' @param iters_all Number of outer Gibbs sweeps over all analytes (the mixing /
#'   convergence control). Each sweep updates every analyte once.
#' @param iters_each Inner refinement draws per analyte per sweep. This matters
#'   for stochastic conditional models (e.g. cross-validated `glmnet`); for the
#'   deterministic default it is redundant, hence the default of `1`.
#' @param initial Initialisation of censored cells: `"bounds"` (default; draw
#'   from each analyte's observed mean/sd, truncated to the cell interval) or
#'   `"qrilc"` (quantile-regression left-censored imputation, requires the
#'   'imputeLCMD' package, then clamped into bounds).
#' @param imp_model Conditional model: `"ridge"` (default, base-R ridge
#'   regression), `"lm"`, `"glmnet"` (elastic net, requires 'glmnet'), or a
#'   custom `function(yo, Xo, Xm)` returning `list(mean, sd)`.
#' @param n_cores Reserved for future parallel sweeps; currently only serial
#'   execution is implemented (a value `> 1` emits a message and runs serially).
#' @param verbose If `TRUE`, report progress per sweep.
#'
#' @return A numeric matrix of the same dimensions as the input, with censored
#'   cells filled in, on the same (working / log) scale as the `cens_bounds`
#'   input. Observed cells are returned unchanged.
#'
#' @references Wei R. et al. (2018) GSimp: A Gibbs sampler based left-censored
#'   missing value imputation approach for metabolomics studies. PLoS Comput.
#'   Biol. 14(1):e1005973. \doi{10.1371/journal.pcbi.1005973}
#'
#' @seealso [build_bounds()] for constructing the input.
#'
#' @examples
#' a <- as_interval_data(c(0, 2, 5, 8, 3), c(1, 4, 5, 8, 3))
#' b <- as_interval_data(c(3, 0, 7, 2, 6), c(3, 1, 7, 4, 6))
#' bnds <- build_bounds(list(A = a, B = b))
#' filled <- gsimp_impute(bnds, iters_all = 5)
#' anyNA(filled)
#' @export
gsimp_impute <- function(x, iters_all = 10, iters_each = 1,
                         initial = "bounds", imp_model = "ridge",
                         n_cores = 1, verbose = FALSE) {
  if (!is_cens_bounds(x)) {
    stop("`x` must be a <cens_bounds> object; see `build_bounds()`.",
         call. = FALSE)
  }
  if (n_cores > 1) {
    message("n_cores > 1 is not yet supported; running serially.")
  }

  data <- x$data_wide
  lo <- x$lo_mat
  hi <- x$hi_mat
  miss <- x$to_impute
  p <- ncol(data)

  if (!any(miss)) return(data)

  model_fn <- resolve_imp_model(imp_model)
  cur <- initialise_cells(data, lo, hi, miss, initial)

  # Update better-observed analytes first, so early sweeps use good predictors.
  col_order <- order(colSums(miss))

  for (sweep_i in seq_len(iters_all)) {
    for (j in col_order) {
      mj <- which(miss[, j])
      if (length(mj) == 0) next
      oj <- which(!miss[, j])
      Xall <- cur[, -j, drop = FALSE]

      if (length(oj) < 2 || p < 2) {
        # Too little to condition on: redraw from the analyte's own summary.
        m <- if (length(oj) >= 2) mean(cur[oj, j]) else mean(cur[, j])
        s <- if (length(oj) >= 2) stats::sd(cur[oj, j]) else stats::sd(cur[, j])
        if (!is.finite(s) || s <= 0) s <- 1
        cur[mj, j] <- rnorm_trunc(length(mj), m, s, lo[mj, j], hi[mj, j])
        next
      }

      for (it in seq_len(iters_each)) {
        pr <- model_fn(cur[oj, j], Xall[oj, , drop = FALSE],
                       Xall[mj, , drop = FALSE])
        cur[mj, j] <- rnorm_trunc(length(mj), pr$mean, pr$sd,
                                  lo[mj, j], hi[mj, j])
      }
    }
    if (verbose) message(sprintf("gsimp_impute: sweep %d/%d done",
                                 sweep_i, iters_all))
  }

  cur
}
