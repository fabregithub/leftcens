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

#' Censored (Tobit) conditional model
#'
#' Unlike the observed-only models, this fits the analyte's conditional
#' distribution on ALL rows using an interval-censored Gaussian regression
#' (`survival::survreg`, `dist = "gaussian"`): observed rows enter as exact
#' values, censored rows as their `(lo, hi)` interval. This removes the
#' selection bias that arises from fitting only on the (upper-truncated)
#' observed values.
#'
#' When there are more conditioning analytes than the fit can support (wide
#' data, `p >= n`), the predictors are reduced to their leading principal
#' components before the censored fit, rather than falling back to ridge. Ridge
#' is used only when there is nothing to condition on, too few observed rows, or
#' the censored fit cannot be formed.
#'
#' Takes the full per-column context list (see `gsimp_impute`'s loop), not the
#' `(yo, Xo, Xm)` triple the observed-only models use.
#' @param max_pc Maximum number of principal components retained for the
#'   wide-data fit.
#' @keywords internal
#' @noRd
tobit_predict <- function(ctx, max_pc = 20L) {
  X <- ctx$X_all
  n <- nrow(X)
  pX <- ncol(X)
  ridge_fallback <- function() {
    ridge_predict(ctx$y_obs, ctx$X_obs, ctx$X_mis)
  }
  if (pX == 0 || sum(ctx$obs) < 3) return(ridge_fallback())

  # Interval-censored response on the working (log) scale: exact where observed,
  # (lo, hi) where censored. survreg's interval2 uses NA for an open endpoint.
  left <- ifelse(ctx$obs, ctx$y_all, ctx$lo_j)
  right <- ifelse(ctx$obs, ctx$y_all, ctx$hi_j)
  left[is.infinite(left) & left < 0] <- NA     # -Inf lower -> left-censored
  right[is.infinite(right) & right > 0] <- NA  #  Inf upper -> right-censored
  resp <- survival::Surv(left, right, type = "interval2")

  # Design matrix: the predictors directly when the fit can support them,
  # otherwise their leading principal components (handles p >= n wide data).
  # Cap the dimension well below n so the censored fit is not overparameterised.
  k <- min(pX, max(1L, n %/% 3L), max_pc)
  if (pX <= k) {
    D_all <- X
  } else {
    pcs <- tryCatch(
      stats::prcomp(X, center = TRUE, scale. = FALSE, rank. = k),
      error = function(e) NULL
    )
    if (is.null(pcs)) return(ridge_fallback())
    D_all <- pcs$x[, seq_len(min(k, ncol(pcs$x))), drop = FALSE]
  }

  cn <- paste0("v", seq_len(ncol(D_all)))
  fit_df <- as.data.frame(D_all)
  names(fit_df) <- cn

  fit <- tryCatch(
    suppressWarnings(
      survival::survreg(resp ~ ., data = fit_df, dist = "gaussian")
    ),
    error = function(e) NULL
  )
  if (is.null(fit) || !is.finite(fit$scale) || fit$scale <= 0) {
    return(ridge_fallback())
  }

  new_df <- as.data.frame(D_all[ctx$mis_idx, , drop = FALSE])
  names(new_df) <- cn
  lp <- as.numeric(stats::predict(fit, newdata = new_df, type = "lp"))
  if (any(!is.finite(lp))) return(ridge_fallback())
  list(mean = lp, sd = fit$scale)
}

#' Resolve `imp_model` to a function taking the per-column context list.
#'
#' Observed-only models (`ridge`/`lm`/`glmnet`) and user-supplied functions keep
#' the `(yo, Xo, Xm)` contract and are wrapped to read those from the context;
#' `tobit` receives the full context (it needs the censored rows' bounds).
#' @keywords internal
#' @noRd
resolve_imp_model <- function(imp_model) {
  if (is.function(imp_model)) {
    return(function(ctx) imp_model(ctx$y_obs, ctx$X_obs, ctx$X_mis))
  }
  switch(
    imp_model,
    ridge = function(ctx) ridge_predict(ctx$y_obs, ctx$X_obs, ctx$X_mis),
    lm = function(ctx) lm_predict(ctx$y_obs, ctx$X_obs, ctx$X_mis),
    glmnet = function(ctx) glmnet_predict(ctx$y_obs, ctx$X_obs, ctx$X_mis),
    tobit = tobit_predict,
    stop("`imp_model` must be \"ridge\", \"lm\", \"glmnet\", \"tobit\", ",
         "or a function.", call. = FALSE)
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
#' @param imp_model Conditional model. `"copula"` (default) is the skew-robust
#'   choice: it fits a flexible (sinh-arcsinh) margin per analyte, maps the data
#'   to latent-Gaussian scores, runs an interval-censored Gaussian regression
#'   there, and back-transforms -- keeping the mean near-unbiased under
#'   right-skewed margins while matching `"tobit"` on log-normal data. It is a
#'   *proper* multiple-imputation engine (see `margin_draw`), and the slowest
#'   option. `"tobit"` is an interval-censored Gaussian regression via
#'   [survival::survreg()] fit on observed *and* censored rows; fast-ish and
#'   reliable *when the margins are log-normal*, but biased under skew -- use it
#'   for known-symmetric data or very large/wide problems. The observed-only
#'   alternatives are `"ridge"` (base-R ridge regression, fastest and always
#'   solvable), `"lm"`, and `"glmnet"` (elastic net, requires 'glmnet'). A custom
#'   `function(yo, Xo, Xm)` returning `list(mean, sd)` is also accepted. `"tobit"`
#'   (and the latent step of `"copula"`) reduces the predictor dimension by PCA
#'   when there are more analytes than samples, and falls back to `"ridge"` if the
#'   censored fit cannot be formed.
#' @param margin_draw Only used by `imp_model = "copula"`. When `TRUE` (default),
#'   each imputation draws the per-analyte margin parameters from their asymptotic
#'   posterior, so pooling M such imputations propagates margin/tail-estimation
#'   uncertainty into the multiple-imputation variance. Set `FALSE` to plug in the
#'   MLE margin, giving the near-unbiased *point* imputation. Recommended MI
#'   workflow: the plug-in (`FALSE`) for the point estimate, drawn (`TRUE`)
#'   imputations for the between-imputation variance.
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
                         initial = "bounds", imp_model = "copula",
                         n_cores = 1, verbose = FALSE, margin_draw = TRUE) {
  if (!is_cens_bounds(x)) {
    stop("`x` must be a <cens_bounds> object; see `build_bounds()`.",
         call. = FALSE)
  }
  if (n_cores > 1) {
    message("n_cores > 1 is not yet supported; running serially.")
  }

  # Gaussian-copula path: transform to latent-Gaussian margins, impute there with
  # tobit, back-transform. See R/impute-copula.R.
  if (identical(imp_model, "copula")) {
    return(copula_impute(x, iters_all = iters_all, iters_each = iters_each,
                         initial = initial, verbose = verbose,
                         margin_draw = margin_draw))
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

      obs_j <- !miss[, j]
      for (it in seq_len(iters_each)) {
        ctx <- list(
          y_obs = cur[oj, j],
          X_obs = Xall[oj, , drop = FALSE],
          X_mis = Xall[mj, , drop = FALSE],
          X_all = Xall,
          y_all = cur[, j],
          obs = obs_j,
          lo_j = lo[, j],
          hi_j = hi[, j],
          mis_idx = mj
        )
        pr <- model_fn(ctx)
        cur[mj, j] <- rnorm_trunc(length(mj), pr$mean, pr$sd,
                                  lo[mj, j], hi[mj, j])
      }
    }
    if (verbose) message(sprintf("gsimp_impute: sweep %d/%d done",
                                 sweep_i, iters_all))
  }

  cur
}
