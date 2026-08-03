# Gaussian-copula ("normal-scores") imputation for interval-censored data.
#
# S1 in validation/PHASE3_SKEW_PLAN.md. Motivation: the default tobit conditional
# assumes a Gaussian margin, so under right-skew it draws non-detects from a
# too-heavy left tail and biases the mean downward (D1 / G2 findings). This model
# instead treats the data as a Gaussian COPULA with flexible skewed margins:
#
#     X_j = F_j^{-1}( Phi(Z_j) ),   Z = (Z_1..Z_p) ~ N(0, Sigma).
#
# Skew lives in the margins F_j; dependence lives in Sigma. We estimate each F_j
# from the (interval-censored) data, map observed values to latent normal scores
# Z_j = Phi^{-1}(F_j(X_j)) and censoring bounds to latent thresholds, then run the
# EXISTING tobit Gibbs sampler in Z-space -- where the linear-Gaussian conditional
# and the truncated-normal draw are now correctly specified -- and back-transform
# X_j = F_j^{-1}(Phi(Z_j)). Correlation then helps, because the latent conditional
# variance shrinks as dependence grows, pinning a censored Z_j into a narrower band.
#
# Margin family: a 3-parameter sinh-arcsinh transform-to-normality
#   X = mu + sigma * sinh( asinh(Z) + eps ),   Z ~ N(0,1)
# (location mu, scale sigma, skewness eps), fit per analyte by interval-censored
# MLE. eps = 0 reduces to a plain Gaussian margin, so the method degrades
# gracefully to tobit when no skew is detected or the fit fails.
#
# NOTE (evaluation): sinh-arcsinh margins matched to a sinh-arcsinh DGP is an
# "inverse crime" -- favourable by construction. The honest test is G5 (a
# different margin family in the DGP); see the plan. The family here is swappable.

# ---- the sinh-arcsinh margin transform -------------------------------------

# data -> latent normal score (monotone increasing; exact inverse of z_to_x)
x_to_z <- function(x, mu, sigma, eps) sinh(asinh((x - mu) / sigma) - eps)
# latent normal score -> data (the marginal quantile map F^{-1}(Phi(.)))
z_to_x <- function(z, mu, sigma, eps) mu + sigma * sinh(asinh(z) + eps)

#' Fit a 3-parameter sinh-arcsinh margin by interval-censored MLE.
#'
#' @param xo Observed (exact) values for the analyte.
#' @param lo_c,hi_c Lower/upper bounds of the censored cells (`lo_c` may be
#'   `-Inf` for non-detects; `hi_c` may be `Inf` for right-censored).
#' @return list(mu, sigma, eps, par, V): the point estimate plus the fitted
#'   parameter vector `par = (mu, log sigma, eps)` and its asymptotic covariance
#'   `V` (inverse observed information). `V` is `NULL` when it cannot be formed,
#'   signalling [draw_margin()] to fall back to the plug-in (no parameter draw).
#'   Falls back to (mean, sd, 0) -- a plain Gaussian margin -- if there is too
#'   little data or the optimiser fails.
#' @keywords internal
#' @noRd
fit_shash_margin <- function(xo, lo_c, hi_c) {
  fallback <- list(mu = mean(xo), sigma = stats::sd(xo), eps = 0,
                   par = NULL, V = NULL)
  if (!is.finite(fallback$sigma) || fallback$sigma <= 0) fallback$sigma <- 1
  fallback$par <- c(fallback$mu, log(fallback$sigma), 0)
  if (length(xo) < 5L) return(fallback)

  nll <- function(par) {
    mu <- par[1]; sigma <- exp(par[2]); eps <- par[3]
    y <- (xo - mu) / sigma
    a <- asinh(y)
    z <- sinh(a - eps)
    # log density of an exact observation (normal density x Jacobian)
    ll <- sum(stats::dnorm(z, log = TRUE) + log(cosh(a - eps)) -
                par[2] - 0.5 * log1p(y^2))
    if (length(hi_c)) {                       # interval-censored contributions
      zlo <- ifelse(is.infinite(lo_c) & lo_c < 0, -Inf,
                    sinh(asinh((lo_c - mu) / sigma) - eps))
      zhi <- ifelse(is.infinite(hi_c) & hi_c > 0, Inf,
                    sinh(asinh((hi_c - mu) / sigma) - eps))
      pr <- stats::pnorm(zhi) - stats::pnorm(zlo)
      ll <- ll + sum(log(pmax(pr, 1e-12)))
    }
    if (!is.finite(ll)) return(1e10)
    -ll
  }

  start <- c(fallback$mu, log(fallback$sigma), 0)
  fit <- tryCatch(
    stats::optim(start, nll, method = "Nelder-Mead",
                 control = list(maxit = 500, reltol = 1e-8)),
    error = function(e) NULL
  )
  if (is.null(fit) || fit$convergence != 0 || !all(is.finite(fit$par))) {
    return(fallback)
  }
  sigma <- exp(fit$par[2])
  if (!is.finite(sigma) || sigma <= 0 || abs(fit$par[3]) > 10) return(fallback)

  # Asymptotic covariance of (mu, log sigma, eps): inverse of the observed
  # information (Hessian of the negative log-likelihood). Kept only if finite and
  # positive-definite; otherwise V = NULL -> plug-in imputation (no draw).
  V <- tryCatch({
    H <- stats::optimHess(fit$par, nll)
    Vh <- solve(H)
    chol(Vh)                                  # errors unless positive-definite
    Vh
  }, error = function(e) NULL)
  if (!is.null(V) && !all(is.finite(V))) V <- NULL

  list(mu = fit$par[1], sigma = sigma, eps = fit$par[3], par = fit$par, V = V)
}

#' Draw a margin's parameters from their asymptotic posterior (proper MI).
#'
#' Each multiple imputation should condition on a *drawn* margin, not the shared
#' MLE, so the between-imputation variance reflects margin- and (extrapolated)
#' tail-estimation uncertainty -- the fix for the copula's coverage dip under
#' heavy censoring. Returns the plug-in estimate when `m$V` is `NULL` or the draw
#' is degenerate, so behaviour is unchanged where the fit is untrustworthy.
#' @keywords internal
#' @noRd
draw_margin <- function(m) {
  plugin <- list(mu = m$mu, sigma = m$sigma, eps = m$eps)
  if (is.null(m$V)) return(plugin)
  d <- tryCatch(as.numeric(m$par + t(chol(m$V)) %*% stats::rnorm(3)),
                error = function(e) NULL)
  if (is.null(d) || !all(is.finite(d)) || abs(d[3]) > 10) return(plugin)
  s <- exp(d[2])
  if (!is.finite(s) || s <= 0) return(plugin)
  list(mu = d[1], sigma = s, eps = d[3])
}

# ---- the copula imputer ----------------------------------------------------

#' Gaussian-copula imputation (reached via `gsimp_impute(imp_model = "copula")`).
#'
#' Estimates a skewed margin per analyte, transforms the data and bounds to
#' latent-Gaussian space, runs the tobit Gibbs sampler there, and back-transforms.
#' When `margin_draw = TRUE` each call draws the margin parameters from their
#' asymptotic posterior ([draw_margin()]) rather than reusing the MLE, so that
#' repeated calls (Rubin's-rules pooling) carry margin- and tail-estimation
#' uncertainty into the between-imputation variance; such a call is one stochastic
#' draw (seed for reproducibility). With `margin_draw = FALSE` it plugs in the MLE
#' margin -- the near-unbiased *point* imputation. The recommended MI workflow uses
#' the plug-in for the point estimate and drawn-margin imputations for the
#' between-imputation variance, which keeps the coverage fix without the small
#' Jensen (nonlinear back-transform) bias the drawn point picks up at heavy
#' censoring x strong skew.
#' @keywords internal
#' @noRd
copula_impute <- function(x, iters_all = 10, iters_each = 1,
                          initial = "bounds", verbose = FALSE,
                          margin_draw = TRUE) {
  data <- x$data_wide
  lo <- x$lo_mat
  hi <- x$hi_mat
  miss <- x$to_impute
  p <- ncol(data)

  margins <- vector("list", p)
  zobj <- x                                   # a cens_bounds copy in latent space
  for (j in seq_len(p)) {
    obsj <- !miss[, j]
    fit <- fit_shash_margin(data[obsj, j], lo[miss[, j], j], hi[miss[, j], j])
    m <- if (margin_draw) draw_margin(fit)    # proper-MI draw ...
         else list(mu = fit$mu, sigma = fit$sigma, eps = fit$eps)  # ... or MLE plug-in
    margins[[j]] <- m                         # reuse it for the back-transform
    # Transform observed values and the (monotone) censoring bounds. NA at
    # observed-cell bounds and +/-Inf endpoints all map through unchanged.
    zobj$data_wide[obsj, j] <- x_to_z(data[obsj, j], m$mu, m$sigma, m$eps)
    zobj$lo_mat[, j] <- x_to_z(lo[, j], m$mu, m$sigma, m$eps)
    zobj$hi_mat[, j] <- x_to_z(hi[, j], m$mu, m$sigma, m$eps)
  }

  # Impute in latent space with the (now correctly specified) tobit conditional.
  zfilled <- gsimp_impute(zobj, iters_all = iters_all, iters_each = iters_each,
                          initial = initial, imp_model = "tobit",
                          verbose = verbose)

  out <- zfilled
  for (j in seq_len(p)) {
    m <- margins[[j]]
    out[, j] <- z_to_x(zfilled[, j], m$mu, m$sigma, m$eps)   # back to working scale
  }
  out
}
