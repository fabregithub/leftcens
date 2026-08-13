# =============================================================================
# Sinh-arcsinh (SHASH) transform-to-normality margin.
# -----------------------------------------------------------------------------
# A 3-parameter skewed margin fit by interval-censored MLE, plus the monotone
# transforms between the data scale and a latent standard-normal scale:
#
#     X = mu + sigma * sinh( asinh(Z) + eps ),   Z ~ N(0, 1)
#
# (location mu, scale sigma, skewness eps; eps = 0 is a plain Gaussian). This is
# the skew-handling core of the copula imputer (see `gsimp_impute(imp_model =
# "copula")`); it is exported so that general-predictor imputation
# (`impute_censored_conditional()`) and downstream code can map a skewed,
# left-/interval-censored variable to a latent scale where a linear-Gaussian
# conditional model is correctly specified, then map back.
# =============================================================================

#' Sinh-arcsinh transforms between the data scale and a latent normal scale
#'
#' `x_to_z()` maps data to its latent standard-normal score and `z_to_x()` is the
#' exact inverse (the marginal quantile map \eqn{F^{-1}(\Phi(\cdot))}). Both are
#' monotone increasing, so censoring bounds transform through unchanged (and
#' `-Inf` / `Inf` endpoints map to themselves).
#'
#' @param x,z Numeric vectors on the data scale (`x`) or latent normal scale (`z`).
#' @param mu,sigma,eps Margin parameters: location, scale (> 0), and skewness
#'   (`eps = 0` gives a plain Gaussian margin). Typically from [fit_shash_margin()]
#'   or a draw from [draw_margin()].
#' @return A numeric vector the same length as the input.
#' @seealso [fit_shash_margin()], [draw_margin()], [impute_censored_conditional()]
#' @examples
#' z <- x_to_z(c(1, 2, 5), mu = 2, sigma = 1, eps = 0.5)
#' z_to_x(z, mu = 2, sigma = 1, eps = 0.5)   # recovers c(1, 2, 5)
#' @export
x_to_z <- function(x, mu, sigma, eps) sinh(asinh((x - mu) / sigma) - eps)

#' @rdname x_to_z
#' @export
z_to_x <- function(z, mu, sigma, eps) mu + sigma * sinh(asinh(z) + eps)

#' Fit a sinh-arcsinh margin by interval-censored maximum likelihood
#'
#' Estimates the location, scale, and skewness of a 3-parameter sinh-arcsinh
#' margin from a mix of exactly observed values and interval-censored cells
#' (including left-censored non-detects). The fit is skew-aware, so non-detects
#' are later drawn from a correctly shaped tail rather than a Gaussian one.
#'
#' @param observed Numeric vector of exactly observed values.
#' @param lower,upper Lower/upper bounds of the censored cells (same length as
#'   each other). `lower` may be `-Inf` for left-censored non-detects and `upper`
#'   may be `Inf` for right-censored values. Pass length-0 vectors when there are
#'   no censored cells.
#' @return A list with the point estimate (`mu`, `sigma`, `eps`), the fitted
#'   parameter vector `par = c(mu, log sigma, eps)`, and its asymptotic covariance
#'   `V` (inverse observed information), or `V = NULL` when it cannot be formed
#'   (which makes [draw_margin()] fall back to the plug-in estimate). Falls back to
#'   a Gaussian `(mean, sd, 0)` when there is too little data or the optimiser
#'   fails, so it degrades gracefully to tobit.
#' @seealso [draw_margin()] for a proper-MI parameter draw; [x_to_z()] /
#'   [z_to_x()] for the transforms; [impute_censored_conditional()].
#' @examples
#' set.seed(1)
#' x <- 2 + sinh(asinh(rnorm(200)) + 0.6)      # right-skewed
#' lod <- quantile(x, 0.3)
#' obs <- x[x >= lod]
#' fit <- fit_shash_margin(obs, lower = rep(-Inf, sum(x < lod)),
#'                         upper = rep(lod, sum(x < lod)))
#' fit[c("mu", "sigma", "eps")]
#' @export
fit_shash_margin <- function(observed, lower, upper) {
  xo <- observed; lo_c <- lower; hi_c <- upper
  # Robust location/scale even with 0-1 observed values (an entirely censored
  # analyte): centre on the finite censoring bounds when there is nothing observed.
  mu0 <- if (length(xo) >= 1L) mean(xo) else NA_real_
  if (!is.finite(mu0)) {
    fin <- c(lo_c[is.finite(lo_c)], hi_c[is.finite(hi_c)])
    mu0 <- if (length(fin)) mean(fin) else 0
  }
  sd0 <- if (length(xo) >= 2L) stats::sd(xo) else NA_real_
  if (!is.finite(sd0) || sd0 <= 0) sd0 <- 1
  fallback <- list(mu = mu0, sigma = sd0, eps = 0,
                   par = c(mu0, log(sd0), 0), V = NULL)
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

#' Draw a sinh-arcsinh margin's parameters from their asymptotic posterior
#'
#' For *proper* multiple imputation each completed dataset should condition on a
#' drawn margin, not the shared MLE, so that the between-imputation variance
#' reflects margin- and (extrapolated) tail-estimation uncertainty. Returns the
#' plug-in estimate when the fit carries no usable covariance (`margin$V` is
#' `NULL`) or the draw is degenerate, so behaviour is unchanged where the fit is
#' untrustworthy.
#'
#' @param margin A fit from [fit_shash_margin()].
#' @return A list `(mu, sigma, eps)`: one draw of the margin parameters (or the
#'   plug-in estimate).
#' @seealso [fit_shash_margin()].
#' @examples
#' set.seed(1)
#' x <- 2 + sinh(asinh(rnorm(200)) + 0.6)
#' fit <- fit_shash_margin(x, lower = numeric(0), upper = numeric(0))
#' draw_margin(fit)
#' @export
draw_margin <- function(margin) {
  m <- margin
  plugin <- list(mu = m$mu, sigma = m$sigma, eps = m$eps)
  if (is.null(m$V)) return(plugin)
  d <- tryCatch(as.numeric(m$par + t(chol(m$V)) %*% stats::rnorm(3)),
                error = function(e) NULL)
  if (is.null(d) || !all(is.finite(d)) || abs(d[3]) > 10) return(plugin)
  s <- exp(d[2])
  if (!is.finite(s) || s <= 0) return(plugin)
  list(mu = d[1], sigma = s, eps = d[3])
}
