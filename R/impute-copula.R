# Gaussian-copula ("normal-scores") imputation for interval-censored data.
#
# S1 in validation/PHASE3_SKEW_PLAN.md. Motivation: the tobit conditional
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

# The sinh-arcsinh margin fit (`fit_shash_margin`), its proper-MI parameter draw
# (`draw_margin`), and the latent-scale transforms (`x_to_z` / `z_to_x`) now live
# in R/shash-margin.R (exported API). They are used below unchanged.

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
