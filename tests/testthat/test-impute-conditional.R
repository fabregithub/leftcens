make_censored <- function(n = 400, nd = 0.35, skew = 0.6, seed = 1) {
  set.seed(seed)
  z <- stats::rnorm(n)
  x1 <- sinh(asinh(0.6 * z + sqrt(1 - 0.36) * stats::rnorm(n)) + skew)  # skewed, z-correlated
  lod <- as.numeric(stats::quantile(x1, nd))
  cens <- x1 < lod
  y <- x1; y[cens] <- NA
  list(y = y, x = data.frame(z = z), x1 = x1, cens = cens, lod = lod,
       lower = ifelse(cens, -Inf, NA), upper = ifelse(cens, lod, NA))
}

test_that("returns an n x m matrix with observed cells carried through", {
  d <- make_censored()
  imp <- impute_censored_conditional(d$y, d$x, d$lower, d$upper, m = 4)
  expect_equal(dim(imp), c(length(d$y), 4L))
  # observed rows identical across all imputations
  obs <- !d$cens
  expect_equal(imp[obs, 1], d$y[obs])
  expect_equal(imp[obs, 3], d$y[obs])
})

test_that("imputed non-detects respect the LOD (left-censoring bound)", {
  d <- make_censored()
  imp <- impute_censored_conditional(d$y, d$x, d$lower, d$upper, m = 5)
  expect_true(all(imp[d$cens, ] <= d$lod + 1e-8))
})

test_that("interval censoring keeps draws inside the interval", {
  set.seed(4)
  n <- 300; z <- rnorm(n); x1 <- 1 + z + rnorm(n)
  lo <- 0.5; hi <- 1.5
  cens <- x1 > lo & x1 < hi
  y <- x1; y[cens] <- NA
  lower <- ifelse(cens, lo, NA); upper <- ifelse(cens, hi, NA)
  imp <- impute_censored_conditional(y, data.frame(z), lower, upper, m = 3,
                                     margin = "gaussian")
  v <- imp[cens, ]
  expect_true(all(v > lo - 1e-8 & v < hi + 1e-8))
})

test_that("conditioning on Y recovers the exposure-response slope under skew", {
  set.seed(7)
  n <- 1200; z <- rnorm(n)
  x1 <- sinh(asinh(0.5 * z + sqrt(1 - 0.25) * rnorm(n)) + 0.7)   # right-skewed
  y <- 0.5 * x1 + 0.4 * z + rnorm(n)                             # true slope 0.5
  lod <- as.numeric(quantile(x1, 0.4)); cens <- x1 < lod
  yc <- x1; yc[cens] <- NA
  lower <- ifelse(cens, -Inf, NA); upper <- ifelse(cens, lod, NA)
  imp <- impute_censored_conditional(yc, data.frame(y = y, z = z), lower, upper, m = 20)
  b <- vapply(seq_len(20), function(i)
    coef(lm(y ~ imp[, i] + z))[2], numeric(1))
  expect_equal(mean(b), 0.5, tolerance = 0.05)                   # near-unbiased
})

test_that("no censored cells returns replicated y; input checks fire", {
  y <- rnorm(50)
  imp <- impute_censored_conditional(y, data.frame(z = rnorm(50)),
                                     rep(NA, 50), rep(NA, 50), m = 3)
  expect_equal(imp[, 2], y)

  expect_error(impute_censored_conditional(y, data.frame(z = rnorm(49)),
                                           rep(NA, 50), rep(NA, 50)),
               "one row per element")
  expect_error(impute_censored_conditional(y, data.frame(z = rnorm(50)),
                                           rep(NA, 49), rep(NA, 50)),
               "same length")
})
