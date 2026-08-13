test_that("x_to_z and z_to_x are exact inverses", {
  x <- c(-2, 0, 1.5, 4, 10)
  for (eps in c(-0.6, 0, 0.6)) {
    z <- x_to_z(x, mu = 1, sigma = 2, eps = eps)
    expect_equal(z_to_x(z, mu = 1, sigma = 2, eps = eps), x, tolerance = 1e-10)
  }
})

test_that("the transforms are monotone increasing and preserve infinite ends", {
  z <- x_to_z(c(-Inf, 1, 2, 3, Inf), mu = 0.5, sigma = 1.3, eps = 0.4)
  expect_identical(z[1], -Inf)
  expect_identical(z[5], Inf)
  expect_true(all(diff(z[2:4]) > 0))
})

test_that("fit_shash_margin recovers skewness and degrades to Gaussian", {
  set.seed(1)
  # right-skewed data (eps > 0)
  xr <- 2 + sinh(asinh(rnorm(2000)) + 0.7)
  fr <- fit_shash_margin(xr, lower = numeric(0), upper = numeric(0))
  expect_gt(fr$eps, 0.3)

  # symmetric data -> eps near 0 (graceful degradation to tobit)
  xs <- rnorm(2000, mean = 5, sd = 2)
  fs <- fit_shash_margin(xs, lower = numeric(0), upper = numeric(0))
  expect_lt(abs(fs$eps), 0.15)
  expect_equal(fs$mu, 5, tolerance = 0.2)
  expect_equal(fs$sigma, 2, tolerance = 0.2)
})

test_that("fit_shash_margin uses interval-censored contributions", {
  set.seed(2)
  x <- 1 + sinh(asinh(rnorm(1000)) + 0.5)
  lod <- stats::quantile(x, 0.3)
  obs <- x[x >= lod]
  ncens <- sum(x < lod)
  fit <- fit_shash_margin(obs, lower = rep(-Inf, ncens), upper = rep(lod, ncens))
  expect_true(all(is.finite(c(fit$mu, fit$sigma, fit$eps))))
  expect_gt(fit$sigma, 0)
})

test_that("draw_margin returns a valid parameter set", {
  set.seed(3)
  fit <- fit_shash_margin(2 + sinh(asinh(rnorm(500)) + 0.6),
                          lower = numeric(0), upper = numeric(0))
  d <- draw_margin(fit)
  expect_named(d, c("mu", "sigma", "eps"))
  expect_gt(d$sigma, 0)

  # No covariance -> plug-in (deterministic) return.
  plug <- draw_margin(list(mu = 1, sigma = 2, eps = 0.3, par = c(1, log(2), 0.3), V = NULL))
  expect_equal(plug, list(mu = 1, sigma = 2, eps = 0.3))
})
