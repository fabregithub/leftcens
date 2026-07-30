# Tests for the pre-flight reliability tool (preflight.R).

test_that("preflight_reliability returns the expected structure", {
  pf <- preflight_reliability(n = 60, nd_frac = c(0.2, 0.4, 0.55),
                              n_rep = 3, M = 3, iters_all = 5,
                              imp_model = "ridge", seed = 1)
  expect_s3_class(pf, "preflight")
  expect_equal(nrow(pf$by_analyte), 3)
  expect_true(all(c("analyte", "nd_frac", "mean_bias", "median_bias",
                    "rmse", "coverage", "reliable") %in% names(pf$by_analyte)))
  expect_equal(pf$by_analyte$nd_frac, c(0.2, 0.4, 0.55))
  expect_true(all(pf$by_analyte$coverage >= 0 & pf$by_analyte$coverage <= 1))
})

test_that("an analyte at >= 50% ND is never flagged reliable", {
  pf <- preflight_reliability(n = 60, nd_frac = c(0.25, 0.55, 0.70),
                              n_rep = 3, M = 3, iters_all = 5,
                              imp_model = "ridge", seed = 2)
  # analytes 2 and 3 (>=50% ND) must be unreliable regardless of the model
  expect_false(pf$by_analyte$reliable[2])
  expect_false(pf$by_analyte$reliable[3])
})

test_that("scalar nd_frac is recycled across p analytes", {
  pf <- preflight_reliability(n = 50, nd_frac = 0.3, p = 4,
                              n_rep = 2, M = 2, iters_all = 4,
                              imp_model = "ridge", seed = 3)
  expect_equal(nrow(pf$by_analyte), 4)
  expect_equal(pf$by_analyte$nd_frac, rep(0.3, 4))
})

test_that("results are reproducible with a seed", {
  a <- preflight_reliability(n = 50, nd_frac = c(0.3, 0.4), n_rep = 3, M = 3,
                             iters_all = 5, imp_model = "ridge", seed = 42)
  b <- preflight_reliability(n = 50, nd_frac = c(0.3, 0.4), n_rep = 3, M = 3,
                             iters_all = 5, imp_model = "ridge", seed = 42)
  expect_equal(a$by_analyte, b$by_analyte)
})

test_that("input validation errors", {
  expect_error(preflight_reliability(60, nd_frac = 1.2), "\\[0, 1\\)")
  expect_error(preflight_reliability(60, nd_frac = 0.3, rho = 1.5), "rho")
  expect_error(preflight_reliability(60, nd_frac = 0.6, dnq_frac = 0.5),
               "\\[0, 1\\)")
  expect_error(preflight_reliability(60, nd_frac = 0.3, p = 2, n_rep = 0),
               "positive")
})

test_that("preflight_from_data derives the design from cens_data", {
  set.seed(10)
  v1 <- exp(stats::rnorm(80)); v2 <- exp(stats::rnorm(80))
  lim1 <- stats::quantile(v1, 0.3); lim2 <- stats::quantile(v2, 0.45)
  a <- as_interval_data(ifelse(v1 < lim1, 0, v1), ifelse(v1 < lim1, lim1, v1))
  b <- as_interval_data(ifelse(v2 < lim2, 0, v2), ifelse(v2 < lim2, lim2, v2))

  pf <- preflight_from_data(list(a = a, b = b), n_rep = 3, M = 3,
                            iters_all = 5, imp_model = "ridge", seed = 5)
  expect_s3_class(pf, "preflight")
  expect_equal(nrow(pf$by_analyte), 2)
  # derived non-detect fractions should be close to the ~0.30 / ~0.45 imposed
  expect_equal(pf$by_analyte$nd_frac, c(0.30, 0.45), tolerance = 0.05)
  expect_equal(pf$config$n, 80)
})

test_that("print reports a verdict", {
  pf <- preflight_reliability(n = 50, nd_frac = c(0.3, 0.6), n_rep = 2, M = 2,
                              iters_all = 4, imp_model = "ridge", seed = 7)
  expect_output(print(pf), "Verdict")
  expect_output(print(pf), "preflight")
  expect_invisible(print(pf))
})

test_that("with tobit, a low-ND analyte shows small bias", {
  skip_on_cran()
  pf <- preflight_reliability(n = 120, nd_frac = 0.30, p = 3, rho = 0.5,
                              n_rep = 5, M = 5, iters_all = 10,
                              imp_model = "tobit", seed = 11)
  expect_lt(max(abs(pf$by_analyte$mean_bias)), 0.15)
  expect_true(all(abs(pf$by_analyte$median_bias) < 0.1))
})
