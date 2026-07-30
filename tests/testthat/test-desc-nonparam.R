# Tests for the non-parametric NPMLE descriptive module (desc-nonparam.R).

# A small, hand-traceable interval-censored dataset on the original scale.
make_x <- function(log_transform = TRUE) {
  as_interval_data(
    left  = c(0, 0, 2, 3, 5, 5, 8),
    right = c(1, 2, 4, 3, 7, 5, 8),
    log_transform = log_transform
  )
}

test_that("desc_np returns a cens_np_fit with the expected metadata", {
  fit <- desc_np(make_x())
  expect_s3_class(fit, "cens_np_fit")
  expect_equal(fit$method, "turnbull")
  expect_equal(fit$backend, "survival")
  expect_equal(fit$n, 7)
  expect_equal(sum(fit$category), 7)
})

test_that("desc_np requires a cens_data object", {
  expect_error(desc_np(data.frame(a = 1)), "cens_data")
})

test_that("both methods run and produce ordered quantiles", {
  for (m in c("turnbull", "wang")) {
    fit <- desc_np(make_x(), method = m)
    q <- quantile(fit, probs = c(.1, .25, .5, .75, .9), ql = 0)  # raw estimates
    expect_length(q, 5)
    expect_false(any(is.na(q)))
    expect_true(all(diff(q) >= 0))          # monotone non-decreasing
    expect_true(all(q >= 0))                # concentrations are non-negative
  }
})

test_that("turnbull quantiles match a direct survfit call (wrapper cross-check)", {
  lo <- c(0, 0, 2, 3, 5, 5, 8)
  hi <- c(1, 2, 4, 3, 7, 5, 8)
  direct <- survival::survfit(survival::Surv(lo, hi, type = "interval2") ~ 1)
  probs <- c(.1, .25, .5, .75, .9)
  expected <- as.numeric(quantile(direct, probs = probs, conf.int = FALSE))

  fit <- desc_np(as_interval_data(lo, hi, log_transform = FALSE),
                 method = "turnbull")
  expect_equal(as.numeric(quantile(fit, probs = probs, ql = 0)), expected)
})

test_that("fit is invariant to the log_transform flag (monotone-invariance)", {
  probs <- c(.1, .25, .5, .75, .9)
  q_log <- quantile(desc_np(make_x(TRUE)), probs = probs, ql = 0)
  q_lin <- quantile(desc_np(make_x(FALSE)), probs = probs, ql = 0)
  expect_equal(as.numeric(q_log), as.numeric(q_lin), tolerance = 1e-8)
})

test_that("quantiles below the quantitation limit are NA by default", {
  # non-detect (0,1), DNQ (2,4) & (5,7), quantified 3, 5, 8 -> QL = 7
  x <- as_interval_data(left  = c(0, 0, 2, 3, 5, 5, 8),
                        right = c(1, 2, 4, 3, 7, 5, 8),
                        log_transform = FALSE)
  fit <- desc_np(x, method = "turnbull")
  expect_equal(fit$quantitation_limit, 7)
  expect_equal(fit$detection_limit, 2)      # highest non-detect upper bound

  q <- quantile(fit, probs = c(.1, .5, .9))
  raw <- quantile(fit, probs = c(.1, .5, .9), ql = 0)
  # any raw estimate below the QL becomes NA; those at/above it are kept
  expect_true(all(is.na(q[raw < 7])))
  expect_equal(q[raw >= 7], raw[raw >= 7])
})

test_that("a custom ql threshold (e.g. the detection limit) is honoured", {
  x <- as_interval_data(c(0, 0, 2, 3, 5, 8), c(1, 2, 4, 3, 5, 8),
                        log_transform = FALSE)
  fit <- desc_np(x, method = "turnbull")
  q_dl <- quantile(fit, probs = c(.25, .75), ql = fit$detection_limit)
  q_raw <- quantile(fit, probs = c(.25, .75), ql = 0)
  expect_true(all(is.na(q_dl[q_raw < fit$detection_limit])))
})

test_that("no censoring leaves quantitation limit NA and masks nothing", {
  x <- as_interval_data(c(1, 2, 3, 4), c(1, 2, 3, 4), log_transform = FALSE)
  fit <- desc_np(x, method = "turnbull")
  expect_true(is.na(fit$quantitation_limit))
  expect_false(anyNA(quantile(fit, probs = c(.25, .5, .75))))
})

test_that("missing bounds are dropped and counted", {
  x <- as_interval_data(left = c(0, 1, NA), right = c(1, 3, NA),
                        log_transform = FALSE)
  fit <- desc_np(x)
  expect_equal(fit$n, 2)
  expect_equal(fit$dropped, 1)
})

test_that("all-missing data errors informatively", {
  x <- as_interval_data(left = NA_real_, right = NA_real_, log_transform = FALSE)
  expect_error(desc_np(x), "non-missing")
})

test_that("quantile names carry the percentage", {
  q <- quantile(desc_np(make_x()), probs = c(.5, .9))
  expect_equal(names(q), c("50%", "90%"))
})

test_that("print and plot work without error", {
  fit <- desc_np(make_x())
  expect_output(print(fit), "cens_np_fit")
  expect_invisible(print(fit))
  tmp <- tempfile(fileext = ".pdf")
  pdf(tmp)
  on.exit({ dev.off(); unlink(tmp) }, add = TRUE)
  expect_invisible(plot(fit))
  expect_invisible(plot(desc_np(make_x(), method = "wang")))
})
