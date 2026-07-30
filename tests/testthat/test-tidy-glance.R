# Tests for the broom-style tidy() / glance() methods.

test_that("tidy(cens_np_fit) returns a per-quantile data frame", {
  fit <- desc_np(as_interval_data(c(0, 1, 5, 8), c(1, 3, 5, 8),
                                  log_transform = FALSE))
  td <- tidy(fit, probs = c(0.25, 0.5, 0.9))
  expect_s3_class(td, "data.frame")
  expect_equal(nrow(td), 3)
  expect_equal(names(td), c("term", "probability", "estimate"))
  expect_equal(td$term, c("25%", "50%", "90%"))
})

test_that("tidy(cens_np_fit) honours ql masking", {
  x <- as_interval_data(c(0, 0, 2, 3, 5, 5, 8), c(1, 2, 4, 3, 7, 5, 8),
                        log_transform = FALSE)                    # QL = 7
  fit <- desc_np(x)
  td <- tidy(fit, probs = c(0.1, 0.9))
  raw <- tidy(fit, probs = c(0.1, 0.9), ql = 0)
  expect_true(all(is.na(td$estimate[raw$estimate < 7])))
})

test_that("glance(cens_np_fit) returns a one-row model summary", {
  fit <- desc_np(as_interval_data(c(0, 0, 2, 5, 8), c(1, 2, 4, 5, 8),
                                  log_transform = FALSE))
  gl <- glance(fit)
  expect_equal(nrow(gl), 1)
  expect_true(all(c("method", "nobs", "n_non_detect", "detection_limit",
                    "quantitation_limit") %in% names(gl)))
  expect_equal(gl$nobs, 5)
})

test_that("tidy(cens_sp_fit) returns a coefficient table with exponentiate", {
  x <- as_interval_data(c(0, 1, 5, 8, 0, 2), c(1, 3, 5, 8, 1, 4),
                        log_transform = FALSE)
  grp <- factor(c("a", "a", "a", "b", "b", "b"))
  sp <- desc_sp(x, covariates = grp, model = "ph")
  td <- tidy(sp)
  expect_equal(names(td), c("term", "estimate"))
  te <- tidy(sp, exponentiate = TRUE)
  expect_equal(te$estimate, exp(td$estimate))
})

test_that("tidy(cens_sp_fit) adds inference columns when bootstrapped", {
  skip_on_cran()
  x <- as_interval_data(c(0, 1, 5, 8, 0, 2, 6, 1), c(1, 3, 5, 8, 1, 4, 6, 3),
                        log_transform = FALSE)
  grp <- factor(c("a", "a", "a", "b", "b", "b", "a", "b"))
  td <- tidy(desc_sp(x, covariates = grp, model = "ph", bs_samples = 20))
  expect_true(all(c("std.error", "statistic", "p.value") %in% names(td)))
  expect_true(all(td$p.value >= 0 & td$p.value <= 1))
})

test_that("glance(cens_sp_fit) returns a one-row summary", {
  x <- as_interval_data(c(0, 1, 5, 8, 0, 2), c(1, 3, 5, 8, 1, 4),
                        log_transform = FALSE)
  grp <- factor(c("a", "a", "a", "b", "b", "b"))
  gl <- glance(desc_sp(x, covariates = grp, model = "po"))
  expect_equal(nrow(gl), 1)
  expect_true(all(c("model", "nobs", "npar", "logLik") %in% names(gl)))
  expect_equal(gl$model, "po")
  expect_equal(gl$npar, 1)
})
