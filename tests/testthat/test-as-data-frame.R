# Tests for the as.data.frame() methods that make descriptive summaries
# ready for table-formatting packages.

test_that("as.data.frame(cens_np_fit) returns a tidy quantile table", {
  fit <- desc_np(as_interval_data(c(0, 1, 5, 8), c(1, 3, 5, 8),
                                  log_transform = FALSE))
  df <- as.data.frame(fit, probs = c(0.25, 0.5, 0.9))
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 3)
  expect_equal(names(df), c("probability", "quantile", "estimate"))
  expect_equal(df$probability, c(0.25, 0.5, 0.9))
  expect_equal(df$quantile, c("25%", "50%", "90%"))
  expect_type(df$estimate, "double")
  expect_equal(attr(df, "quantitation_limit"), fit$quantitation_limit)
})

test_that("as.data.frame(cens_np_fit) honours the ql masking", {
  x <- as_interval_data(c(0, 0, 2, 3, 5, 5, 8), c(1, 2, 4, 3, 7, 5, 8),
                        log_transform = FALSE)                     # QL = 7
  fit <- desc_np(x)
  masked <- as.data.frame(fit, probs = c(0.1, 0.9))
  raw <- as.data.frame(fit, probs = c(0.1, 0.9), ql = 0)
  expect_true(all(is.na(masked$estimate[raw$estimate < 7])))
  expect_false(anyNA(raw$estimate))
})

test_that("as.data.frame(cens_np_fit) can be formatted by a table package", {
  fit <- desc_np(as_interval_data(c(0, 1, 5, 8), c(1, 3, 5, 8),
                                  log_transform = FALSE))
  # knitr is a Suggest and the canonical minimal table formatter
  skip_if_not_installed("knitr")
  tab <- knitr::kable(as.data.frame(fit))
  expect_true(any(grepl("estimate", tab)))
})

test_that("as.data.frame(cens_sp_fit) returns a coefficient table", {
  x <- as_interval_data(c(0, 1, 5, 8, 0, 2), c(1, 3, 5, 8, 1, 4),
                        log_transform = FALSE)
  grp <- factor(c("a", "a", "a", "b", "b", "b"))
  df <- as.data.frame(desc_sp(x, covariates = grp, model = "ph"))
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 1)
  expect_true(all(c("term", "coef", "hazard_ratio") %in% names(df)))
  expect_equal(df$hazard_ratio, exp(df$coef))
})

test_that("as.data.frame(cens_sp_fit) adds SE / p when bootstrapped", {
  skip_on_cran()
  x <- as_interval_data(c(0, 1, 5, 8, 0, 2, 6, 1), c(1, 3, 5, 8, 1, 4, 6, 3),
                        log_transform = FALSE)
  grp <- factor(c("a", "a", "a", "b", "b", "b", "a", "b"))
  df <- as.data.frame(desc_sp(x, covariates = grp, model = "po",
                              bs_samples = 20))
  expect_true(all(c("std_error", "p_value", "odds_ratio") %in% names(df)))
  expect_true(all(df$std_error > 0))
  expect_true(all(df$p_value >= 0 & df$p_value <= 1))
})

test_that("as.data.frame(cens_sp_fit) is empty for an intercept-only fit", {
  x <- as_interval_data(c(0, 1, 5, 8), c(1, 3, 5, 8), log_transform = FALSE)
  df <- as.data.frame(desc_sp(x, covariates = NULL))
  expect_equal(nrow(df), 0)
})
