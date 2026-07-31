# Tests for multi-analyte descriptive collections (cens_np_fits / cens_sp_fits).

make_list <- function() {
  as_cens_list(data.frame(
    A.left = c(0, 0, 2, 3, 5, 8), A.right = c(1, 1, 4, 3, 5, 8),
    B.left = c(3, 0, 0, 7, 2, 6), B.right = c(3, 1, 1, 7, 4, 6)
  ))
}

test_that("desc_np on a single cens_data is unchanged (backward compatible)", {
  fit <- desc_np(as_interval_data(c(0, 1, 5, 8), c(1, 3, 5, 8)))
  expect_s3_class(fit, "cens_np_fit")
  expect_false(inherits(fit, "cens_np_fits"))
})

test_that("desc_np on a list returns a cens_np_fits collection", {
  fits <- desc_np(make_list())
  expect_s3_class(fits, "cens_np_fits")
  expect_equal(names(fits), c("A", "B"))
  expect_true(all(vapply(fits, inherits, logical(1), "cens_np_fit")))
})

test_that("quantile.cens_np_fits returns an analyte-by-probability matrix", {
  q <- quantile(desc_np(make_list()), probs = c(0.5, 0.9))
  expect_true(is.matrix(q))
  expect_equal(rownames(q), c("A", "B"))
  expect_equal(colnames(q), c("50%", "90%"))
})

test_that("as.data.frame / tidy / glance stack analytes with an analyte column", {
  fits <- desc_np(make_list())
  df <- as.data.frame(fits, probs = c(0.5, 0.9))
  expect_true("analyte" %in% names(df))
  expect_equal(sort(unique(df$analyte)), c("A", "B"))
  expect_equal(nrow(df), 4)                       # 2 analytes x 2 probs

  td <- tidy(fits, probs = c(0.5, 0.9))
  expect_true("analyte" %in% names(td))
  gl <- glance(fits)
  expect_equal(nrow(gl), 2)
  expect_true(all(c("analyte", "quantitation_limit") %in% names(gl)))
})

test_that("print and plot work for the np collection", {
  fits <- desc_np(make_list())
  expect_output(print(fits), "cens_np_fits")
  tmp <- tempfile(fileext = ".pdf"); pdf(tmp)
  on.exit({ dev.off(); unlink(tmp) }, add = TRUE)
  expect_invisible(plot(fits))
})

test_that("desc_sp on a list returns a cens_sp_fits collection with combined table", {
  cl <- make_list()
  grp <- factor(rep(c("a", "b"), 3))
  fits <- desc_sp(cl, covariates = grp, model = "ph")
  expect_s3_class(fits, "cens_sp_fits")
  expect_output(print(fits), "cens_sp_fits")

  df <- as.data.frame(fits)
  expect_true(all(c("analyte", "term", "coef", "hazard_ratio") %in% names(df)))
  expect_equal(sort(unique(df$analyte)), c("A", "B"))

  gl <- glance(fits)
  expect_equal(nrow(gl), 2)
  expect_true(all(c("analyte", "model", "nobs") %in% names(gl)))
})
