# Tests for the per-cell imputation bounds builder (impute-bounds.R).

test_that("build_bounds accepts a single cens_data as one column", {
  a <- as_interval_data(left = c(0, 2, 5), right = c(1, 4, 5),
                        log_transform = FALSE)
  b <- build_bounds(a)
  expect_s3_class(b, "cens_bounds")
  expect_true(is_cens_bounds(b))
  expect_equal(dim(b$lo_mat), c(3, 1))
  expect_equal(colnames(b$lo_mat), "V1")
})

test_that("list names become column names", {
  a <- as_interval_data(0, 1, log_transform = FALSE)
  b <- as_interval_data(3, 3, log_transform = FALSE)
  bnds <- build_bounds(list(A = a, B = b))
  expect_equal(colnames(bnds$data_wide), c("A", "B"))
  expect_equal(dim(bnds$data_wide), c(1, 2))
})

test_that("quantified cells are observed; censored cells are imputed", {
  # non_detect (0,1), DNQ (2,4), quantified (5,5)
  a <- as_interval_data(left = c(0, 2, 5), right = c(1, 4, 5),
                        log_transform = FALSE)
  b <- build_bounds(a)

  expect_equal(as.vector(b$to_impute), c(TRUE, TRUE, FALSE))

  # data_wide: value only at the quantified cell
  expect_equal(as.vector(b$data_wide), c(NA, NA, 5))
  # lo/hi are NA exactly where observed
  expect_equal(is.na(as.vector(b$lo_mat)), c(FALSE, FALSE, TRUE))
  expect_equal(is.na(as.vector(b$hi_mat)), c(FALSE, FALSE, TRUE))
})

test_that("non-detects become left-censored (-Inf, log MDL) on the log scale", {
  a <- as_interval_data(left = 0, right = 2, log_transform = TRUE) # non_detect
  b <- build_bounds(a)
  expect_identical(as.vector(b$lo_mat), -Inf)
  expect_equal(as.vector(b$hi_mat), log(2))
  expect_true(is.na(as.vector(b$data_wide)))  # censored -> not observed
})

test_that("detected-not-quantified cells get a finite interval", {
  a <- as_interval_data(left = 2, right = 4, log_transform = TRUE) # DNQ
  b <- build_bounds(a)
  expect_equal(as.vector(b$lo_mat), log(2))
  expect_equal(as.vector(b$hi_mat), log(4))
  expect_true(all(is.finite(c(b$lo_mat, b$hi_mat))))
})

test_that("data_wide missingness pattern exactly matches to_impute", {
  a <- as_interval_data(c(0, 2, 5, 0), c(1, 4, 5, 1), log_transform = TRUE)
  b <- build_bounds(a)
  expect_equal(is.na(b$data_wide), b$to_impute)
})

test_that("cells with NA bounds are imputed unbounded", {
  a <- as_interval_data(left = c(0, NA), right = c(1, NA),
                        log_transform = TRUE)
  b <- build_bounds(a)
  expect_equal(as.vector(b$to_impute), c(TRUE, TRUE))
  expect_identical(unname(b$lo_mat[2, 1]), -Inf)
  expect_identical(unname(b$hi_mat[2, 1]), Inf)
})

test_that("multi-column bounds line up cell by cell", {
  a <- as_interval_data(left = c(0, 2, 5), right = c(1, 4, 5),
                        log_transform = FALSE)  # ND, DNQ, Q
  b <- as_interval_data(left = c(3, 0, 7), right = c(3, 1, 7),
                        log_transform = FALSE)  # Q,  ND, Q
  bnds <- build_bounds(list(A = a, B = b))
  expect_equal(as.vector(bnds$to_impute[, "A"]), c(TRUE, TRUE, FALSE))
  expect_equal(as.vector(bnds$to_impute[, "B"]), c(FALSE, TRUE, FALSE))
  expect_equal(unname(bnds$data_wide[1, "B"]), 3)
  expect_equal(unname(bnds$data_wide[3, "A"]), 5)
})

test_that("unequal column lengths error", {
  a <- as_interval_data(c(0, 1), c(1, 2), log_transform = FALSE)
  b <- as_interval_data(0, 1, log_transform = FALSE)
  expect_error(build_bounds(list(a, b)), "same length")
})

test_that("mixed log_transform flags error", {
  a <- as_interval_data(0, 1, log_transform = TRUE)
  b <- as_interval_data(0, 1, log_transform = FALSE)
  expect_error(build_bounds(list(a, b)), "log_transform")
})

test_that("non-cens_data input errors", {
  expect_error(build_bounds(data.frame(a = 1)), "cens_data")
  expect_error(build_bounds(list(1, 2)), "cens_data")
})

test_that("print reports dimensions and imputation counts", {
  a <- as_interval_data(c(0, 2, 5), c(1, 4, 5), log_transform = TRUE)
  b <- build_bounds(a)
  expect_output(print(b), "cens_bounds")
  expect_output(print(b), "cells to impute: 2 of 3")
  expect_output(print(b), "left-censored")
  expect_invisible(print(b))
})
