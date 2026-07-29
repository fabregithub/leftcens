# Tests for the shared interval-censored data representation (core-cens.R).

test_that("as_interval_data classifies the three tiers correctly", {
  x <- as_interval_data(
    left  = c(0, 1, 5),
    right = c(1, 3, 5),
    log_transform = FALSE
  )
  expect_s3_class(x, "cens_data")
  expect_true(is_cens_data(x))
  expect_equal(
    as.character(x$category),
    c("non_detect", "detected_not_quantified", "quantified")
  )
  expect_equal(levels(x$category), cens_categories)
})

test_that("left == right is quantified even away from zero", {
  x <- as_interval_data(left = c(0, 4), right = c(0, 4), log_transform = FALSE)
  expect_equal(as.character(x$category), c("quantified", "quantified"))
})

test_that("left == 0 is a non-detect, not an error, under log_transform", {
  x <- as_interval_data(left = 0, right = 2, log_transform = TRUE)
  expect_equal(as.character(x$category), "non_detect")
  expect_identical(x$left, -Inf)        # log(0)
  expect_equal(x$right, log(2))
  expect_true(x$log_transform)
})

test_that("bounds are stored on the original scale when log_transform = FALSE", {
  x <- as_interval_data(left = 2, right = 8, log_transform = FALSE)
  expect_equal(x$left, 2)
  expect_equal(x$right, 8)
  expect_false(x$log_transform)
})

test_that("left > right is an error and names the offending index", {
  expect_error(
    as_interval_data(left = c(1, 9), right = c(2, 3), log_transform = FALSE),
    "index.*2"
  )
})

test_that("negative bounds error only under log_transform", {
  expect_error(
    as_interval_data(left = -1, right = 2, log_transform = TRUE),
    "non-negative"
  )
  expect_silent(as_interval_data(left = -1, right = 2, log_transform = FALSE))
})

test_that("length-1 recycling works and mismatched lengths error", {
  x <- as_interval_data(left = 0, right = c(1, 2, 3), log_transform = FALSE)
  expect_length(x$left, 3)
  expect_equal(as.character(unique(x$category)), "non_detect")
  expect_error(
    as_interval_data(left = c(0, 1), right = c(1, 2, 3), log_transform = FALSE),
    "length 1"
  )
})

test_that("zero-length input yields an empty cens_data", {
  x <- as_interval_data(left = numeric(0), right = numeric(0))
  expect_s3_class(x, "cens_data")
  expect_length(x$left, 0)
})

test_that("input validation rejects non-numeric and bad log_transform", {
  expect_error(as_interval_data("a", 1), "numeric")
  expect_error(as_interval_data(1, 2, log_transform = NA), "single")
  expect_error(as_interval_data(1, 2, log_transform = c(TRUE, FALSE)), "single")
})

test_that("as_interval_data_from_single collapses to the two-tier case", {
  x <- as_interval_data_from_single(
    value = c(2.0, 5.1),
    censored = c(TRUE, FALSE),
    log_transform = FALSE
  )
  expect_equal(as.character(x$category), c("non_detect", "quantified"))
  # censored -> (0, value); detected -> (value, value)
  expect_equal(x$left, c(0, 5.1))
  expect_equal(x$right, c(2.0, 5.1))
})

test_that("as_interval_data_from_single validates its inputs", {
  expect_error(as_interval_data_from_single("a", TRUE), "numeric")
  expect_error(as_interval_data_from_single(1, "yes"), "logical")
})

test_that("validate_cens_data catches a corrupted object", {
  x <- as_interval_data(left = 1, right = 2, log_transform = FALSE)
  x$right <- 0  # now left > right
  expect_error(validate_cens_data(x), "left.*<=.*right")
})

test_that("print returns its input invisibly", {
  x <- as_interval_data(left = c(0, 1, 5), right = c(1, 3, 5))
  expect_output(print(x), "cens_data")
  expect_invisible(print(x))
})
