# Tests for as_cens_list(): wide data frame -> named list of cens_data.

test_that("left_right format builds one cens_data per analyte", {
  d <- data.frame(ID = 1:3,
                  A.left = c(0, 2, 5), A.right = c(1, 4, 5),
                  B.left = c(3, 0, 7), B.right = c(3, 1, 7))
  cl <- as_cens_list(d)
  expect_s3_class(cl, "cens_list")
  expect_true(is_cens_list(cl))
  expect_equal(names(cl), c("A", "B"))
  expect_true(all(vapply(cl, is_cens_data, logical(1))))
  expect_equal(as.character(cl$A$category),
               c("non_detect", "detected_not_quantified", "quantified"))
})

test_that("flagged format maps qualifiers to the two tiers", {
  f <- data.frame(ID = 1:4,
                  A = c(2.1, 0.5, 3.4, 0.5), A.flag = c("", "<", "", "ND"))
  cl <- as_cens_list(f, format = "flagged")
  expect_equal(names(cl), "A")
  expect_equal(as.character(cl$A$category),
               c("quantified", "non_detect", "quantified", "non_detect"),
               ignore_attr = TRUE)
})

test_that("explicit analytes and custom suffixes/flags are honoured", {
  d <- data.frame(Pb_lo = c(0, 2), Pb_hi = c(1, 2))
  cl <- as_cens_list(d, analytes = "Pb", left_suffix = "_lo",
                     right_suffix = "_hi")
  expect_equal(names(cl), "Pb")

  f <- data.frame(Hg = c(1, 5), Hg.q = c("U", ""))
  cl2 <- as_cens_list(f, format = "flagged", flag_suffix = ".q",
                      nd_flags = "U")
  expect_equal(as.character(cl2$Hg$category), c("non_detect", "quantified"))
})

test_that("as_cens_list output feeds build_bounds and desc_np directly", {
  d <- data.frame(A.left = c(0, 2, 5), A.right = c(1, 4, 5),
                  B.left = c(3, 0, 7), B.right = c(3, 1, 7))
  cl <- as_cens_list(d)
  expect_s3_class(build_bounds(cl), "cens_bounds")
  expect_s3_class(desc_np(cl), "cens_np_fits")
})

test_that("as_cens_list errors informatively", {
  expect_error(as_cens_list(data.frame(x = 1, y = 2)), "No analytes found")
  expect_error(
    as_cens_list(data.frame(A.left = 1), analytes = "A"),
    "Missing column"
  )
})

test_that("print.cens_list summarises analytes", {
  d <- data.frame(A.left = c(0, 2), A.right = c(1, 2))
  expect_output(print(as_cens_list(d)), "cens_list")
  expect_output(print(as_cens_list(d)), "A")
})
