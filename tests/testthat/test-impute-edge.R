# Edge-case tests for gsimp_impute(): degenerate and stress configurations that
# the main tests don't exercise. All small and deterministic so they stay in the
# CI suite.

# Helper: every imputed cell lies inside its assigned interval.
expect_within_bounds <- function(filled, bnds) {
  m <- bnds$to_impute
  hi_ok <- all(filled[m] <= bnds$hi_mat[m] + 1e-8)
  lo_ok <- all(is.infinite(bnds$lo_mat[m]) | filled[m] >= bnds$lo_mat[m] - 1e-8)
  testthat::expect_true(hi_ok && lo_ok)
}

test_that("an entirely censored analyte still imputes within bounds", {
  set.seed(101)
  a <- as_interval_data(c(0, 2, 5, 8, 3, 1), c(1, 4, 5, 8, 3, 2)) # has observed
  b <- as_interval_data(rep(0, 6), rep(1, 6))                     # all non-detect
  bnds <- build_bounds(list(A = a, B = b))
  expect_true(all(bnds$to_impute[, "B"]))                          # sanity

  filled <- gsimp_impute(bnds, iters_all = 5)
  expect_false(anyNA(filled))
  expect_within_bounds(filled, bnds)
  # every cell of the all-censored analyte sits below its limit (log(1) = 0)
  expect_true(all(filled[, "B"] <= 0 + 1e-8))
})

test_that("more analytes than samples (p > n) stays solvable with ridge", {
  set.seed(102)
  n <- 6L; p <- 10L
  cols <- lapply(seq_len(p), function(j) {
    v <- exp(stats::rnorm(n))
    lim <- stats::quantile(v, 0.3)
    as_interval_data(ifelse(v < lim, 0, v), ifelse(v < lim, lim, v))
  })
  bnds <- build_bounds(cols)
  filled <- gsimp_impute(bnds, iters_all = 5, imp_model = "ridge")
  expect_equal(dim(filled), c(n, p))
  expect_false(anyNA(filled))
  expect_within_bounds(filled, bnds)
})

test_that("uncorrelated analytes impute sensibly (no signal, no blow-up)", {
  set.seed(103)
  n <- 80L; p <- 3L
  X <- matrix(stats::rnorm(n * p), n, p)          # independent columns
  conc <- exp(X)
  cols <- lapply(seq_len(p), function(j) {
    lim <- stats::quantile(conc[, j], 0.30)
    as_interval_data(ifelse(conc[, j] < lim, 0, conc[, j]),
                     ifelse(conc[, j] < lim, lim, conc[, j]))
  })
  bnds <- build_bounds(cols)
  filled <- gsimp_impute(bnds, iters_all = 10)
  expect_false(anyNA(filled))
  expect_true(all(is.finite(filled)))
  expect_within_bounds(filled, bnds)
})

test_that("extreme censoring (~85%) still runs", {
  set.seed(104)
  n <- 60L; p <- 3L
  X <- matrix(stats::rnorm(n * p), n, p) %*% chol(matrix(0.5, p, p) + diag(0.5, p))
  conc <- exp(X)
  cols <- lapply(seq_len(p), function(j) {
    lim <- stats::quantile(conc[, j], 0.85)       # 85% below the limit
    as_interval_data(ifelse(conc[, j] < lim, 0, conc[, j]),
                     ifelse(conc[, j] < lim, lim, conc[, j]))
  })
  bnds <- build_bounds(cols)
  expect_gt(mean(bnds$to_impute), 0.75)
  filled <- gsimp_impute(bnds, iters_all = 10)
  expect_false(anyNA(filled))
  expect_within_bounds(filled, bnds)
})

test_that("a single-column cens_bounds imputes with no predictors", {
  set.seed(105)
  v <- exp(stats::rnorm(40))
  lim <- stats::quantile(v, 0.35)
  cd <- as_interval_data(ifelse(v < lim, 0, v), ifelse(v < lim, lim, v))
  bnds <- build_bounds(cd)                         # p = 1, no conditioning cols
  filled <- gsimp_impute(bnds, iters_all = 10)
  expect_false(anyNA(filled))
  expect_within_bounds(filled, bnds)
})
