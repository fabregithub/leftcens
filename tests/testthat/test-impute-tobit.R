# Tests for the censored/Tobit conditional model (imp_model = "tobit").

within_bounds <- function(filled, bnds) {
  m <- bnds$to_impute
  all(filled[m] <= bnds$hi_mat[m] + 1e-8) &&
    all(is.infinite(bnds$lo_mat[m]) | filled[m] >= bnds$lo_mat[m] - 1e-8)
}

# Small correlated three-tier dataset builder (pure left-censoring).
sim_bounds <- function(n = 300, p = 6, rho = 0.5, nd = 0.35, seed = 1) {
  set.seed(seed)
  R <- matrix(rho, p, p); diag(R) <- 1
  z <- sweep(matrix(stats::rnorm(n * p), n, p) %*% chol(R), 2,
             seq(0, 1.5, length.out = p), "+")
  conc <- exp(z)
  cds <- lapply(seq_len(p), function(j) {
    lim <- stats::quantile(conc[, j], nd)
    as_interval_data(ifelse(conc[, j] < lim, 0, conc[, j]),
                     ifelse(conc[, j] < lim, lim, conc[, j]))
  })
  list(bounds = build_bounds(cds), z = z)
}

test_that("tobit runs and imputes within bounds", {
  d <- sim_bounds(n = 150, p = 5, nd = 0.3, seed = 2)
  filled <- gsimp_impute(d$bounds, iters_all = 10, imp_model = "tobit")
  expect_false(anyNA(filled))
  expect_true(within_bounds(filled, d$bounds))
})

test_that("tobit sharply reduces mean bias vs ridge at moderate censoring", {
  d <- sim_bounds(n = 300, p = 6, rho = 0.5, nd = 0.35, seed = 3)
  bias <- function(f) mean(colMeans(f) - colMeans(d$z))
  b_ridge <- bias(gsimp_impute(d$bounds, iters_all = 20, imp_model = "ridge"))
  b_tobit <- bias(gsimp_impute(d$bounds, iters_all = 20, imp_model = "tobit"))
  # ridge is upward-biased (~0.13 here); tobit should be near zero and clearly
  # smaller in magnitude.
  expect_gt(b_ridge, 0.05)
  expect_lt(abs(b_tobit), 0.05)
  expect_lt(abs(b_tobit), abs(b_ridge))
})

test_that("tobit falls back to ridge on wide (p >= n) designs without error", {
  set.seed(4)
  n <- 6L; p <- 10L
  cds <- lapply(seq_len(p), function(j) {
    v <- exp(stats::rnorm(n)); lim <- stats::quantile(v, 0.3)
    as_interval_data(ifelse(v < lim, 0, v), ifelse(v < lim, lim, v))
  })
  bnds <- build_bounds(cds)
  filled <- gsimp_impute(bnds, iters_all = 5, imp_model = "tobit")
  expect_false(anyNA(filled))
  expect_true(within_bounds(filled, bnds))
})

test_that("unknown imp_model still errors informatively", {
  d <- sim_bounds(n = 50, p = 3, seed = 5)
  expect_error(gsimp_impute(d$bounds, imp_model = "nope"), "tobit")
})
