# Tests for the Gibbs-sampler imputation module (impute-gsimp.R).

# ---- rnorm_trunc ----------------------------------------------------------

test_that("rnorm_trunc respects finite interval bounds", {
  set.seed(1)
  x <- rnorm_trunc(2000, mean = 0, sd = 3, lower = -1, upper = 2)
  expect_true(all(x >= -1 & x <= 2))
  expect_gt(mean(x), -1)
  expect_lt(mean(x), 2)
})

test_that("rnorm_trunc handles left-censoring (-Inf lower) without underflow", {
  set.seed(2)
  x <- rnorm_trunc(2000, mean = 5, sd = 1, lower = -Inf, upper = log(2))
  expect_true(all(is.finite(x)))
  expect_true(all(x <= log(2)))
})

test_that("rnorm_trunc handles right-unbounded and point intervals", {
  set.seed(3)
  x <- rnorm_trunc(1000, mean = 0, sd = 1, lower = 1, upper = Inf)
  expect_true(all(x >= 1) && all(is.finite(x)))
  pt <- rnorm_trunc(5, mean = 9, sd = 1, lower = 4, upper = 4)
  expect_equal(pt, rep(4, 5))
})

test_that("rnorm_trunc recycles vectorised bounds and guards sd <= 0", {
  set.seed(4)
  x <- rnorm_trunc(4, mean = 0, sd = 1,
                   lower = c(-1, 0, -Inf, 2), upper = c(1, 5, 0, Inf))
  expect_true(x[1] >= -1 && x[1] <= 1)
  expect_true(x[2] >= 0 && x[2] <= 5)
  expect_true(x[3] <= 0)
  expect_true(x[4] >= 2)
  # sd = 0 collapses to the clamped mean
  expect_equal(rnorm_trunc(1, mean = 3, sd = 0, lower = 0, upper = 10), 3)
})

# ---- gsimp_impute basics --------------------------------------------------

make_bounds <- function(log_transform = TRUE) {
  a <- as_interval_data(c(0, 2, 5, 8, 3), c(1, 4, 5, 8, 3),
                        log_transform = log_transform)
  b <- as_interval_data(c(3, 0, 7, 2, 6), c(3, 1, 7, 4, 6),
                        log_transform = log_transform)
  build_bounds(list(A = a, B = b))
}

test_that("gsimp_impute fills all censored cells and preserves observed ones", {
  set.seed(10)
  bnds <- make_bounds()
  filled <- gsimp_impute(bnds, iters_all = 5)
  expect_equal(dim(filled), dim(bnds$data_wide))
  expect_false(anyNA(filled))
  # observed cells unchanged
  obs <- !bnds$to_impute
  expect_equal(filled[obs], bnds$data_wide[obs])
})

test_that("imputed values fall within their per-cell bounds", {
  set.seed(11)
  bnds <- make_bounds()
  filled <- gsimp_impute(bnds, iters_all = 8)
  m <- bnds$to_impute
  expect_true(all(filled[m] <= bnds$hi_mat[m] + 1e-8))
  lo_ok <- is.infinite(bnds$lo_mat[m]) | filled[m] >= bnds$lo_mat[m] - 1e-8
  expect_true(all(lo_ok))
})

test_that("no missing cells returns the data unchanged", {
  a <- as_interval_data(c(1, 2, 3), c(1, 2, 3), log_transform = TRUE) # all quantified
  bnds <- build_bounds(a)
  expect_equal(gsimp_impute(bnds), bnds$data_wide)
})

test_that("input and argument validation", {
  expect_error(gsimp_impute(list(a = 1)), "cens_bounds")
  bnds <- make_bounds()
  expect_error(gsimp_impute(bnds, imp_model = "nope"), "imp_model")
  expect_error(gsimp_impute(bnds, initial = "nope"), "initial")
})

test_that("all imp_model choices run, including a custom function", {
  set.seed(12)
  bnds <- make_bounds()
  for (m in c("ridge", "lm")) {
    expect_false(anyNA(gsimp_impute(bnds, iters_all = 3, imp_model = m)))
  }
  # custom model: predict column mean, fixed sd
  my_model <- function(yo, Xo, Xm) list(mean = rep(mean(yo), nrow(Xm)),
                                        sd = stats::sd(yo))
  expect_false(anyNA(gsimp_impute(bnds, iters_all = 3, imp_model = my_model)))
})

test_that("glmnet backend runs when available", {
  skip_if_not_installed("glmnet")
  skip_on_cran()
  set.seed(13)
  # needs >= 2 predictor columns
  cds <- lapply(1:3, function(j) {
    v <- exp(rnorm(30))
    lim <- stats::quantile(v, 0.3)
    left <- ifelse(v < lim, 0, v)
    right <- ifelse(v < lim, lim, v)
    as_interval_data(left, right, log_transform = TRUE)
  })
  bnds <- build_bounds(cds)
  # cv.glmnet emits expected small-fold notices on this toy size; not our concern
  filled <- suppressWarnings(
    gsimp_impute(bnds, iters_all = 3, imp_model = "glmnet")
  )
  expect_false(anyNA(filled))
})

# ---- recovery test (the section-6 acceptance criterion) -------------------

test_that("gsimp_impute recovers ground truth from three-tier censoring", {
  set.seed(123)
  n <- 200L
  p <- 4L
  mu <- c(0.5, 1.0, 1.5, 0.8)          # log-scale analyte means
  R <- matrix(0.6, p, p); diag(R) <- 1 # correlated analytes
  z_true <- matrix(stats::rnorm(n * p), n, p) %*% chol(R)
  z_true <- sweep(z_true, 2, mu, "+")  # ground-truth log-concentrations
  conc <- exp(z_true)                  # original scale

  # Censor each analyte into the three tiers at its own MDL / LCMRL.
  mdl <- apply(conc, 2, stats::quantile, 0.25)
  lcmrl <- apply(conc, 2, stats::quantile, 0.40)
  cds <- lapply(seq_len(p), function(j) {
    cj <- conc[, j]
    left <- ifelse(cj < mdl[j], 0, ifelse(cj < lcmrl[j], mdl[j], cj))
    right <- ifelse(cj < mdl[j], mdl[j], ifelse(cj < lcmrl[j], lcmrl[j], cj))
    as_interval_data(left, right, log_transform = TRUE)
  })
  bnds <- build_bounds(cds)

  # A good chunk of cells should be censored (else the test is vacuous).
  expect_gt(mean(bnds$to_impute), 0.30)

  filled <- gsimp_impute(bnds, iters_all = 30, imp_model = "ridge")

  # (a) every imputed value lands inside its assigned interval
  m <- bnds$to_impute
  expect_true(all(filled[m] <= bnds$hi_mat[m] + 1e-8))
  lo_ok <- is.infinite(bnds$lo_mat[m]) | filled[m] >= bnds$lo_mat[m] - 1e-8
  expect_true(all(lo_ok))

  # (b) column means recover the (log-scale) ground truth within tolerance,
  #     and beat naive substitution at the detection limit (upward-biased).
  true_means <- colMeans(z_true)
  gs_means <- colMeans(filled)
  naive <- bnds$data_wide
  naive[m] <- bnds$hi_mat[m]            # substitute the upper bound
  naive_means <- colMeans(naive)
  expect_lt(max(abs(gs_means - true_means)), 0.30)
  expect_lt(sum(abs(gs_means - true_means)),
            sum(abs(naive_means - true_means)))

  # (c) correlation structure is recovered better than naive substitution
  off <- upper.tri(diag(p))
  err_gs <- mean(abs((cor(filled) - cor(z_true))[off]))
  err_naive <- mean(abs((cor(naive) - cor(z_true))[off]))
  expect_lt(err_gs, 0.20)
  expect_lt(err_gs, err_naive)
})

test_that("single-limit (two-tier) left-censoring is handled as a special case", {
  set.seed(321)
  v <- exp(stats::rnorm(120, mean = 1, sd = 0.7))
  dl <- stats::quantile(v, 0.35)            # one detection limit, no DNQ band
  cd <- as_interval_data_from_single(
    value = ifelse(v < dl, dl, v),
    censored = v < dl,
    log_transform = TRUE
  )
  bnds <- build_bounds(cd)
  filled <- gsimp_impute(bnds, iters_all = 15)
  m <- bnds$to_impute
  # every imputed non-detect sits below the detection limit
  expect_true(all(filled[m] <= log(dl) + 1e-8))
  expect_false(anyNA(filled))
})
