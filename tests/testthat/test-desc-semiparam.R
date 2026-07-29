# Tests for the semi-parametric descriptive module (desc-semiparam.R).

# A dataset with a clear group effect, so the covariate coefficient is
# well-identified and its sign is predictable.
make_grouped <- function(log_transform = FALSE, seed = 42) {
  set.seed(seed)
  n <- 80
  grp <- factor(rep(c("low", "high"), each = n / 2))
  truth <- rlnorm(n, meanlog = ifelse(grp == "high", 1.4, 0.4), sdlog = 0.5)
  lo <- ifelse(truth < 1, 0, floor(truth))
  hi <- ifelse(truth < 1, 1, ceiling(truth))
  list(
    x = as_interval_data(lo, hi, log_transform = log_transform),
    grp = grp
  )
}

test_that("desc_sp returns a cens_sp_fit with expected metadata", {
  d <- make_grouped()
  fit <- desc_sp(d$x, covariates = d$grp, model = "ph")
  expect_s3_class(fit, "cens_sp_fit")
  expect_equal(fit$model, "ph")
  expect_equal(fit$covariates, "covariate")
  expect_equal(fit$n, 80)
})

test_that("desc_sp requires a cens_data object", {
  expect_error(desc_sp(1:5), "cens_data")
})

test_that("intercept-only model fits and has no coefficients", {
  d <- make_grouped()
  fit <- desc_sp(d$x, covariates = NULL)
  expect_length(coef(fit), 0)
  expect_output(print(fit), "Intercept-only")
})

test_that("covariate effect has the expected sign (higher group -> larger values)", {
  d <- make_grouped()
  fit <- desc_sp(d$x, covariates = d$grp, model = "ph")
  cf <- coef(fit)
  expect_length(cf, 1)
  # 'low' is the reference; 'high' has larger concentrations, hence *lower*
  # hazard of the (small) event time -> negative PH coefficient.
  expect_true(cf[["covariatelow"]] > 0 || cf[["covariatehigh"]] < 0)
})

test_that("both models run", {
  d <- make_grouped()
  for (m in c("ph", "po")) {
    fit <- desc_sp(d$x, covariates = d$grp, model = m)
    expect_s3_class(fit, "cens_sp_fit")
    expect_length(coef(fit), 1)
  }
})

test_that("a data frame of covariates is accepted with its own names", {
  d <- make_grouped()
  cov <- data.frame(region = d$grp)
  fit <- desc_sp(d$x, covariates = cov, model = "ph")
  expect_equal(fit$covariates, "region")
  expect_equal(names(coef(fit)), "regionlow")
})

test_that("mismatched covariate length errors", {
  d <- make_grouped()
  expect_error(
    desc_sp(d$x, covariates = factor(c("a", "b"))),
    "one row/element per observation"
  )
})

test_that("covariates are subset in step with dropped bounds", {
  x <- as_interval_data(left = c(0, 1, NA, 3), right = c(1, 3, NA, 3),
                        log_transform = FALSE)
  grp <- factor(c("a", "b", "a", "b"))
  fit <- desc_sp(x, covariates = grp)
  expect_equal(fit$n, 3)
  expect_equal(fit$dropped, 1)
})

test_that("bootstrap standard errors are produced when requested", {
  skip_on_cran()
  d <- make_grouped()
  fit <- desc_sp(d$x, covariates = d$grp, model = "ph", bs_samples = 20)
  expect_gt(length(fit$fit$bsMat), 0)
})

test_that("plot works for baseline and per-group curves", {
  d <- make_grouped()
  fit <- desc_sp(d$x, covariates = d$grp, model = "ph")
  tmp <- tempfile(fileext = ".pdf")
  pdf(tmp)
  on.exit({ dev.off(); unlink(tmp) }, add = TRUE)
  expect_invisible(plot(fit))  # baseline
  nd <- data.frame(covariate = factor(c("low", "high"),
                                      levels = c("high", "low")))
  expect_invisible(plot(fit, newdata = nd))  # per-group
})
