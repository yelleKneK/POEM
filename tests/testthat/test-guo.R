# Tests for the real-data-motivated heterogeneous setting (simulate_guo_mediation
# and the guo_calibration constants).

test_that("guo_calibration coefficients give the documented total indirect effect", {
  expect_equal(length(guo_calibration$alpha_m), 11L)
  expect_equal(length(guo_calibration$Gamma_x), 11L)
  expect_equal(dim(guo_calibration$Gamma_z), c(11L, 9L))
  # beta_per_c1 is EXACTLY sum(Gamma_x * alpha_m) from the shipped (rounded)
  # coefficients -- this is the value the package actually produces.
  expect_equal(guo_calibration$beta_per_c1,
               sum(guo_calibration$Gamma_x * guo_calibration$alpha_m))
  # It prints as -1.597 to four decimals (the headline number a user sees);
  # the article's -1.5977 comes from the unrounded coefficients.
  expect_equal(round(guo_calibration$beta_per_c1, 4), -1.597)
  expect_equal(guo_calibration$beta_per_c1, -1.5977, tolerance = 1e-3)
})

test_that("simulate_guo_mediation beta equals c1 * beta_per_c1 exactly", {
  for (k in c(0.3, 0.5, 1)) {
    d <- simulate_guo_mediation(c1 = k, seed = 1)
    expect_equal(d$beta, k * guo_calibration$beta_per_c1)
  }
})

test_that("simulate_guo_mediation matches the calibrated dimensions and beta", {
  d <- simulate_guo_mediation(c1 = 0.5, seed = 113)
  expect_equal(dim(d$M), c(85L, 1008L))
  expect_equal(dim(d$X), c(85L, 1L))
  expect_null(d$Z)
  expect_equal(d$active_mediators, 1:11)
  expect_equal(d$beta, guo_calibration$beta_per_c1 * 0.5, tolerance = 1e-8)
})

test_that("there is no mediation at the null, and confounders are supported", {
  expect_length(simulate_guo_mediation(c1 = 0, seed = 1)$active_mediators, 0L)
  dc <- simulate_guo_mediation(c1 = 0.5, confounders = TRUE, seed = 7)
  expect_equal(dim(dc$Z), c(85L, 8L))   # 8 covariates returned (intercept absorbed)
  expect_equal(dc$beta, guo_calibration$beta_per_c1 * 0.5, tolerance = 1e-8)
})

test_that("alpha_m variants and a custom vector are accepted", {
  d1 <- simulate_guo_mediation(c1 = 0.5, alpha_m = "contrasting_like", seed = 1)
  expect_equal(d1$alpha_m[1:4], c(1, -0.5, 0.4, -0.3))
  d2 <- simulate_guo_mediation(c1 = 0.5, alpha_m = rep(0.5, 11), seed = 1)
  expect_equal(unique(d2$alpha_m[1:11]), 0.5)
  expect_error(simulate_guo_mediation(c1 = 0.5, alpha_m = 1:3),
               "must have length 11")
})

test_that("the power-enhanced screen fires on the heterogeneous Guo signal", {
  set.seed(113)
  d <- simulate_guo_mediation(c1 = 1, n = 150)
  f <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous",
                    lambda_grid = seq(0.1, 10, length.out = 50))
  # The screen identifies active mediators (J_m > 0), so PE strictly improves
  # on the benchmark Wald test under this mixed-sign heterogeneous signal.
  expect_gte(length(attr(f, "active_mediators")), 1L)
  expect_lt(f$value[f$term == "pval_pe"], f$value[f$term == "pval_hdmm"])
})
