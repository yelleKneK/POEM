# pe_lambda_grid(): the article's grids at the article's design, rescaled by
# sqrt(log p / n) elsewhere, and the validation of its arguments.

test_that("at the article's design the grids are the article's scripts' grids", {
  # code-main-HDMM-Linear.R: seq(0.2, 0.39, 20) full, seq(0.27, 0.46, 20) reduced;
  # code-main-HDGMM-Logistic.R: seq(0.04, 0.2, 15); code-main-HDGMM-Poisson.R:
  # seq(0.7, 5, 15). Exact arithmetic up to floating point.
  expect_equal(pe_lambda_grid(300, 500, "continuous"),
               seq(0.2, 0.39, length.out = 20), tolerance = 1e-12)
  expect_equal(pe_lambda_grid(300, 500, "continuous", "reduced"),
               seq(0.27, 0.46, length.out = 20), tolerance = 1e-12)
  expect_equal(pe_lambda_grid(300, 500, "binary"),
               seq(0.04, 0.2, length.out = 15), tolerance = 1e-12)
  expect_equal(pe_lambda_grid(300, 500, "count"),
               seq(0.7, 5, length.out = 15), tolerance = 1e-12)
})

test_that("the endpoints scale with sqrt(log p / n)", {
  ratio <- sqrt(log(200) / 100) / sqrt(log(500) / 300)
  expect_equal(range(pe_lambda_grid(100, 200, "continuous")),
               c(0.2, 0.39) * ratio, tolerance = 1e-12)
  # the reduced grid of a binary or count outcome is the full grid
  expect_identical(pe_lambda_grid(100, 200, "binary", "reduced"),
                   pe_lambda_grid(100, 200, "binary", "full"))
  expect_length(pe_lambda_grid(100, 200, "count", length.out = 7), 7L)
  expect_true(all(diff(pe_lambda_grid(50, 20, "continuous")) > 0))
})

test_that("pe_lambda_grid validates its arguments", {
  expect_error(pe_lambda_grid(1, 20), "`n` must be")
  expect_error(pe_lambda_grid(50, NA), "`p` must be")
  expect_error(pe_lambda_grid(50, 20, length.out = 1), "`length.out`")
  expect_error(pe_lambda_grid(50, 20, length.out = 2.5), "`length.out`")
  expect_error(pe_lambda_grid(50, 20, outcome = "gaussian"), "should be one of")
})

test_that("pe_mediation() searches that grid by default and records the choice", {
  d <- simulate_mediation_data(n = 150, p = 40, outcome = "continuous",
                               pattern = "homogeneous", c1 = 1, seed = 113)
  fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
  tun <- attr(fit, "tuning")
  expect_equal(tun$lambda_grid, pe_lambda_grid(150, 40, "continuous"))
  expect_equal(tun$lambda_grid_reduced,
               pe_lambda_grid(150, 40, "continuous", "reduced"))
  expect_true(tun$lambda_selected %in% tun$lambda_grid)
  expect_true(tun$lambda_selected_reduced %in% tun$lambda_grid_reduced)
  expect_identical(tun$at_lower_end, tun$lambda_selected <= min(tun$lambda_grid))
  # a user grid serves the reduced model too unless a reduced grid is given
  fit2 <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous",
                       lambda_grid = c(0.1, 0.2, 0.3))
  expect_equal(attr(fit2, "tuning")$lambda_grid_reduced, c(0.1, 0.2, 0.3))
  fit3 <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous",
                       lambda_grid = c(0.1, 0.2, 0.3),
                       lambda_grid_reduced = c(0.2, 0.4))
  expect_equal(attr(fit3, "tuning")$lambda_grid_reduced, c(0.2, 0.4))
  expect_output(print(fit), "Tuning parameter \\(HBIC\\): lambda = ")
})

test_that("binary and count fits default to the review-era grid and record it", {
  d <- simulate_mediation_data(n = 200, p = 30, outcome = "binary",
                               pattern = "homogeneous", c1 = 1, c2 = 1,
                               seed = 113)
  fit <- pe_mediation(d$X, d$Y, d$M, outcome = "binary")
  tun <- attr(fit, "tuning")
  expect_equal(tun$lambda_grid, seq(0.05, 1, length.out = 20))
  expect_null(tun$lambda_grid_reduced)
  expect_true(tun$lambda_selected %in% tun$lambda_grid)
  # the article's logistic grid is available on request
  fit2 <- pe_mediation(d$X, d$Y, d$M, outcome = "binary",
                       lambda_grid = pe_lambda_grid(200, 30, "binary"))
  expect_equal(attr(fit2, "tuning")$lambda_grid, pe_lambda_grid(200, 30, "binary"))
})
