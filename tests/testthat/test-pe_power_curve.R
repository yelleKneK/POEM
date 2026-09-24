test_that("the power curve returns the documented structure", {
  set.seed(113)
  pc <- pe_power_curve(n = 80, p = 30, outcome = "continuous",
                       pattern = "contrasting", c1_grid = c(0, 0.6),
                       n_rep = 8)
  expect_s3_class(pc, "poem_tbl")
  expect_equal(nrow(pc), 2L)
  expect_true(all(c("c1", "rejection_hdmm", "rejection_pe", "n_valid") %in%
                    names(pc)))
  # Rejection rates are proportions in [0, 1].
  expect_true(all(pc$rejection_hdmm >= 0 & pc$rejection_hdmm <= 1))
  expect_true(all(pc$rejection_pe >= 0 & pc$rejection_pe <= 1))
})

test_that("under contrasting mediation the PE test out-rejects the Wald test as c1 grows", {
  skip_on_cran()
  set.seed(113)
  # The article's design (n = 300, p = 500) at c1 = 1: its Figure 1(b) shows
  # the PE test rejecting in about 76 percent of replications and the
  # benchmark in about 16 (the benchmark rises above its size there because
  # the penalized fit sometimes keeps only one member of a cancelling pair);
  # the default grid reproduces both (0.80 and 0.16 to 0.22 in 60- and
  # 200-replication runs, 2026-09-25). At 40 replications the margin of 0.3
  # sits about two Monte Carlo standard errors below that gap.
  pc <- pe_power_curve(n = 300, p = 500, outcome = "continuous",
                       pattern = "contrasting", c1_grid = c(0, 1),
                       n_rep = 40)
  null_row <- pc[pc$c1 == 0, ]
  alt_row  <- pc[pc$c1 == 1, ]
  # Size: both tests are near nominal under the global null (generous band).
  expect_lt(null_row$rejection_pe, 0.2)
  # Power: the PE test rejects far more often than the Wald test under the
  # contrasting alternative, which is the article's headline result.
  expect_gt(alt_row$rejection_pe, alt_row$rejection_hdmm + 0.3)
})

test_that("the power curve validates alpha_level", {
  expect_error(pe_power_curve(n = 50, p = 20, c1_grid = 0, n_rep = 2,
                              alpha_level = 0), "in \\(0, 1\\)")
})

test_that("the study functions forward the tuning grid and count empty fits", {
  set.seed(113)
  big <- c(5, 10)                                  # selects nothing at any n
  pc <- pe_power_curve(n = 80, p = 30, outcome = "continuous",
                       pattern = "contrasting", c1_grid = c(0, 1), n_rep = 4,
                       lambda_grid = big)
  expect_true("n_empty" %in% names(pc))
  expect_equal(pc$n_empty, c(4L, 4L))
  expect_equal(pc$n_valid, c(4L, 4L))
  expect_equal(pc$rejection_pe, c(0, 0))          # an empty fit never rejects
  expect_equal(attr(pc, "lambda_grid"), big)
  # The default records the rescaled article grid.
  pc2 <- pe_power_curve(n = 80, p = 30, outcome = "continuous",
                        pattern = "contrasting", c1_grid = 1, n_rep = 2)
  expect_equal(attr(pc2, "lambda_grid"), pe_lambda_grid(80, 30, "continuous"))
  ss <- ss_power_pe_mediation(outcome = "continuous", pattern = "contrasting",
                              p = 30, n_grid = c(60, 80), n_rep = 2,
                              lambda_grid = big, seed = 113)
  expect_equal(ss$n_empty, c(2L, 2L))
  st <- pe_simulation_study(n = 80, p = 30, outcome = "continuous",
                            c1_grid = 1, n_rep = 2, lambda_grid = big,
                            seed = 113)
  expect_true(all(st$n_empty == 2L))
  expect_equal(attr(st, "lambda_grid"), big)
})

test_that("the study functions reject invalid designs instead of returning NaN", {
  expect_error(pe_power_curve(n = 60, p = 20, n_rep = 3, error_level = 2),
               "`error_level` must be")
  expect_error(pe_power_curve(n = 60, p = 20, n_rep = 0), "`n_rep` must be")
  expect_error(pe_power_curve(n = 0, p = 20, n_rep = 3), "`n` must be")
  expect_error(pe_power_curve(n = 60, p = 20, n_rep = 3, c1_grid = c(NA, 1)),
               "`c1_grid` must be")
  expect_error(pe_power_curve(n = 60, p = 20, n_rep = 3,
                              outcome_args = list(seed = 1)),
               "may not contain `seed`")
  expect_error(pe_power_curve(n = 60, p = 20, n_rep = 3,
                              outcome_args = list(n = 40)),
               "may not contain `n`")
  expect_error(pe_power_curve(n = 60, p = 20, n_rep = 3, lambda_grid = c(0, 1)),
               "`lambda_grid` must be")
  expect_error(pe_simulation_study(n = 60, p = 20, n_rep = 3, error_level = 0),
               "`error_level` must be")
  expect_error(ss_power_pe_mediation(p = 20, n_grid = c(50, 60), n_rep = 3,
                                     c1 = c(0.5, 1)),
               "`c1` must be a single")
  expect_error(ss_power_pe_mediation(p = 20, n_grid = c(50, 60), n_rep = 3,
                                     error_level = 2),
               "`error_level` must be")
  expect_error(pe_identification_study(n = 60, p = 20, n_rep = 3,
                                       error_level = 2),
               "`error_level` must be")
})
