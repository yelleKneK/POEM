# Tests for summary.poem_tbl, the cross-group digest of any grouped POEM comparison table.

test_that("summary.poem_tbl surfaces where PE beats the benchmark (default table)", {
  res <- WHO_mediation_analysis("imr", groupings = c("global", "region"),
                                lambda_grid = seq(0.1, 10, length.out = 30))
  expect_s3_class(res, "poem_tbl")
  s <- summary(res)
  expect_s3_class(s, "summary.poem_tbl")
  expect_identical(s$outcome, "imr")
  expect_gte(s$n_pe, s$n_hdmm)                       # PE detects at least as many
  # SEAR and WPR are the PE-only rescues in the article's IMR table.
  expect_true(all(c("SEAR", "WPR") %in% s$pe_only$group))
  expect_true(all(s$pe_only$pval_hdmm > 0.05))       # benchmark missed them
  expect_true(all(s$pe_only$pval_pe <= 0.05))        # PE caught them
  expect_true("gge_gdp" %in% names(s$top_mediators))
  expect_error(summary(res, alpha_level = 0), "`alpha_level` must be")
  expect_output(print(s), "alpha_level = 0.05")
})

test_that("summary.poem_tbl also works on the extended (full_table) layout", {
  res <- WHO_mediation_analysis("imr", groupings = c("global", "region"),
                                full_table = TRUE,
                                lambda_grid = seq(0.1, 10, length.out = 30))
  s <- summary(res)
  expect_true(s$full_table)
  expect_true(all(c("SEAR", "WPR") %in% s$pe_only$group))
  # print method runs without error
  expect_output(print(s), "POEM comparison")
})
