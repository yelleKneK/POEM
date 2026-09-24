test_that("WHO_mediation_analysis reproduces the global infant-mortality finding", {
  skip_on_cran()                       # fits a penalized model on 2002 x 57
  res <- WHO_mediation_analysis("imr", groupings = "global",
                                lambda_grid = seq(0.1, 10, length.out = 30))
  expect_s3_class(res, "data.frame")
  expect_true(all(c("grouping", "group", "n_countries", "pval_hdmm", "pval_pe",
                    "n_active", "active_mediators") %in% names(res)))
  expect_equal(res$group, "ALL")
  expect_equal(res$n_countries, 91L)
  # The power-enhanced test rejects and identifies general government
  # expenditure as a percent of GDP, matching the article.
  expect_lt(res$pval_pe, 0.001)
  expect_match(res$active_mediators, "gge_gdp")
})

test_that("WHO_mediation_analysis runs the region grouping and tolerates empty cells", {
  skip_on_cran()
  res <- WHO_mediation_analysis("imr", groupings = c("global", "region"),
                                lambda_grid = seq(0.1, 10, length.out = 30))
  expect_setequal(res$group, c("ALL", "AFR", "AMR", "EMR", "EUR", "SEAR", "WPR"))
  # Western Pacific is a case where the benchmark test sees nothing but the PE
  # test detects active mediators (the article's headline pattern).
  wpr <- res[res$group == "WPR", ]
  expect_gt(wpr$pval_hdmm, 0.5)
  expect_lt(wpr$pval_pe, 0.05)
})

test_that("WHO_mediation_analysis validates its arguments", {
  expect_error(WHO_mediation_analysis("nope"), "should be one of")
  expect_error(WHO_mediation_analysis("imr", groupings = "country"),
               "should be one of")
  expect_error(WHO_mediation_analysis("imr", error_level = 2),
               "`error_level` must be")
  expect_error(WHO_mediation_analysis("imr", lambda_grid = c(0, 1)),
               "`lambda_grid` must be")
  expect_error(WHO_mediation_analysis("imr", lambda_grid_global = numeric(0)),
               "`lambda_grid_global` must be")
  expect_error(WHO_mediation_analysis("imr", cores = 0), "`cores` must be")
  expect_error(WHO_mediation_analysis("imr", verbose = NA), "`verbose` must be")
  expect_error(WHO_mediation_analysis("imr", full_table = "yes"),
               "`full_table` must be")
  # A wrong `data` stops up front instead of every cell coming back NA.
  expect_error(WHO_mediation_analysis("imr", data = iris), "lacks the column")
  expect_error(WHO_mediation_analysis("imr", data = "x"), "must be a data.frame")
})

test_that("WHO_mediation_analysis names the groups it could not fit", {
  skip_on_cran()
  # Remove one region's rows: its cell has no observations, is skipped,
  # reported as NA, named in a message, and recorded in the notes.
  sub <- WHO_health_mediation[WHO_health_mediation$region != "SEAR", ]
  # (The coarse five-point grid can also make a small region's fit fail; that
  # is reported as a warning and is not what this test is about.)
  expect_message(
    res <- suppressWarnings(
      WHO_mediation_analysis("imr", groupings = "region", data = sub,
                             lambda_grid = seq(0.5, 10, length.out = 5))),
    "SEAR")
  expect_true(is.na(res$pval_pe[res$group == "SEAR"]))
  notes <- attr(res, "notes")
  expect_s3_class(notes, "data.frame")
  expect_true("SEAR" %in% notes$group[notes$stage == "design"])
  # A run where every group fits carries no notes.
  full <- WHO_mediation_analysis("imr", groupings = "global",
                                 lambda_grid = seq(0.5, 10, length.out = 5))
  expect_null(attr(full, "notes"))
})
