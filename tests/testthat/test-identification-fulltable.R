# Tests for the mediator-identification simulation study and the per-method
# (Bonferroni / BH / BY) global-statistic reporting added for the article's
# FWER/FDR and "full-table" analyses.

test_that("pe_identification_study returns FWER/FDR/precision/recall by method", {
  set.seed(113)
  id <- pe_identification_study(n = 120, p = 40, outcome = "continuous",
                                pattern = "contrasting", c1_grid = c(0, 1),
                                n_rep = 8, methods = c("Bonferroni", "BH"))
  expect_s3_class(id, "poem_tbl")
  expect_setequal(names(id),
    c("c1", "method", "fwer", "fdr", "precision", "recall", "n_valid",
      "n_empty"))
  expect_equal(nrow(id), 4L)                       # 2 c1 values x 2 methods
  expect_setequal(id$method, c("Bonferroni", "BH"))
  expect_true(all(id$n_valid <= 8L))
  # Under the article's convention (truth = "outcome_effect", the default)
  # the five mediators with a nonzero outcome coefficient count as active even
  # at c1 = 0, so recall and precision are defined proportions everywhere.
  expect_true(all(id$recall >= 0 & id$recall <= 1))
  expect_true(all(id$precision >= 0 & id$precision <= 1))
  expect_true(all(id$fwer >= 0 & id$fwer <= 1))
  expect_true(all(id$fdr >= 0 & id$fdr <= 1))
  expect_identical(attr(id, "truth"), "outcome_effect")
  # Under truth = "mediation" no mediator is active at c1 = 0 (the exposure
  # paths are zero), so recall is undefined there and every selection is a
  # false positive.
  set.seed(113)
  idm <- pe_identification_study(n = 120, p = 40, outcome = "continuous",
                                 pattern = "contrasting", c1_grid = c(0, 1),
                                 n_rep = 8, methods = "Bonferroni",
                                 truth = "mediation")
  expect_true(all(is.na(idm$recall[idm$c1 == 0])))
  expect_true(all(idm$recall[idm$c1 == 1] >= 0 & idm$recall[idm$c1 == 1] <= 1))
  expect_identical(attr(idm, "truth"), "mediation")
})

test_that("pe_selection reports the PE statistic and p-value for each method", {
  d <- simulate_mediation_data(n = 200, p = 60, outcome = "continuous",
                               pattern = "contrasting", c1 = 1, seed = 113)
  fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous",
                      report_all_methods = TRUE)
  sel <- pe_selection(fit)
  expect_setequal(names(sel),
    c("method", "n_active", "j_pe", "stat_pe", "pval_pe", "active_mediators"))
  expect_setequal(sel$method, c("Bonferroni", "BH", "BY"))
  # The primary method (Bonferroni) p-value matches the fit's own pval_pe.
  prim <- fit$value[fit$term == "pval_pe"]
  expect_equal(sel$pval_pe[sel$method == "Bonferroni"], prim, tolerance = 1e-10)
  expect_true(all(sel$pval_pe >= 0 & sel$pval_pe <= 1))
  # stat_pe = benchmark S_n + J_m, so it is at least the benchmark statistic.
  expect_true(all(sel$stat_pe >= sel$j_pe))
})

test_that("pe_selection still works for a single-method fit", {
  d <- simulate_mediation_data(n = 200, p = 60, outcome = "continuous",
                               pattern = "contrasting", c1 = 1, seed = 113)
  fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous", method = "BH")
  sel <- pe_selection(fit)
  expect_equal(nrow(sel), 1L)
  expect_equal(sel$method, "BH")
  expect_equal(sel$pval_pe, fit$value[fit$term == "pval_pe"], tolerance = 1e-10)
})

test_that("WHO_mediation_analysis(full_table = TRUE) emits per-method columns", {
  res <- WHO_mediation_analysis("imr", groupings = "global", full_table = TRUE,
                                lambda_grid = seq(0.1, 10, length.out = 30))
  expect_setequal(names(res),
    c("grouping", "group", "n_countries", "pval_hdmm",
      "pval_pe_bonferroni", "active_bonferroni",
      "pval_pe_bh", "active_bh", "pval_pe_by", "active_by"))
  expect_equal(res$group, "ALL")
  # The global infant-mortality model recovers gge_gdp under Bonferroni.
  expect_true(grepl("gge_gdp", res$active_bonferroni))
  expect_lt(res$pval_pe_bonferroni, 1e-3)
  expect_gt(res$pval_hdmm, res$pval_pe_bonferroni)   # PE is far smaller
})
