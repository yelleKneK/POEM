# Benchmark numerical-correctness tests. WHO_health_mediation is shipped as a
# benchmark data set and WHO_mediation_analysis() as the benchmark fit; these
# tests pin the fit to the active mediators reported by the article
# (in press), so a future change to the estimation internals that
# silently moved the published findings would fail here. The grid below
# (seq(0.1, 10, length.out = 40)) is the one used to certify the benchmark.

bench_lambda <- seq(0.1, 10, length.out = 40)

# Convenience: the active-mediator string for one outcome/region/income cell.
bench_active <- function(outcome, grouping, group, ...) {
  res <- WHO_mediation_analysis(outcome, groupings = grouping,
                                lambda_grid = bench_lambda, ...)
  res$active_mediators[res$group == group]
}

test_that("benchmark fit reproduces the infant-mortality findings of the article", {
  skip_on_cran()
  res <- WHO_mediation_analysis("imr", groupings = c("global", "region"),
                                lambda_grid = bench_lambda)
  active <- function(g) strsplit(res$active_mediators[res$group == g], ", ")[[1]]

  # The article's Table 1 gives the complete active set per cell, so each is
  # pinned exactly (a spurious extra mediator would fail here).
  # Global model: general government expenditure as a percent of GDP.
  expect_setequal(active("ALL"), "gge_gdp")
  # Western Pacific: the benchmark Wald test sees nothing, the PE test detects
  # compulsory health insurance and private expenditure.
  expect_gt(res$pval_hdmm[res$group == "WPR"], 0.5)
  expect_setequal(active("WPR"), c("chi_che", "pvtd_gdp"))
  # South-East Asia: external expenditure in constant 2021 US dollars per capita.
  expect_setequal(active("SEAR"), "ext_usd2021_pc")
  # Africa and the Americas: no active mediator detected.
  expect_identical(active("AFR"), "none")
  expect_identical(active("AMR"), "none")
})

test_that("benchmark fit reproduces the income-group infant-mortality findings", {
  skip_on_cran()
  res <- WHO_mediation_analysis("imr", groupings = "income",
                                lambda_grid = bench_lambda)
  active <- function(g) strsplit(res$active_mediators[res$group == g], ", ")[[1]]
  expect_setequal(active("Low"), "pvtd_usd2021")        # private expenditure
  expect_setequal(active("High"), c("oops_che", "shi_che"))  # out-of-pocket, social insurance
  expect_identical(active("Lower-middle"), "none")
  expect_identical(active("Upper-middle"), "none")
})

test_that("benchmark fit reproduces the global undernourishment finding", {
  skip_on_cran()
  pou <- bench_active("pou", "global", "ALL")
  # Social health insurance as a percent of current health expenditure: the
  # benchmark Wald test misses it, the PE test detects it (Yu and Kelley, in
  # press, Table 2).
  expect_identical(pou, "shi_che")
})
