test_that("the shipped WHO data has the documented structure", {
  data(WHO_health_mediation, package = "POEM")
  data(WHO_indicator_codebook, package = "POEM")
  expect_equal(nrow(WHO_health_mediation), 2002L)
  expect_equal(length(unique(WHO_health_mediation$code)), 91L)
  expect_equal(range(WHO_health_mediation$year), c(2000L, 2021L))
  # The 57 indicators are present, complete, and match the codebook.
  ind <- setdiff(names(WHO_health_mediation),
                 c("code", "country", "region", "income", "year",
                   "gdp_growth", "imr", "u5mr", "leb", "lbw", "pou"))
  expect_equal(length(ind), 57L)
  expect_false(anyNA(WHO_health_mediation[, ind]))
  expect_setequal(ind, WHO_indicator_codebook$indicator)
  # Region tallies match the article.
  reg <- tapply(WHO_health_mediation$code, WHO_health_mediation$region,
                function(x) length(unique(x)))
  expect_equal(reg[["AFR"]], 32L); expect_equal(reg[["EUR"]], 9L)
  expect_equal(reg[["SEAR"]], 5L)
})

test_that("WHO_mediation_design assembles a usable design", {
  des <- WHO_mediation_design("imr")
  expect_equal(des$n, 2002L)
  expect_equal(ncol(des$M), 57L)
  expect_equal(ncol(des$X), 1L)
  expect_equal(length(des$Y), des$n)
  # Confounders: 5 region + 3 income indicators + year = 9 columns.
  expect_equal(ncol(des$Z), 9L)
  expect_false(anyNA(des$M))
})

test_that("subsetting drops constant covariates and unused factor levels", {
  eur <- WHO_mediation_design("imr", region = "EUR")
  expect_equal(eur$n, 9L * 22L)
  # Region is constant within EUR, so Z has no region columns; the design is
  # full rank (no all-zero columns from unused factor levels).
  expect_false(any(grepl("region", colnames(eur$Z))))
  expect_equal(qr(cbind(eur$X, eur$Z))$rank, ncol(eur$Z) + 1L)
})

test_that("outcomes with narrower coverage yield fewer observations", {
  expect_equal(WHO_mediation_design("pou")$n, 79L * 21L)   # 2001-2021
  expect_equal(WHO_mediation_design("lbw")$n, 73L * 21L)   # 2000-2020
})

test_that("WHO_mediation_design validates its arguments", {
  expect_error(WHO_mediation_design("imr", region = "XYZ"), "Unknown region")
  expect_error(WHO_mediation_design("imr", income = "rich"), "Unknown income")
  expect_error(WHO_mediation_design("nope"), "should be one of")
})
