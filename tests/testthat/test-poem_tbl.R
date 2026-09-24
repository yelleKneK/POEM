test_that("the display rounds but the stored values keep full precision", {
  set.seed(113)
  d <- simulate_mediation_data(n = 100, p = 40, outcome = "continuous",
                               pattern = "homogeneous", c1 = 0.6)
  fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
  # The stored statistic is a full-precision double, not the rounded display.
  s <- fit$value[fit$term == "stat_pe"]
  expect_true(is.numeric(s))
  expect_false(s == signif(s, 1))      # almost surely not a round number
  # The format method returns character columns; the object itself does not.
  disp <- format(fit)
  expect_type(disp$value, "character")
  expect_type(fit$value, "double")
})

test_that("whole-number rows print without a decimal part", {
  set.seed(113)
  d <- simulate_mediation_data(n = 100, p = 40, outcome = "continuous",
                               pattern = "homogeneous", c1 = 0.6)
  fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
  disp <- format(fit)
  expect_identical(disp$value[fit$term == "n_observations"], "100")
  expect_identical(disp$value[fit$term == "df"], "1")
})

test_that("p-values print to fixed decimals with a floor", {
  # A tiny p-value prints as the floor label, not as 0.0000.
  tab <- data.frame(term = c("pval_pe", "stat_pe"), value = c(1e-15, 12.3))
  tab <- POEM:::.as_poem_tbl(tab, p_terms = "pval_pe")
  disp <- format(tab)
  expect_identical(disp$value[1], "< 0.0001")
  # A moderate p-value prints to four decimals.
  tab2 <- data.frame(term = "pval_pe", value = 0.0317)
  tab2 <- POEM:::.as_poem_tbl(tab2, p_terms = "pval_pe")
  expect_identical(format(tab2)$value[1], "0.0317")
})

test_that("printing returns the object invisibly and is idempotent in class", {
  set.seed(113)
  d <- simulate_mediation_data(n = 80, p = 30, outcome = "continuous",
                               pattern = "homogeneous")
  fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
  expect_invisible(print(fit))
  twice <- POEM:::.as_poem_tbl(fit)
  expect_identical(class(twice), class(fit))   # .as_poem_tbl is idempotent
})
