# Tests for drop_constant: zero-variance mediator columns are an error by
# default and dropped (with the active set renumbered to the original columns)
# when drop_constant = TRUE.

make_data <- function() {
  d <- simulate_mediation_data(n = 200, p = 60, outcome = "continuous",
                               pattern = "contrasting", c1 = 1, seed = 113)
  d$M[, 10] <- 5      # force two constant (zero-variance) mediator columns
  d$M[, 40] <- -2
  d
}

test_that("constant mediator columns error by default and name the columns", {
  d <- make_data()
  expect_error(
    pe_mediation(d$X, d$Y, d$M, outcome = "continuous"),
    "constant \\(zero-variance\\) columns: 10, 40")
})

test_that("drop_constant = TRUE drops them, warns, and keeps original indexing", {
  d <- make_data()
  expect_warning(
    fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous",
                        drop_constant = TRUE),
    "2 constant .* dropped: 10, 40")
  # Two of the original 60 columns are gone, so 58 candidates remain.
  expect_equal(fit$value[fit$term == "n_candidate_mediators"], 58)
  # Reported indices live in the ORIGINAL 1..60 space, never touching 10 or 40.
  active <- attr(fit, "active_mediators")
  expect_false(any(c(10, 40) %in% active))
  expect_true(all(active >= 1 & active <= 60))
  expect_true(all(pe_mediators(fit)$mediator <= 60))
  # The true active mediators (1..4) are recovered, not shifted by the drop.
  expect_true(all(active %in% 1:4))
})

test_that("drop_constant errors when every mediator is constant", {
  d <- simulate_mediation_data(n = 120, p = 20, outcome = "continuous",
                               pattern = "contrasting", c1 = 1, seed = 1)
  d$M[] <- 1
  expect_error(
    pe_mediation(d$X, d$Y, d$M, outcome = "continuous", drop_constant = TRUE),
    "Every column .* is constant")
})

test_that("drop_constant flows through the binary worker too", {
  d <- simulate_mediation_data(n = 250, p = 50, outcome = "binary",
                               pattern = "contrasting", c1 = 1, c2 = 1, seed = 7)
  d$M[, 5] <- 0
  expect_warning(
    fit <- pe_mediation(d$X, d$Y, d$M, outcome = "binary",
                        drop_constant = TRUE),
    "dropped: 5")
  expect_false(5 %in% attr(fit, "active_mediators"))
})
