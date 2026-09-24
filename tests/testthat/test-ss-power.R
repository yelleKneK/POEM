# Regression tests for ss_power_pe_mediation(): the default (unspecified)
# pattern call must work, every supported pattern must be accepted, and the
# RNG state must be restored on exit in both the existing- and absent-seed
# cases (the @param seed contract).

test_that("ss_power_pe_mediation runs with the pattern left at its default", {
  plan <- ss_power_pe_mediation(p = 10, n_grid = c(40, 60), n_rep = 2,
                                seed = 113)
  expect_s3_class(plan, "poem_tbl")
  expect_equal(nrow(plan), 2L)
})

test_that("ss_power_pe_mediation accepts each supported pattern", {
  for (pat in c("homogeneous", "contrasting", "heterogeneous"))
    expect_s3_class(
      ss_power_pe_mediation(pattern = pat, p = 10, n_grid = 40, n_rep = 2,
                            seed = 1),
      "poem_tbl")
})

test_that("ss_power_pe_mediation restores the RNG state on exit", {
  # Absent-seed case: .Random.seed must remain absent afterwards.
  if (exists(".Random.seed", envir = globalenv()))
    rm(".Random.seed", envir = globalenv())
  invisible(ss_power_pe_mediation(p = 10, n_grid = 40, n_rep = 2, seed = 113))
  expect_false(exists(".Random.seed", envir = globalenv()))

  # Existing-seed case: the caller's stream is left untouched.
  set.seed(999)
  before <- get(".Random.seed", envir = globalenv())
  invisible(ss_power_pe_mediation(p = 10, n_grid = 40, n_rep = 2, seed = 113))
  expect_identical(get(".Random.seed", envir = globalenv()), before)
})
