# The seed contract of every exported function with a `seed` argument: a
# supplied seed reproduces the run bit for bit, the caller's `.Random.seed` is
# byte-identical afterwards whether or not one existed before the call, the
# generator kind is restored, and `seed = NULL` leaves the stream advancing.
# Verified at runtime, never by reading the source (QC master, Randomness).

seeded_calls <- list(
  simulate_mediation_data = function(seed)
    simulate_mediation_data(n = 40, p = 12, outcome = "continuous",
                            pattern = "contrasting", seed = seed),
  simulate_guo_mediation = function(seed)
    simulate_guo_mediation(c1 = 0.5, n = 30, seed = seed),
  pe_power_curve = function(seed)
    pe_power_curve(n = 60, p = 15, outcome = "continuous",
                   pattern = "contrasting", c1_grid = 1, n_rep = 2, seed = seed),
  pe_simulation_study = function(seed)
    pe_simulation_study(n = 60, p = 15, outcome = "continuous", c1_grid = 1,
                        n_rep = 2, seed = seed),
  pe_identification_study = function(seed)
    pe_identification_study(n = 60, p = 15, outcome = "continuous",
                            pattern = "contrasting", c1_grid = 1, n_rep = 2,
                            methods = "Bonferroni", seed = seed),
  ss_power_pe_mediation = function(seed)
    ss_power_pe_mediation(outcome = "continuous", pattern = "contrasting",
                          p = 15, n_grid = c(40, 60), n_rep = 2, seed = seed))

strip <- function(x) { attributes(x) <- NULL; x }

for (nm in names(seeded_calls)) {
  f <- seeded_calls[[nm]]
  test_that(paste0(nm, "() reproduces a seeded run and restores the RNG state"), {
    # A seed reproduces the run bit for bit.
    a <- f(113); b <- f(113)
    expect_identical(a, b)
    # A different seed gives different data (the seed is not ignored); the
    # Monte Carlo tables are too coarse at two replications to assert this.
    if (nm %in% c("simulate_mediation_data", "simulate_guo_mediation"))
      expect_false(identical(strip(unclass(a)), strip(unclass(f(114)))))
    # Existing state: byte-identical afterwards, kind untouched.
    set.seed(999)
    before <- get(".Random.seed", envir = globalenv())
    kind_before <- RNGkind()
    invisible(f(113))
    expect_identical(get(".Random.seed", envir = globalenv()), before)
    expect_identical(RNGkind(), kind_before)
    # Absent state: still absent afterwards.
    if (exists(".Random.seed", envir = globalenv()))
      rm(".Random.seed", envir = globalenv())
    invisible(f(113))
    expect_false(exists(".Random.seed", envir = globalenv()))
    # seed = NULL draws from the user's stream and advances it.
    set.seed(999); invisible(f(NULL))
    after <- get(".Random.seed", envir = globalenv())
    set.seed(999)
    expect_false(identical(after, get(".Random.seed", envir = globalenv())))
  })
}

test_that("a seeded call restores a non-default generator kind", {
  old <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old)), add = TRUE)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(5)
  before <- get(".Random.seed", envir = globalenv())
  invisible(seeded_calls$pe_power_curve(113))
  expect_identical(RNGkind()[1L], "L'Ecuyer-CMRG")
  expect_identical(get(".Random.seed", envir = globalenv()), before)
})

test_that("`seed` must be a single number or NULL", {
  expect_error(simulate_mediation_data(n = 40, p = 12, seed = "a"), "`seed` must be")
  expect_error(simulate_mediation_data(n = 40, p = 12, seed = c(1, 2)), "`seed` must be")
  expect_error(pe_power_curve(n = 40, p = 12, n_rep = 2, seed = NA), "`seed` must be")
})

test_that("a seeded parallel run reproduces itself and restores the RNG state", {
  skip_on_cran()
  skip_on_os("windows")
  par_call <- function(seed, cores)
    pe_power_curve(n = 60, p = 15, outcome = "continuous", pattern = "contrasting",
                   c1_grid = 1, n_rep = 4, cores = cores, seed = seed)
  a <- par_call(113, 2L); b <- par_call(113, 2L)
  expect_identical(a, b)
  set.seed(999)
  before <- get(".Random.seed", envir = globalenv()); kind_before <- RNGkind()
  invisible(par_call(113, 2L))
  expect_identical(get(".Random.seed", envir = globalenv()), before)
  expect_identical(RNGkind(), kind_before)
  if (exists(".Random.seed", envir = globalenv()))
    rm(".Random.seed", envir = globalenv())
  invisible(par_call(113, 2L))
  expect_false(exists(".Random.seed", envir = globalenv()))
  # pe_identification_study takes the same path
  id1 <- pe_identification_study(n = 60, p = 15, outcome = "continuous",
                                 pattern = "contrasting", c1_grid = 1, n_rep = 4,
                                 methods = "Bonferroni", cores = 2L, seed = 113)
  id2 <- pe_identification_study(n = 60, p = 15, outcome = "continuous",
                                 pattern = "contrasting", c1_grid = 1, n_rep = 4,
                                 methods = "Bonferroni", cores = 2L, seed = 113)
  expect_identical(id1, id2)
  # an unseeded parallel call switches the kind for the call only
  RNGkind("Mersenne-Twister"); kind_before <- RNGkind()
  invisible(par_call(NULL, 2L))
  expect_identical(RNGkind(), kind_before)
})
