test_that("the simulator returns the documented structure and truth", {
  set.seed(113)
  d <- simulate_mediation_data(n = 100, p = 40, outcome = "continuous",
                               pattern = "homogeneous", d = 2)
  expect_equal(dim(d$X), c(100L, 1L))
  expect_equal(dim(d$M), c(100L, 40L))
  expect_length(d$Y, 100L)
  expect_equal(dim(d$Z), c(100L, 2L))
  expect_length(d$alpha_m, 40L)
  # beta is the total indirect effect Gamma_x %*% alpha_m, and the active set
  # is exactly the mediators nonzero on both paths.
  expect_equal(d$beta, as.numeric(d$Gamma_x %*% d$alpha_m))
  expect_equal(d$active_mediators,
               which(d$alpha_m != 0 & apply(d$Gamma_x != 0, 2, any)))
})

test_that("the contrasting pattern has mixed-sign active mediators", {
  for (oc in c("continuous", "binary", "count")) {
    set.seed(113)
    d <- simulate_mediation_data(n = 80, p = 30, outcome = oc,
                                 pattern = "contrasting")
    expect_gt(length(d$active_mediators), 1L)
    # Individual indirect effects (path product) of the active mediators carry
    # both signs: this is what "contrasting" means. The continuous and binary
    # presets cancel to an exactly zero total with the default tau; the count
    # preset (from the article) is mixed-sign but does not exactly cancel.
    ind <- (d$Gamma_x %*% diag(d$alpha_m))[1, d$active_mediators]
    expect_true(any(ind > 0) && any(ind < 0), label = paste("signs for", oc))
  }
})

test_that("the continuous and binary contrasting presets cancel to zero total", {
  for (oc in c("continuous", "binary")) {
    set.seed(113)
    d <- simulate_mediation_data(n = 80, p = 30, outcome = oc,
                                 pattern = "contrasting")
    expect_lt(abs(d$beta), 1e-8, label = paste("beta for", oc))
  }
})

test_that("outcome types produce correctly typed responses", {
  set.seed(113)
  b <- simulate_mediation_data(n = 80, p = 30, outcome = "binary",
                               pattern = "homogeneous")
  expect_true(all(b$Y %in% c(0, 1)))
  set.seed(113)
  k <- simulate_mediation_data(n = 80, p = 30, outcome = "count",
                               pattern = "homogeneous")
  expect_true(all(k$Y >= 0 & k$Y == round(k$Y)))
})

test_that("c1 = 0 yields the global null (no exposure-mediator signal)", {
  set.seed(113)
  d <- simulate_mediation_data(n = 80, p = 30, outcome = "continuous",
                               pattern = "homogeneous", c1 = 0)
  expect_true(all(d$Gamma_x == 0))
  expect_equal(d$beta, 0)
})

test_that("seeding is reproducible and restores the caller's RNG state", {
  d1 <- simulate_mediation_data(n = 50, p = 20, outcome = "continuous",
                                pattern = "homogeneous", seed = 7)
  d2 <- simulate_mediation_data(n = 50, p = 20, outcome = "continuous",
                                pattern = "homogeneous", seed = 7)
  expect_identical(d1$M, d2$M)
  expect_identical(d1$Y, d2$Y)
  # The caller's stream is untouched by a seeded call.
  set.seed(42); before <- .Random.seed
  invisible(simulate_mediation_data(n = 30, p = 10, outcome = "continuous",
                                    pattern = "homogeneous", seed = 99))
  expect_identical(.Random.seed, before)
})

test_that("a user-supplied alpha_m overrides the preset", {
  set.seed(113)
  am <- c(2, -2, rep(0, 28))
  d <- simulate_mediation_data(n = 60, p = 30, outcome = "continuous",
                               pattern = "homogeneous", alpha_m = am)
  expect_identical(d$alpha_m, am)
  expect_error(simulate_mediation_data(n = 60, p = 30, alpha_m = c(1, 2)),
               "length")
})

test_that("the simulator validates its design arguments by name", {
  expect_error(simulate_mediation_data(n = 0, p = 20), "`n` must be")
  expect_error(simulate_mediation_data(n = 50, p = 1.5), "`p` must be")
  expect_error(simulate_mediation_data(n = 50, p = 20, c1 = "a"), "`c1` must be")
  expect_error(simulate_mediation_data(n = 50, p = 20, c2 = c(1, 2)), "`c2` must be")
  expect_error(simulate_mediation_data(n = 50, p = 20, sigma_y = 0), "`sigma_y` must be")
  expect_error(simulate_mediation_data(n = 50, p = 20, sigma_y = Inf), "`sigma_y` must be")
  expect_error(simulate_mediation_data(n = 50, p = 20, d = -1), "`d` must be")
  expect_error(simulate_mediation_data(n = 50, p = 20, q = 0), "`q` must be")
  expect_error(simulate_mediation_data(n = 50, p = 20, rho = 1), "`rho` must")
  expect_error(simulate_mediation_data(n = 50, p = 20, rho = NA), "`rho` must")
})
