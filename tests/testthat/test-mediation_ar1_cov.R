test_that("the AR(1) covariance has the right structure", {
  S <- mediation_ar1_cov(5, rho = 0.5)
  expect_equal(dim(S), c(5L, 5L))
  expect_equal(diag(S), rep(1, 5))            # unit variances
  expect_equal(S[1, 2], 0.5)                  # adjacent: rho^1
  expect_equal(S[1, 3], 0.25)                 # two apart: rho^2
  expect_equal(S[2, 5], 0.5^3)                # three apart
  expect_true(isSymmetric(S))
  # A valid covariance is positive definite (all eigenvalues > 0).
  expect_true(all(eigen(S, only.values = TRUE)$values > 0))
})

test_that("the AR(1) covariance validates its arguments", {
  expect_error(mediation_ar1_cov(0), "positive integer")
  expect_error(mediation_ar1_cov(3.5), "positive integer")
  expect_error(mediation_ar1_cov(5, rho = 1), "in \\(-1, 1\\)")
})
