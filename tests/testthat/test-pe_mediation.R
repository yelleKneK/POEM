# Extract one value from a poem_tbl result by term name.
v <- function(tab, t) tab$value[tab$term == t]

test_that("the output has the documented schema and a numeric value column", {
  set.seed(113)
  d <- simulate_mediation_data(n = 120, p = 50, outcome = "continuous",
                               pattern = "homogeneous", c1 = 0.6)
  fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
  expect_s3_class(fit, "poem_tbl")
  expect_true(is.numeric(fit$value))
  expect_true(all(c("stat_hdmm", "pval_hdmm", "stat_pe", "j_pe", "pval_pe",
                    "total_indirect_effect", "n_active_mediators", "df",
                    "n_candidate_mediators", "n_observations") %in% fit$term))
  expect_equal(v(fit, "df"), 1)
  expect_equal(v(fit, "n_candidate_mediators"), 50)
  expect_equal(v(fit, "n_observations"), 120)
})

test_that("M_PE = S_n + J_m and the p-value is the chi-square upper tail", {
  set.seed(113)
  d <- simulate_mediation_data(n = 120, p = 50, outcome = "continuous",
                               pattern = "contrasting", c1 = 0.6)
  fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
  # The power-enhanced statistic is exactly the Wald statistic plus the PE
  # component; this is the defining identity of the method.
  expect_equal(v(fit, "stat_pe"), v(fit, "stat_hdmm") + v(fit, "j_pe"))
  expect_equal(v(fit, "pval_pe"),
               pchisq(v(fit, "stat_pe"), df = v(fit, "df"), lower.tail = FALSE))
  expect_equal(v(fit, "pval_hdmm"),
               pchisq(v(fit, "stat_hdmm"), df = v(fit, "df"), lower.tail = FALSE))
})

test_that("the front end dispatches identically to the workers", {
  set.seed(113)
  d <- simulate_mediation_data(n = 120, p = 50, outcome = "continuous",
                               pattern = "homogeneous", c1 = 0.6)
  a <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
  b <- pe_mediation_linear(d$X, d$Y, d$M)
  expect_identical(a$value, b$value)
  expect_identical(attr(a, "active_mediators"), attr(b, "active_mediators"))
})

test_that("the PE test is never less significant than the Wald test (deterministic)", {
  set.seed(113)
  d <- simulate_mediation_data(n = 120, p = 50, outcome = "continuous",
                               pattern = "homogeneous", c1 = 0.6)
  fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
  # Because M_PE = S_n + J_m with J_m >= 0, the PE p-value is always at most
  # the Wald p-value; this holds for every data set, not just on average.
  expect_lte(v(fit, "pval_pe"), v(fit, "pval_hdmm"))
})

test_that("power enhancement detects contrasting mediation that the Wald test misses", {
  set.seed(113)
  # A strong contrasting signal where the total indirect effect is exactly
  # zero, so the benchmark test has essentially no power but the PE test does.
  d <- simulate_mediation_data(n = 200, p = 60, outcome = "continuous",
                               pattern = "contrasting", c1 = 1.0)
  expect_lt(abs(d$beta), 1e-8)
  fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
  expect_gt(v(fit, "pval_hdmm"), 0.05)        # Wald test does not reject
  expect_lt(v(fit, "pval_pe"), 0.01)          # PE test rejects decisively
  expect_gt(v(fit, "n_active_mediators"), 0)  # and names active mediators
})

test_that("j_pe is nonnegative and stat_pe is at least stat_hdmm", {
  set.seed(113)
  d <- simulate_mediation_data(n = 120, p = 50, outcome = "count",
                               pattern = "homogeneous", c1 = 0.4, c2 = 0.4)
  fit <- pe_mediation(d$X, d$Y, d$M, outcome = "count")
  expect_gte(v(fit, "j_pe"), 0)
  expect_gte(v(fit, "stat_pe"), v(fit, "stat_hdmm"))
})

test_that("input validation rejects malformed arguments", {
  set.seed(113)
  d <- simulate_mediation_data(n = 80, p = 40, outcome = "continuous",
                               pattern = "homogeneous")
  expect_error(pe_mediation(d$X, d$Y[-1], d$M, outcome = "continuous"),
               "same number of rows")
  expect_error(pe_mediation(d$X, d$Y, d$M, outcome = "continuous",
                            error_level = 0), "in \\(0, 1\\)")
  expect_error(pe_mediation(d$X, d$Y, d$M, outcome = "continuous",
                            error_level = 1.5), "in \\(0, 1\\)")
  expect_error(pe_mediation_linear(d$X, d$Y, d$M, lambda_grid = c(-1, 1)),
               "> 0")
  # A binary outcome routed to the linear worker should be rejected.
  expect_error(pe_mediation(d$X, rbinom(length(d$Y), 1, 0.5), d$M,
                            outcome = "continuous"), "binary")
})

test_that("the three multiplicity methods all run and preserve the identity", {
  set.seed(113)
  d <- simulate_mediation_data(n = 150, p = 60, outcome = "continuous",
                               pattern = "contrasting", c1 = 0.6)
  for (m in c("Bonferroni", "BH", "BY")) {
    fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous", method = m)
    expect_equal(v(fit, "stat_pe"), v(fit, "stat_hdmm") + v(fit, "j_pe"),
                 info = m)
  }
})

test_that("an empty penalized selection is audible and keeps the row schema", {
  set.seed(113)
  d <- simulate_mediation_data(n = 100, p = 30, outcome = "continuous",
                               pattern = "homogeneous", c1 = 0)
  # A grid of large penalties selects nothing at every value.
  expect_warning(
    fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous",
                        lambda_grid = c(5, 10)),
    class = "poem_empty_fit")
  expect_true(attr(fit, "empty_fit"))
  expect_equal(nrow(fit), 12L)                    # same rows as a non-empty fit
  expect_equal(v(fit, "stat_hdmm"), 0)
  expect_equal(v(fit, "pval_hdmm"), 1)
  expect_equal(v(fit, "stat_pe"), 0)
  expect_equal(v(fit, "pval_pe"), 1)
  expect_equal(v(fit, "total_indirect_effect"), 0)
  expect_true(is.na(v(fit, "total_indirect_lower")))
  expect_true(is.na(v(fit, "total_indirect_upper")))
  expect_length(attr(fit, "active_mediators"), 0L)
  expect_equal(nrow(pe_mediators(fit)), 0L)
  expect_output(print(fit), "Penalized selection: empty")
  # The same path in the two other families, with the per-method bookkeeping.
  db <- simulate_mediation_data(n = 100, p = 30, outcome = "binary",
                                pattern = "homogeneous", c1 = 0, c2 = 1,
                                seed = 113)
  expect_warning(
    fb <- pe_mediation(db$X, db$Y, db$M, outcome = "binary",
                       lambda_grid = c(5, 10), report_all_methods = TRUE),
    class = "poem_empty_fit")
  expect_true(attr(fb, "empty_fit"))
  expect_equal(attr(fb, "pe_by_method")$pval_pe, c(1, 1, 1))
  expect_true(all(lengths(attr(fb, "selection_by_method")) == 0L))
  dc <- simulate_mediation_data(n = 100, p = 30, outcome = "count",
                                pattern = "homogeneous", c1 = 0, c2 = 0.4,
                                seed = 113)
  expect_warning(
    fc <- pe_mediation(dc$X, dc$Y, dc$M, outcome = "count",
                       lambda_grid = c(50, 100)),
    class = "poem_empty_fit")
  expect_true(attr(fc, "empty_fit"))
})

test_that("a non-empty fit is not flagged empty and carries a finite interval", {
  d <- simulate_mediation_data(n = 150, p = 40, outcome = "continuous",
                               pattern = "homogeneous", c1 = 1, seed = 113)
  fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
  expect_false(attr(fit, "empty_fit"))
  expect_true(is.finite(v(fit, "total_indirect_lower")))
  expect_lt(v(fit, "total_indirect_lower"), v(fit, "total_indirect_effect"))
  expect_gt(v(fit, "total_indirect_upper"), v(fit, "total_indirect_effect"))
})

test_that("missing, non-finite, constant, and single-column inputs stop with the argument named", {
  d <- simulate_mediation_data(n = 80, p = 20, outcome = "continuous",
                               pattern = "homogeneous", c1 = 1, d = 2, seed = 113)
  Xna <- d$X; Xna[3, 1] <- NA
  expect_error(pe_mediation(Xna, d$Y, d$M, outcome = "continuous"), "`X` contains 1 missing")
  Yna <- d$Y; Yna[c(2, 5)] <- c(NaN, Inf)
  expect_error(pe_mediation(d$X, Yna, d$M, outcome = "continuous"), "`Y` contains 2 missing")
  Mna <- d$M; Mna[7, 4] <- NA
  expect_error(pe_mediation(d$X, d$Y, Mna, outcome = "continuous"), "`M` contains 1 missing")
  Zna <- d$Z; Zna[1, 1] <- NA
  expect_error(pe_mediation(d$X, d$Y, d$M, Z = Zna, outcome = "continuous"), "`Z` contains 1 missing")
  # the same guard reaches the binary and count workers
  db <- simulate_mediation_data(n = 80, p = 20, outcome = "binary", c1 = 1, c2 = 1, seed = 113)
  Yb <- db$Y; Yb[1] <- NA
  expect_error(pe_mediation(db$X, Yb, db$M, outcome = "binary"), "`Y` contains 1 missing")
  # constant exposure or confounder
  expect_error(pe_mediation(rep(1, 80), d$Y, d$M, outcome = "continuous"), "`X` has constant column")
  expect_error(pe_mediation(d$X, d$Y, d$M, Z = cbind(d$Z, 2), outcome = "continuous"), "`Z` has constant column")
  # a single mediator
  expect_error(pe_mediation(d$X, d$Y, d$M[, 1, drop = FALSE], outcome = "continuous"), "at least two columns")
})

test_that("column names on M label the active set, the table, and the footer", {
  set.seed(113)
  d <- simulate_mediation_data(n = 200, p = 60, outcome = "continuous",
                               pattern = "contrasting", c1 = 1)
  M <- d$M; colnames(M) <- paste0("med", seq_len(ncol(M)))
  fit <- pe_mediation(d$X, d$Y, M, outcome = "continuous")
  plain <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
  # Names change no number and no position.
  expect_equal(fit$value, plain$value)
  expect_identical(attr(fit, "active_mediators"), attr(plain, "active_mediators"))
  expect_identical(attr(fit, "mediator_names"), colnames(M))
  expect_null(attr(plain, "mediator_names"))
  tab <- pe_mediators(fit)
  expect_identical(names(tab)[1:2], c("mediator", "name"))
  expect_identical(tab$name, colnames(M)[tab$mediator])
  expect_false("name" %in% names(pe_mediators(plain)))
  active <- attr(fit, "active_mediators")
  expect_gt(length(active), 0L)
  expect_output(print(fit), paste(colnames(M)[active], collapse = ", "),
                fixed = TRUE)
  expect_match(pe_selection(fit)$active_mediators, colnames(M)[active][1],
               fixed = TRUE)
  expect_true("name" %in% names(tidy(fit)))
  # A data.frame M carries its names too; blank or duplicated names do not.
  fit_df <- pe_mediation(d$X, d$Y, as.data.frame(M), outcome = "continuous")
  expect_identical(attr(fit_df, "mediator_names"), colnames(M))
  dup <- M; colnames(dup) <- rep("m", ncol(M))
  expect_null(attr(pe_mediation(d$X, d$Y, dup, outcome = "continuous"),
                   "mediator_names"))
  # With a dropped constant column the names still index the original M.
  M2 <- M; M2[, 3] <- 1
  fit2 <- suppressWarnings(pe_mediation(d$X, d$Y, M2, outcome = "continuous",
                                        drop_constant = TRUE))
  tab2 <- pe_mediators(fit2)
  expect_identical(tab2$name, colnames(M2)[tab2$mediator])
  expect_false("med3" %in% tab2$name)
})

test_that("named exposures name the per-exposure rows", {
  set.seed(113)
  n <- 150; p <- 30
  tau <- rbind(c(0.1, 0.2, 0.3, 0.4, 0.5, rep(0, p - 5)),
               c(0.5, 0.4, 0.3, 0.2, 0.1, rep(0, p - 5)))
  d <- simulate_mediation_data(n = n, p = p, outcome = "continuous",
                               pattern = "contrasting", c1 = 1, q = 2, tau = tau)
  X <- d$X; colnames(X) <- c("dose", "age")
  fit <- pe_mediation(X, d$Y, d$M, outcome = "continuous")
  expect_true(all(c("total_indirect_effect_dose", "total_indirect_effect_age",
                    "total_indirect_lower_dose", "total_indirect_upper_age")
                  %in% fit$term))
  expect_identical(attr(fit, "exposure_names"), c("dose", "age"))
  plain <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
  expect_true(all(c("total_indirect_effect_1", "total_indirect_effect_2")
                  %in% plain$term))
  expect_null(attr(plain, "exposure_names"))
  expect_equal(fit$value, plain$value)
  g <- glance(fit)
  expect_true(all(c("total_indirect_effect_dose", "total_indirect_effect_age")
                  %in% names(g)))
  expect_equal(unname(unlist(g[c("total_indirect_effect_dose", "total_indirect_effect_age")])),
               unname(unlist(glance(plain)[c("total_indirect_effect_1", "total_indirect_effect_2")])))
})

test_that("the print and format methods validate their digits arguments", {
  set.seed(113)
  d <- simulate_mediation_data(n = 100, p = 20, outcome = "continuous",
                               pattern = "contrasting", c1 = 1)
  fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
  expect_error(print(fit, digits = 0), "`digits` must be")
  expect_error(format(fit, digits_p = 2.5), "`digits_p` must be")
  expect_error(print(fit, digits = "4"), "`digits` must be")
  expect_output(print(fit, digits = 6, digits_p = 6), "Outcome model")
})
