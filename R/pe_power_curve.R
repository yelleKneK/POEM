#' Monte Carlo size and power curve for the PE mediation tests
#'
#' Reproduces the article's simulation studies: for each value of the
#' signal-strength scale \eqn{c_1} on a grid, it simulates many data sets,
#' applies both the benchmark Wald test and the power-enhanced test, and
#' returns their empirical rejection rates. At \eqn{c_1 = 0} (the global
#' null of no mediation) the rejection rate estimates the Type I error
#' rate; at nonzero \eqn{c_1} it estimates power. The headline finding is
#' visible directly in the table: under a contrasting pattern the
#' benchmark test stays near its nominal size as \eqn{|c_1|} grows while
#' the power-enhanced test climbs toward one.
#'
#' @details
#' A rejection rate from `n_rep` replications carries a Monte Carlo
#' standard error of about \eqn{\sqrt{r(1 - r) / n_{rep}}}; at 100
#' replications a rate near 0.05 is known to about 0.02 and a rate near
#' 0.5 to about 0.05. The article uses 1000 replications. A replication
#' whose penalized fit selects no mediator counts as a non-rejection (both
#' p-values are 1), the article's convention; the number of such
#' replications is reported in `n_empty`. A replication whose fit fails
#' outright is dropped from the denominator and reported in `n_valid`, with
#' a warning that quotes the first error.
#'
#' The tuning grid matters for what these curves show. The article's
#' figures were produced with per-setting grids (see [pe_lambda_grid()]);
#' the default here is that grid rescaled to `n` and `p`, so the
#' article's settings reproduce its curves within Monte Carlo error.
#' Pass `lambda_grid` to study another grid.
#'
#' @param n,p Number of observations and candidate mediators per data set.
#' @param outcome Outcome type passed to [simulate_mediation_data()]:
#'   `"continuous"`, `"binary"`, or `"count"`.
#' @param pattern Mediation pattern: `"homogeneous"` or `"contrasting"`
#'   (`"heterogeneous"` is a synonym for `"contrasting"`).
#' @param c1_grid Numeric vector of signal-strength scales to sweep.
#'   Default `seq(0, 1, by = 0.25)`. Include 0 to estimate the Type I
#'   error rate; the article also uses negative values.
#' @param c2 Direct effect used in the simulation. Default 0.5.
#' @param n_rep Number of Monte Carlo replications per grid point. Default
#'   100. The article uses 1000.
#' @param alpha_level Significance level for the global test. A replication
#'   counts as a rejection when its p-value is at most `alpha_level`. Default
#'   0.05.
#' @param method,error_level Passed to [pe_mediation()]: the multiplicity
#'   method and target error rate for the mediator-screening step.
#' @param lambda_grid,lambda_grid_reduced Passed to [pe_mediation()]: the
#'   tuning grids for the penalized fit. Default `NULL`, the family's
#'   default grid (see [pe_lambda_grid()]).
#' @param outcome_args A list of further arguments forwarded to
#'   [simulate_mediation_data()] (for example `rho`, `q`, `d`, `tau`,
#'   `alpha_m`). Default empty. The design arguments this function sets
#'   itself (`n`, `p`, `outcome`, `pattern`, `c1`, `c2`) and `seed` may not
#'   appear here; a `seed` inside `outcome_args` would make every
#'   replication draw the same data.
#' @param cores Number of CPU cores for the Monte Carlo replications.
#'   Values above 1 fork via the base \pkg{parallel} package (Unix only).
#'   With parallel cores, reproducibility uses parallel-safe streams, so a
#'   seeded parallel run is reproducible across runs at the same `cores` but
#'   need not match a serial run. Default 1.
#' @param seed Optional integer seed, set locally with the caller's random
#'   number generator state restored on exit.
#' @param progress Logical; if `TRUE`, print a line per grid point as it
#'   completes. Default `FALSE`.
#'
#' @return A tidy `data.frame` of class `poem_tbl` with one row per grid
#'   point and columns `c1`, `rejection_hdmm`, and `rejection_pe` (the
#'   empirical rejection rates of the benchmark and power-enhanced tests),
#'   `n_valid` (replications whose fit succeeded, the denominator of the
#'   rates), and `n_empty` (replications whose penalized fit selected no
#'   mediator, counted as non-rejections). The simulation settings are
#'   recorded in the attributes `outcome`, `pattern`, `n`, `p`, `n_rep`,
#'   `alpha_level`, `error_level`, and `lambda_grid`.
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @seealso [pe_mediation()], [simulate_mediation_data()],
#'   [pe_lambda_grid()].
#'
#' @family mediation simulation
#'
#' @examples
#' # n_rep = 5 keeps this example fast; a rate from five replications has a
#' # Monte Carlo standard error of up to 0.22, so read the shape, not the
#' # numbers. A reported study uses the article's design (n = 300, p = 500)
#' # and 1000 replications.
#' set.seed(113)
#' pe_power_curve(n = 200, p = 60, outcome = "continuous",
#'                pattern = "contrasting", c1_grid = c(0, 0.5, 1),
#'                n_rep = 5)
#'
#' @export
pe_power_curve <- function(n, p,
                           outcome = c("continuous", "binary", "count"),
                           pattern = c("homogeneous", "contrasting",
                                       "heterogeneous"),
                           c1_grid = seq(0, 1, by = 0.25), c2 = 0.5,
                           n_rep = 100, alpha_level = 0.05,
                           method = c("Bonferroni", "BH", "BY"),
                           error_level = 0.05, lambda_grid = NULL,
                           lambda_grid_reduced = NULL, outcome_args = list(),
                           cores = 1L, seed = NULL, progress = FALSE) {
  outcome <- match.arg(outcome)
  pattern <- match.arg(pattern)
  method <- match.arg(method)
  .validate_mc_design(n, p, n_rep, c1_grid, c2, cores)
  .validate_alpha_level(alpha_level, "alpha_level")
  .validate_error_level(error_level)
  .validate_outcome_args(outcome_args)
  if (!is.null(lambda_grid)) .validate_lambda_grid(lambda_grid, "lambda_grid")
  if (!is.null(lambda_grid_reduced))
    .validate_lambda_grid(lambda_grid_reduced, "lambda_grid_reduced")

  .poem_local_seed(seed, parallel = cores > 1L)

  reject_hdmm <- numeric(length(c1_grid))
  reject_pe   <- numeric(length(c1_grid))
  n_valid     <- integer(length(c1_grid))
  n_empty     <- integer(length(c1_grid))

  for (g in seq_along(c1_grid)) {
    c1 <- c1_grid[g]
    # One replication: simulate, fit both tests, return their two p-values
    # and whether the penalized selection was empty. An empty selection is a
    # non-rejection (both p-values 1) and is counted rather than warned about
    # here; a fit that errors returns NA with its message and is dropped from
    # the denominator rather than aborting the sweep.
    one_rep <- function(r) {
      dat <- do.call(simulate_mediation_data,
                     c(list(n = n, p = p, outcome = outcome, pattern = pattern,
                            c1 = c1, c2 = c2), outcome_args))
      .pe_mc_fit(dat, outcome = outcome, method = method,
                 error_level = error_level, lambda_grid = lambda_grid,
                 lambda_grid_reduced = lambda_grid_reduced)
    }
    reps <- .pe_lapply(seq_len(n_rep), one_rep, cores = cores)
    p_hdmm <- vapply(reps, function(x) x$p_hdmm, numeric(1))
    p_pe   <- vapply(reps, function(x) x$p_pe, numeric(1))
    empty  <- vapply(reps, function(x) x$empty, logical(1))
    errs   <- unlist(lapply(reps, function(x) x$error))
    ok <- !is.na(p_hdmm) & !is.na(p_pe)
    .pe_mc_report_failures(errs, n_rep, sprintf("c1 = %.3g", c1))
    n_valid[g]     <- sum(ok)
    n_empty[g]     <- sum(empty[ok])
    reject_hdmm[g] <- mean(p_hdmm[ok] <= alpha_level)
    reject_pe[g]   <- mean(p_pe[ok] <= alpha_level)
    if (progress)
      message(sprintf("c1 = %.3g: reject HDMM = %.3f, PE = %.3f (%d/%d valid, %d empty)",
                      c1, reject_hdmm[g], reject_pe[g], n_valid[g], n_rep,
                      n_empty[g]))
  }

  out <- data.frame(c1 = c1_grid, rejection_hdmm = reject_hdmm,
                    rejection_pe = reject_pe, n_valid = n_valid,
                    n_empty = n_empty, stringsAsFactors = FALSE)
  out <- .as_poem_tbl(out)
  attr(out, "outcome")     <- outcome
  attr(out, "pattern")     <- pattern
  attr(out, "n")           <- n
  attr(out, "p")           <- p
  attr(out, "n_rep")       <- n_rep
  attr(out, "alpha_level") <- alpha_level
  attr(out, "error_level") <- error_level
  attr(out, "lambda_grid") <- if (is.null(lambda_grid))
    .pe_default_grid(n, p, outcome, "full") else lambda_grid
  out
}

# One Monte Carlo fit shared by the study functions: fit pe_mediation() on a
# simulated data set, muffling the classed empty-fit warning (counted instead)
# and turning an error into an NA record that carries the message. Returns a
# list with the two p-values, the empty flag, the fit (or NULL), and the error
# message (or NULL). Not exported.
#' @keywords internal
#' @noRd
.pe_mc_fit <- function(dat, outcome, method, error_level, lambda_grid,
                       lambda_grid_reduced, ...) {
  state <- new.env(parent = emptyenv()); state$empty <- FALSE
  fit <- tryCatch(
    withCallingHandlers(
      pe_mediation(dat$X, dat$Y, dat$M, Z = dat$Z, outcome = outcome,
                   method = method, error_level = error_level,
                   lambda_grid = lambda_grid,
                   lambda_grid_reduced = lambda_grid_reduced, ...),
      poem_empty_fit = function(w) {
        state$empty <- TRUE
        invokeRestart("muffleWarning")
      }),
    error = function(e) e)
  if (inherits(fit, "error"))
    return(list(p_hdmm = NA_real_, p_pe = NA_real_, empty = FALSE, fit = NULL,
                error = conditionMessage(fit)))
  list(p_hdmm = fit$value[fit$term == "pval_hdmm"],
       p_pe = fit$value[fit$term == "pval_pe"], empty = state$empty, fit = fit,
       error = NULL)
}

# Report failed Monte Carlo replications: a warning naming the count and the
# first message when some failed, an error when every replication failed (a
# rate over zero valid replications is not a number). Not exported.
#' @keywords internal
#' @noRd
.pe_mc_report_failures <- function(errs, n_rep, where) {
  if (!length(errs)) return(invisible(TRUE))
  if (length(errs) >= n_rep)
    stop(sprintf("Every one of the %d replications at %s failed to fit; the first error was: %s",
                 n_rep, where, errs[1L]), call. = FALSE)
  warning(sprintf("%d of %d replications at %s failed to fit and were dropped; the first error was: %s",
                  length(errs), n_rep, where, errs[1L]), call. = FALSE)
  invisible(TRUE)
}

# Validation shared by the Monte Carlo study functions. Not exported.
#' @keywords internal
#' @noRd
.validate_mc_design <- function(n, p, n_rep, c1_grid, c2, cores) {
  for (nm in c("n", "p", "n_rep")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 1 ||
        v != round(v))
      stop(sprintf("`%s` must be a single positive whole number.", nm),
           call. = FALSE)
  }
  if (!is.numeric(c1_grid) || length(c1_grid) < 1L || any(!is.finite(c1_grid)))
    stop("`c1_grid` must be a non-empty numeric vector of finite values.",
         call. = FALSE)
  if (!is.numeric(c2) || length(c2) != 1L || !is.finite(c2))
    stop("`c2` must be a single finite number.", call. = FALSE)
  if (!is.numeric(cores) || length(cores) != 1L || !is.finite(cores) ||
      cores < 1 || cores != round(cores))
    stop("`cores` must be a single positive whole number.", call. = FALSE)
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
.validate_alpha_level <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x >= 1)
    stop(sprintf("`%s` must be a single number in (0, 1).", name), call. = FALSE)
  invisible(TRUE)
}

# `outcome_args` may carry only arguments the study does not set itself; a
# `seed` in it would hand every replication the same data set. Not exported.
#' @keywords internal
#' @noRd
.validate_outcome_args <- function(outcome_args) {
  if (!is.list(outcome_args))
    stop("`outcome_args` must be a list.", call. = FALSE)
  reserved <- c("n", "p", "outcome", "pattern", "c1", "c2", "seed")
  bad <- intersect(names(outcome_args), reserved)
  if (length(bad))
    stop("`outcome_args` may not contain ", paste0("`", bad, "`", collapse = ", "),
         "; the study sets those itself (and a `seed` there would make every ",
         "replication identical).", call. = FALSE)
  invisible(TRUE)
}
