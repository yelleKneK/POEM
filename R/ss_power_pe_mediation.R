#' Simulation-based sample size planning for the power-enhanced mediation test
#'
#' Estimates, by Monte Carlo simulation, the sample size needed for the
#' power-enhanced test to reach a target power under a given mediation
#' pattern and signal strength. It evaluates the empirical power of the PE
#' test (and, for comparison, the benchmark Wald test) at each candidate
#' sample size on a grid, and reports the smallest grid value that reaches
#' the target. This is the design counterpart of the analysis function
#' [pe_mediation()]: use it to plan a study, then [pe_mediation()] to
#' analyze it.
#'
#' @details
#' This is a simulation-based procedure, not a closed-form power formula:
#' no analytic power expression exists for the PE test, so sample size is
#' instead determined by direct Monte Carlo simulation. Power is estimated
#' by simulating `n_rep` data sets at each candidate
#' `n` with [simulate_mediation_data()] and recording how often the test
#' rejects at level `alpha_level`. Because the estimate is a Monte Carlo
#' proportion, it carries simulation error of roughly
#' \eqn{\sqrt{power(1 - power) / n\_rep}}; raise `n_rep` for a smoother
#' curve and a more stable recommendation. The grid approach (rather than
#' a root search) is deliberate: the power curve is monotone but noisy, so
#' a search can stop early on a lucky draw, whereas the grid shows the
#' whole trajectory.
#'
#' @param outcome Outcome type: `"continuous"`, `"binary"`, or `"count"`.
#' @param pattern Mediation pattern: `"homogeneous"` or `"contrasting"`
#'   (`"heterogeneous"` is a synonym for `"contrasting"`).
#' @param p Number of candidate mediators.
#' @param n_grid Vector of candidate sample sizes to evaluate.
#' @param c1 Signal strength (the exposure-on-mediator scale). Default 1.
#' @param c2 Direct effect. Default 0.5.
#' @param target_power Desired power for the PE test. Default 0.8.
#' @param n_rep Monte Carlo replications per candidate `n`. Default 200.
#' @param alpha_level Significance level. Default 0.05.
#' @param method,error_level Passed to [pe_mediation()].
#' @param lambda_grid,lambda_grid_reduced Passed to [pe_mediation()]; the
#'   default `NULL` uses the family's default grid at each candidate `n`
#'   (for a continuous outcome the article's grid rescaled by
#'   [pe_lambda_grid()], so the grid follows the sample size being planned).
#' @param cores Number of CPU cores (Unix forking). A seeded parallel run is
#'   reproducible across runs at the same `cores` but need not match a
#'   serial run. Default 1.
#' @param seed Optional integer seed, set locally and restored on exit.
#'
#' @return A tidy `data.frame` (class `poem_tbl`) with one row per
#'   candidate sample size and columns `n`, `power_pe`, `power_hdmm`,
#'   `n_valid`, and `n_empty` (see [pe_power_curve()]). The smallest `n`
#'   reaching `target_power` for the PE test is stored in the
#'   `"recommended_n"` attribute (`NA` if no grid value reaches it), with
#'   `"target_power"`, `"outcome"`, `"pattern"`, `"c1"`, `"c2"`, `"n_rep"`,
#'   and `"alpha_level"`.
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @seealso [pe_power_curve()], [pe_mediation()],
#'   [simulate_mediation_data()].
#'
#' @family mediation simulation
#'
#' @examples
#' # n_rep = 5 keeps this example fast; a power estimate from five
#' # replications has a Monte Carlo standard error of up to 0.22, so the
#' # recommended n here is a demonstration, not a plan. Planning a study
#' # deserves n_rep of several hundred.
#' set.seed(113)
#' plan <- ss_power_pe_mediation(outcome = "continuous", pattern = "contrasting",
#'                               p = 50, n_grid = c(80, 120, 160), n_rep = 5)
#' plan
#' attr(plan, "recommended_n")
#'
#' @export
ss_power_pe_mediation <- function(outcome = c("continuous", "binary", "count"),
                                  pattern = c("homogeneous", "contrasting",
                                              "heterogeneous"),
                                  p, n_grid, c1 = 1, c2 = 0.5,
                                  target_power = 0.8, n_rep = 200, alpha_level = 0.05,
                                  method = c("Bonferroni", "BH", "BY"),
                                  error_level = 0.05, lambda_grid = NULL,
                                  lambda_grid_reduced = NULL, cores = 1L,
                                  seed = NULL) {
  outcome <- match.arg(outcome)
  pattern <- match.arg(pattern)
  method <- match.arg(method)
  if (!is.numeric(n_grid) || length(n_grid) < 1L || any(!is.finite(n_grid)) ||
      any(n_grid < 10))
    stop("`n_grid` must be a vector of sample sizes (each >= 10).",
         call. = FALSE)
  if (!is.numeric(c1) || length(c1) != 1L || !is.finite(c1))
    stop("`c1` must be a single finite number (one signal strength to plan for).",
         call. = FALSE)
  .validate_mc_design(n_grid[1L], p, n_rep, c1, c2, cores)
  .validate_alpha_level(alpha_level, "alpha_level")
  .validate_error_level(error_level)
  if (!is.numeric(target_power) || length(target_power) != 1L ||
      !is.finite(target_power) || target_power <= 0 || target_power >= 1)
    stop("`target_power` must be a single number in (0, 1).", call. = FALSE)

  .poem_local_seed(seed)

  # Evaluate power at each candidate n by reusing the single-point power curve
  # (c1 fixed); the c1 = 0 size case is not needed here, so the grid is one c1.
  n_grid <- sort(unique(as.integer(n_grid)))
  rows <- lapply(n_grid, function(nn) {
    pc <- pe_power_curve(n = nn, p = p, outcome = outcome, pattern = pattern,
                         c1_grid = c1, c2 = c2, n_rep = n_rep, alpha_level = alpha_level,
                         method = method, error_level = error_level,
                         lambda_grid = lambda_grid,
                         lambda_grid_reduced = lambda_grid_reduced,
                         cores = cores)
    data.frame(n = nn, power_pe = pc$rejection_pe[1],
               power_hdmm = pc$rejection_hdmm[1], n_valid = pc$n_valid[1],
               n_empty = pc$n_empty[1], stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)

  # Recommended n: the smallest grid value whose PE power reaches the target.
  reached <- which(out$power_pe >= target_power)
  recommended_n <- if (length(reached)) out$n[min(reached)] else NA_integer_

  out <- .as_poem_tbl(out)
  attr(out, "recommended_n") <- recommended_n
  attr(out, "target_power")  <- target_power
  attr(out, "outcome")       <- outcome
  attr(out, "pattern")       <- pattern
  attr(out, "c1")            <- c1
  attr(out, "c2")            <- c2
  attr(out, "n_rep")         <- n_rep
  attr(out, "alpha_level")   <- alpha_level
  out
}
