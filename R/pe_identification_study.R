#' Monte Carlo study of individual-mediator identification (FWER / FDR)
#'
#' Evaluates how well the power-enhanced screen recovers the *individual*
#' active mediators, the experiment reported in the article's supplement.
#' For each signal-strength scale \eqn{c_1} on a grid it simulates many data
#' sets, identifies the active set under each multiplicity method, and scores
#' those selections against the known truth, returning the empirical
#' familywise error rate, false discovery rate, precision, and recall. This
#' complements [pe_power_curve()], which evaluates the *global* test rather
#' than which mediators are flagged.
#'
#' @details
#' The four metrics follow the article's supplement. In each replication
#' the selected set is compared with the truth set; the familywise error
#' rate is the proportion of replications selecting at least one mediator
#' outside the truth set, the false discovery rate is the mean false
#' discovery proportion (zero when nothing is selected), precision is the
#' mean proportion of selected mediators that are true (scored zero when
#' nothing is selected), and recall is the mean proportion of true mediators
#' selected. Each is a Monte Carlo proportion with standard error about
#' \eqn{\sqrt{r(1 - r) / n_{rep}}}; the article uses 1000 replications.
#'
#' What counts as a true mediator is set by `truth`. The article scores a
#' mediator as active when its outcome coefficient \eqn{\alpha_{m,j}} is
#' nonzero (`truth = "outcome_effect"`, the default), which is what its
#' supplement tables use and what makes precision and recall defined at
#' \eqn{c_1 = 0}, where the exposure-on-mediator paths are all zero; there
#' the familywise error rate is the probability of selecting a mediator
#' with no outcome effect. Under `truth = "mediation"` a mediator is active
#' only when both paths are nonzero (the simulator's `active_mediators`), so
#' at \eqn{c_1 = 0} no mediator is active, any selection is a false
#' positive, and recall is `NA`. For nonzero \eqn{c_1} the two definitions
#' agree under the article's designs, whose exposure-on-mediator
#' loadings are all nonzero.
#'
#' The tuning grid governs these rates as much as the screen does (see
#' [pe_lambda_grid()]); the default reproduces the article's settings.
#'
#' @param n,p Number of observations and candidate mediators per data set.
#' @param outcome Outcome type passed to [simulate_mediation_data()]:
#'   `"continuous"`, `"binary"`, or `"count"`.
#' @param pattern Mediation pattern: `"homogeneous"` or `"contrasting"`
#'   (`"heterogeneous"` is a synonym for `"contrasting"`).
#' @param c1_grid Numeric vector of signal-strength scales to sweep.
#'   Default `seq(0, 1, by = 0.25)`. Include 0 to estimate the null
#'   familywise error rate.
#' @param c2 Direct effect used in the simulation. Default 0.5.
#' @param n_rep Number of Monte Carlo replications per grid point. Default
#'   200. The article uses 1000.
#' @param methods Multiplicity methods to evaluate, any of `"Bonferroni"`
#'   (familywise error rate), `"BH"`, and `"BY"` (false discovery rate).
#'   Default all three.
#' @param error_level Target error rate for the mediator-screening step
#'   (the FWER level for Bonferroni, the FDR level for BH and BY). Default
#'   0.05.
#' @param truth Which mediators count as truly active: `"outcome_effect"`
#'   (nonzero outcome coefficient, the article's convention) or
#'   `"mediation"` (both paths nonzero). See Details.
#' @param lambda_grid,lambda_grid_reduced Passed to [pe_mediation()]: the
#'   tuning grids for the penalized fit. Default `NULL`, the family's
#'   default grid (see [pe_lambda_grid()]).
#' @param outcome_args A list of further arguments forwarded to
#'   [simulate_mediation_data()] (for example `rho`, `q`, `d`, `tau`,
#'   `alpha_m`). Default empty. The design arguments this function sets
#'   itself (`n`, `p`, `outcome`, `pattern`, `c1`, `c2`) and `seed` may not
#'   appear here.
#' @param cores Number of CPU cores for the replications. Values above 1
#'   fork via the base \pkg{parallel} package (Unix only); a seeded parallel
#'   run is reproducible across runs at the same `cores` but need not match a
#'   serial run. Default 1.
#' @param seed Optional integer seed, set locally with the caller's random
#'   number generator state restored on exit.
#' @param progress Logical; if `TRUE`, print a line per grid point. Default
#'   `FALSE`.
#'
#' @return A tidy `data.frame` of class `poem_tbl` with one row per
#'   (\eqn{c_1}, method) combination and columns `c1`, `method`, `fwer`,
#'   `fdr`, `precision`, `recall` (as defined in Details), `n_valid`
#'   (replications whose fit succeeded), and `n_empty` (replications whose
#'   penalized fit selected no mediator). The settings are recorded in the
#'   attributes `outcome`, `pattern`, `n`, `p`, `n_rep`, `error_level`,
#'   `truth`, and `lambda_grid`.
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @references
#' Yu, X., and Kelley, K. (in press). Power Enhancement in
#' High-Dimensional Heterogeneous Mediation Analysis. \emph{Journal of the
#' American Statistical Association}.
#'
#' @seealso [pe_power_curve()] for the global test, [pe_selection()] for the
#'   per-method active set of a single fit, [pe_lambda_grid()].
#'
#' @family mediation simulation
#'
#' @examples
#' # n_rep = 5 keeps this example fast; a rate from five replications has a
#' # Monte Carlo standard error of up to 0.22, so read the shape, not the
#' # numbers. A reported study uses the article's design (n = 300, p = 500)
#' # and 1000 replications.
#' set.seed(113)
#' pe_identification_study(n = 200, p = 60, outcome = "continuous",
#'                         pattern = "contrasting", c1_grid = c(0, 1),
#'                         n_rep = 5)
#'
#' @export
pe_identification_study <- function(n, p,
                                    outcome = c("continuous", "binary",
                                                "count"),
                                    pattern = c("homogeneous", "contrasting",
                                                "heterogeneous"),
                                    c1_grid = seq(0, 1, by = 0.25), c2 = 0.5,
                                    n_rep = 200,
                                    methods = c("Bonferroni", "BH", "BY"),
                                    error_level = 0.05,
                                    truth = c("outcome_effect", "mediation"),
                                    lambda_grid = NULL,
                                    lambda_grid_reduced = NULL,
                                    outcome_args = list(),
                                    cores = 1L, seed = NULL, progress = FALSE) {
  outcome <- match.arg(outcome)
  pattern <- match.arg(pattern)
  methods <- match.arg(methods, c("Bonferroni", "BH", "BY"), several.ok = TRUE)
  truth <- match.arg(truth)
  .validate_mc_design(n, p, n_rep, c1_grid, c2, cores)
  .validate_error_level(error_level)
  .validate_outcome_args(outcome_args)
  if (!is.null(lambda_grid)) .validate_lambda_grid(lambda_grid, "lambda_grid")
  if (!is.null(lambda_grid_reduced))
    .validate_lambda_grid(lambda_grid_reduced, "lambda_grid_reduced")

  .poem_local_seed(seed, parallel = cores > 1L)

  # The active set selected under each method, for one fit. report_all_methods
  # populates selection_by_method; the empty-fit path leaves it NULL (no
  # mediator selected under any method).
  method_sets <- function(fit) {
    sbm <- attr(fit, "selection_by_method")
    if (is.null(sbm)) {
      act <- attr(fit, "active_mediators")
      sbm <- stats::setNames(rep(list(integer(0)), length(methods)), methods)
      prim <- attr(fit, "method")
      if (prim %in% methods) sbm[[prim]] <- act
    }
    sbm[methods]
  }
  na_row <- function() data.frame(method = methods, any_fp = NA_real_,
                                  fdp = NA_real_, precision = NA_real_,
                                  recall = NA_real_, valid = FALSE,
                                  empty = FALSE, stringsAsFactors = FALSE)

  res <- vector("list", length(c1_grid))
  for (g in seq_along(c1_grid)) {
    c1 <- c1_grid[g]
    # One replication: simulate, fit with all requested methods, and score each
    # method's selection against the truth set. A failed fit returns a row of
    # NAs flagged invalid and is dropped from the denominators; an empty
    # penalized selection is scored (nothing selected) and counted.
    one_rep <- function(r) {
      dat <- do.call(simulate_mediation_data,
                     c(list(n = n, p = p, outcome = outcome, pattern = pattern,
                            c1 = c1, c2 = c2), outcome_args))
      true_set <- if (truth == "outcome_effect") which(dat$alpha_m != 0)
                  else dat$active_mediators
      mc <- .pe_mc_fit(dat, outcome = outcome, method = methods[1L],
                       error_level = error_level, lambda_grid = lambda_grid,
                       lambda_grid_reduced = lambda_grid_reduced,
                       report_all_methods = length(methods) > 1L)
      if (is.null(mc$fit)) {
        out <- na_row(); attr(out, "error") <- mc$error; return(out)
      }
      sets <- method_sets(mc$fit)
      do.call(rbind, lapply(methods, function(m) {
        sel <- sets[[m]]
        tp <- length(intersect(sel, true_set))
        fp <- length(setdiff(sel, true_set))
        n_sel <- length(sel); n_tru <- length(true_set)
        data.frame(
          method = m,
          any_fp = as.numeric(fp >= 1L),
          fdp = if (n_sel > 0L) fp / n_sel else 0,
          precision = if (n_sel > 0L) tp / n_sel else 0,
          recall = if (n_tru > 0L) tp / n_tru else NA_real_,
          valid = TRUE, empty = mc$empty, stringsAsFactors = FALSE)
      }))
    }
    rep_list <- .pe_lapply(seq_len(n_rep), one_rep, cores = cores)
    errs <- unlist(lapply(rep_list, function(x) attr(x, "error")))
    .pe_mc_report_failures(errs, n_rep, sprintf("c1 = %.3g", c1))
    reps <- do.call(rbind, rep_list)

    # Aggregate over the valid replications, per method.
    avg <- function(x) if (any(!is.na(x))) mean(x, na.rm = TRUE) else NA_real_
    rows <- lapply(methods, function(m) {
      sub <- reps[reps$method == m & reps$valid, , drop = FALSE]
      data.frame(c1 = c1, method = m,
                 fwer = avg(sub$any_fp), fdr = avg(sub$fdp),
                 precision = avg(sub$precision), recall = avg(sub$recall),
                 n_valid = nrow(sub), n_empty = sum(sub$empty),
                 stringsAsFactors = FALSE)
    })
    res[[g]] <- do.call(rbind, rows)
    if (progress)
      message(sprintf("c1 = %.3g done (%d methods x %d reps)",
                      c1, length(methods), n_rep))
  }

  out <- do.call(rbind, res)
  rownames(out) <- NULL
  out <- .as_poem_tbl(out)
  attr(out, "outcome")     <- outcome
  attr(out, "pattern")     <- pattern
  attr(out, "n")           <- n
  attr(out, "p")           <- p
  attr(out, "n_rep")       <- n_rep
  attr(out, "error_level") <- error_level
  attr(out, "truth")       <- truth
  attr(out, "lambda_grid") <- if (is.null(lambda_grid))
    .pe_default_grid(n, p, outcome, "full") else lambda_grid
  out
}
