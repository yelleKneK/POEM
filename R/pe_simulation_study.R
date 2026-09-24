#' Run the power study under several mediation patterns at once
#'
#' A convenience wrapper around [pe_power_curve()] that runs the size and
#' power study under more than one mediation pattern in a single call and
#' stacks the results, so a homogeneous and a contrasting curve (the two
#' panels of a figure in the article come back
#' together. Each pattern is swept over the same signal grid.
#'
#' @param n,p Observations and candidate mediators per simulated data set.
#' @param outcome Outcome type: `"continuous"`, `"binary"`, or `"count"`.
#' @param patterns Character vector of mediation patterns to run. Default
#'   `c("homogeneous", "contrasting")`.
#' @param c1_grid,c2,n_rep,alpha_level,method,error_level,lambda_grid,lambda_grid_reduced,cores,seed
#'   Passed to [pe_power_curve()].
#'
#' @return A tidy `data.frame` (class `poem_tbl`) with a leading `pattern`
#'   column and the [pe_power_curve()] columns (`c1`, `rejection_hdmm`,
#'   `rejection_pe`, `n_valid`, `n_empty`) for each pattern. The settings
#'   are recorded in the attributes `outcome`, `n`, `p`, `n_rep`, `alpha_level`,
#'   `error_level`, and `lambda_grid`.
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @seealso [pe_power_curve()] for a single pattern, [plot.poem_tbl()] to
#'   draw one pattern's curves.
#'
#' @family mediation simulation
#'
#' @examples
#' # n_rep = 5 keeps this example fast; a rate from five replications has a
#' # Monte Carlo standard error of up to 0.22, so read the shape, not the
#' # numbers. A reported study uses the article's design (n = 300, p = 500)
#' # and 1000 replications.
#' set.seed(113)
#' pe_simulation_study(n = 120, p = 50, outcome = "continuous",
#'                     c1_grid = c(0, 1), n_rep = 5)
#'
#' @export
pe_simulation_study <- function(n, p,
                                outcome = c("continuous", "binary", "count"),
                                patterns = c("homogeneous", "contrasting"),
                                c1_grid = seq(0, 1, by = 0.25), c2 = 0.5,
                                n_rep = 100, alpha_level = 0.05,
                                method = c("Bonferroni", "BH", "BY"),
                                error_level = 0.05, lambda_grid = NULL,
                                lambda_grid_reduced = NULL, cores = 1L,
                                seed = NULL) {
  outcome <- match.arg(outcome)
  method <- match.arg(method)
  patterns <- match.arg(patterns, c("homogeneous", "contrasting",
                                    "heterogeneous"), several.ok = TRUE)
  .validate_mc_design(n, p, n_rep, c1_grid, c2, cores)
  .validate_alpha_level(alpha_level, "alpha_level")
  .validate_error_level(error_level)

  # Seed once here (rather than per pattern) so the whole study is reproducible
  # and the patterns are not handed the same stream.
  .poem_local_seed(seed)

  parts <- lapply(patterns, function(pat) {
    pc <- pe_power_curve(n = n, p = p, outcome = outcome, pattern = pat,
                         c1_grid = c1_grid, c2 = c2, n_rep = n_rep,
                         alpha_level = alpha_level, method = method,
                         error_level = error_level, lambda_grid = lambda_grid,
                         lambda_grid_reduced = lambda_grid_reduced,
                         cores = cores)
    cbind(pattern = pat, as.data.frame(unclass(pc)), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out <- .as_poem_tbl(out)
  attr(out, "outcome")     <- outcome
  attr(out, "n")           <- n
  attr(out, "p")           <- p
  attr(out, "n_rep")       <- n_rep
  attr(out, "alpha_level") <- alpha_level
  attr(out, "error_level") <- error_level
  attr(out, "lambda_grid") <- if (is.null(lambda_grid))
    .pe_default_grid(n, p, outcome, "full") else lambda_grid
  out
}
