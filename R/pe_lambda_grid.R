#' Default tuning-parameter grid for the penalized mediator fit
#'
#' Builds the grid of candidate SCAD tuning parameters that
#' [pe_mediation()] searches by the high-dimensional BIC when no
#' `lambda_grid` is supplied. The grid is the one the article's simulation
#' scripts used at their design (\eqn{n = 300}, \eqn{p = 500}), rescaled to
#' the sample size and number of mediators at hand by the rate the theory
#' requires of the tuning parameter, \eqn{\sqrt{\log p / n}}.
#'
#' @details
#' The SCAD penalty in Yu and Kelley (in press) has a single tuning
#' parameter \eqn{\lambda}; the fit is repeated at each value of a candidate
#' grid and the value minimizing the high-dimensional BIC (HBIC) is kept.
#' Two facts make the grid, and especially its lower end, part of the
#' estimator rather than a detail. First, the article's theory requires
#' \eqn{\lambda} to be of larger order than \eqn{\sqrt{\log p / n}}
#' (and than \eqn{\sqrt{s_0 / n}}, with \eqn{s_0} the number of active
#' mediators); a grid whose values fall below that regime lets HBIC admit
#' spurious mediators to the penalized selection, and the screening step
#' then reports some of them as active. Second, the HBIC minimum drifts
#' toward small \eqn{\lambda} and often sits at the smallest value offered,
#' so the grid's lower end frequently decides the fit.
#'
#' The article's scripts set the grid by hand for each outcome family at
#' \eqn{n = 300}, \eqn{p = 500}: for the continuous outcome, 20 values from
#' 0.2 to 0.39 (full model) and from 0.27 to 0.46 (reduced model); for the
#' binary outcome, 15 values from 0.04 to 0.2; for the count outcome, 15
#' values from 0.7 to 5. At that design \eqn{\sqrt{\log p / n} = 0.1439},
#' so the continuous grid's lower end is about 1.4 times the rate. This
#' function returns exactly those grids at the article's design and scales
#' their endpoints by \eqn{\sqrt{\log p / n}} elsewhere, which keeps the
#' lower end in the regime the theory requires as \eqn{n} and \eqn{p}
#' change. The multipliers differ across families because the penalized
#' likelihood of a logistic or Poisson outcome is on a different scale from
#' penalized least squares.
#'
#' Which grid [pe_mediation()] searches by default depends on the outcome.
#' For a continuous outcome the default is this grid. The review-era
#' implementation searched a fixed grid of 20 values from 0.05 to 1 for
#' every family and design; at the article's linear design that grid lets
#' the identified set's familywise error rate rise from the reported 0.000
#' to about 0.17 under contrasting mediation (0.35 to 0.39 with a finer grid
#' over the same low range), and the article's linear tables do not
#' reproduce under it. The guarantee has a price at small designs, where
#' the rescaled floor removes more signal: at \eqn{n = 100}, \eqn{p = 40}
#' under contrasting mediation with \eqn{c_1 = 1} the power-enhanced test
#' rejects in 17 percent of replications with this grid against 77 percent
#' with the review-era grid, while the identified set's familywise error
#' rate is 0.000 against 0.07; at \eqn{n = 200}, \eqn{p = 100} the figures
#' are 58 against 97 percent and 0.000 against 0.13 (200 replications
#' each). Pass `lambda_grid = seq(0.05, 1, length.out = 20)` to trade the
#' guarantee for that power, or to reproduce results computed with the
#' review-era grid; the print footer reports which value HBIC chose.
#'
#' For binary and count outcomes the default of [pe_mediation()] remains
#' the review-era grid. The article's scripts for those families carry one
#' grid each (0.04 to 0.2 and 0.7 to 5), which reproduce the article's
#' homogeneous panels, while its contrasting panels reproduce under the
#' review-era grid and not under the family grids (the Poisson grid leaves
#' about a third of contrasting-count fits with no selected mediator); the
#' grids behind each panel were set by hand and not recorded. This function
#' still returns the family scripts' grids for those outcomes, for use with
#' the homogeneous designs, and the default will follow the article's grids
#' once they are on record.
#'
#' @param n,p Number of observations and candidate mediators.
#' @param outcome Outcome type: `"continuous"`, `"binary"`, or `"count"`.
#' @param model Which fit the grid is for: `"full"` (the mediator model
#'   used for selection and the power-enhanced screen) or `"reduced"` (the
#'   reduced model the continuous-outcome Wald test refits; for binary and
#'   count outcomes the reduced grid is not used and `"reduced"` returns
#'   the same grid as `"full"`).
#' @param length.out Number of grid values. Default `NULL` uses the
#'   article's counts: 20 for a continuous outcome, 15 for binary and
#'   count outcomes.
#'
#' @return A numeric vector of increasing candidate tuning parameters.
#'
#' @references
#' Yu, X., and Kelley, K. (in press). Power Enhancement in
#' High-Dimensional Heterogeneous Mediation Analysis. \emph{Journal of the
#' American Statistical Association}.
#'
#' Wang, L., Kim, Y., and Li, R. (2013). Calibrating nonconvex penalized
#' regression in ultra-high dimension. \emph{The Annals of Statistics,
#' 41}(5), 2505--2536. \doi{10.1214/13-AOS1159}
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @seealso [pe_mediation()], whose `lambda_grid` and
#'   `lambda_grid_reduced` arguments default to this grid.
#'
#' @family mediation tests
#'
#' @examples
#' # The article's linear-model grid at its own design
#' pe_lambda_grid(n = 300, p = 500, outcome = "continuous")
#' # The same grid rescaled to a smaller study
#' pe_lambda_grid(n = 100, p = 200, outcome = "continuous")
#' pe_lambda_grid(n = 300, p = 500, outcome = "count")
#'
#' @export
pe_lambda_grid <- function(n, p, outcome = c("continuous", "binary", "count"),
                           model = c("full", "reduced"), length.out = NULL) {
  outcome <- match.arg(outcome)
  model <- match.arg(model)
  for (nm in c("n", "p")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 2)
      stop(sprintf("`%s` must be a single number of at least 2.", nm),
           call. = FALSE)
  }
  if (is.null(length.out))
    length.out <- if (outcome == "continuous") 20L else 15L
  if (!is.numeric(length.out) || length(length.out) != 1L ||
      !is.finite(length.out) || length.out < 2 || length.out != round(length.out))
    stop("`length.out` must be a single integer of at least 2.", call. = FALSE)

  # The article's grids at its simulation design (n = 300, p = 500), and the
  # rate sqrt(log p / n) at that design; the endpoints scale with the rate.
  ref_rate <- sqrt(log(500) / 300)
  ends <- switch(outcome,
                 continuous = if (model == "full") c(0.2, 0.39) else c(0.27, 0.46),
                 binary     = c(0.04, 0.2),
                 count      = c(0.7, 5))
  rate <- sqrt(log(p) / n)
  seq(ends[1L] / ref_rate * rate, ends[2L] / ref_rate * rate,
      length.out = length.out)
}
