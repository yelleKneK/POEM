#' Per-mediator detail behind a power-enhanced mediation test
#'
#' Returns the per-mediator table that explains a [pe_mediation()] result:
#' one row for each mediator the penalized fit selected as a candidate,
#' with its two path statistics, the screening p-value, and whether it
#' was identified as active. This is the "why" behind the active set, so a
#' user can see which mediators were close to the screening threshold and
#' which were screened out.
#'
#' @param fit A result from [pe_mediation()] or one of its workers.
#'
#' @return A tidy `data.frame` (class `poem_tbl`) with columns:
#' \describe{
#'   \item{mediator}{Column position of the mediator in the original `M`.}
#'   \item{name}{The mediator's column name in `M`; present only when `M`
#'     has column names.}
#'   \item{t_outcome}{The standardized mediator-on-outcome statistic
#'     (the \eqn{M \to Y} path, \eqn{\hat\alpha_{m,j} /
#'     \hat\sigma_{m,j}}).}
#'   \item{t_exposure}{The strongest standardized exposure-on-mediator
#'     statistic across exposures (the \eqn{X \to M} path).}
#'   \item{screen_p}{The screening p-value, the larger of the two path
#'     p-values (the smaller across exposures when there are several).
#'     A mediator is active when this clears the multiplicity threshold.}
#'   \item{selected}{Logical; whether the mediator was identified as
#'     active.}
#' }
#'   The table is empty when the fit selected no candidate mediators.
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @seealso [pe_mediation()], [pe_selection()].
#'
#' @family mediation tests
#'
#' @examples
#' set.seed(113)
#' d <- simulate_mediation_data(n = 200, p = 60, outcome = "continuous",
#'                              pattern = "contrasting", c1 = 1)
#' fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
#' pe_mediators(fit)
#'
#' @export
pe_mediators <- function(fit) {
  tab <- attr(fit, "mediator_table")
  if (is.null(tab))
    stop("`fit` does not carry a mediator table; was it produced by ",
         "pe_mediation()?", call. = FALSE)
  .as_poem_tbl(tab)
}

#' Active mediators identified under each multiplicity method
#'
#' Summarizes which mediators a [pe_mediation()] fit identifies as active.
#' When the fit was produced with `report_all_methods = TRUE`, this
#' reports the active set under each of the three multiplicity methods
#' (Bonferroni for familywise error rate control, Benjamini-Hochberg and
#' Benjamini-Yekutieli for false discovery rate control), so their
#' conservativeness can be compared on the same fit. Otherwise it reports
#' the single method that was used.
#'
#' @param fit A result from [pe_mediation()] or one of its workers.
#'
#' @return A tidy `data.frame` of class `poem_tbl` (so its p-values print
#'   to fixed decimals) with one row per multiplicity method and columns
#'   `method`, `n_active` (number of active mediators), `j_pe` (the power
#'   enhancement component \eqn{J_m} under that method's screen), `stat_pe`
#'   and `pval_pe` (the resulting power-enhanced statistic \eqn{M_{PE}} and
#'   its p-value), and `active_mediators` (their column names when `M` has
#'   them, otherwise their column positions, comma-separated, or
#'   `"none"`). Because each multiplicity method screens
#'   a different active set into \eqn{J_m}, the `pval_pe` column gives the
#'   PE / PE_BH / PE_BY global p-values side by side, as the article's
#'   extended data-analysis tables report them.
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @seealso [pe_mediation()] (its `report_all_methods` argument),
#'   [pe_mediators()].
#'
#' @family mediation tests
#'
#' @examples
#' set.seed(113)
#' d <- simulate_mediation_data(n = 200, p = 60, outcome = "continuous",
#'                              pattern = "contrasting", c1 = 1)
#' fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous",
#'                     report_all_methods = TRUE)
#' pe_selection(fit)
#'
#' @export
pe_selection <- function(fit) {
  by_method <- attr(fit, "selection_by_method")
  # When only one method was run, build a single-row summary from the stored
  # active set and method label.
  if (is.null(by_method)) {
    active <- attr(fit, "active_mediators")
    by_method <- stats::setNames(list(active), attr(fit, "method"))
  }
  meths <- names(by_method)
  # The per-method PE statistics travel as the pe_by_method attribute; align
  # them to the methods in by_method (NA if, for an old fit, they are absent).
  pbm <- attr(fit, "pe_by_method")
  pick <- function(col) if (is.null(pbm)) rep(NA_real_, length(meths)) else
    pbm[[col]][match(meths, pbm$method)]
  # Label the active set by column name when `M` had names.
  nm <- attr(fit, "mediator_names")
  lab <- function(a) if (!length(a)) "none" else
    paste(if (is.null(nm)) a else nm[a], collapse = ", ")
  data.frame(
    method = meths,
    n_active = vapply(by_method, length, integer(1)),
    j_pe = pick("j_pe"),
    stat_pe = pick("stat_pe"),
    pval_pe = pick("pval_pe"),
    active_mediators = vapply(by_method, lab, character(1)),
    row.names = NULL, stringsAsFactors = FALSE) |> .as_poem_tbl()
}
