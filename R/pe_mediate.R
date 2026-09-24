#' Power-enhanced mediation test from a data frame (formula interface)
#'
#' A convenience front end to [pe_mediation()] for users who have a data
#' frame rather than ready-made matrices. Give it the outcome and exposure
#' as a formula, name the mediator columns, and (optionally) name
#' confounder columns; it assembles the exposure, outcome, mediator, and
#' confounder matrices (turning factor confounders into indicator
#' variables) and runs the test. This is the gentler entry point; the
#' matrix interface [pe_mediation()] is the workhorse.
#'
#' @param formula A two-sided formula `outcome ~ exposure` naming the
#'   outcome and one or more exposure columns in `data`.
#' @param data A `data.frame` containing the outcome, exposure, mediator,
#'   and any confounder columns.
#' @param mediators Character vector of mediator column names in `data`
#'   (the candidate mediators `M`).
#' @param confounders Optional character vector of confounder column names.
#'   Factor and character columns are expanded into indicator variables;
#'   numeric columns enter as is. Default `NULL`.
#' @param outcome Outcome type: `"continuous"`, `"binary"`, or `"count"`.
#' @param ... Further arguments passed to [pe_mediation()] (for example
#'   `method`, `error_level`, `conf_level`, `report_all_methods`,
#'   `lambda_grid`).
#'
#' @return The tidy `poem_tbl` returned by [pe_mediation()].
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @seealso [pe_mediation()] for the matrix interface,
#'   [WHO_mediation_design()] for the worked health-expenditure design.
#'
#' @family mediation tests
#'
#' @examples
#' # Build a small data frame and test through the formula interface.
#' set.seed(113)
#' d <- simulate_mediation_data(n = 200, p = 60, outcome = "continuous",
#'                              pattern = "contrasting", c1 = 1)
#' df <- data.frame(y = d$Y, x = d$X[, 1], d$M)
#' meds <- grep("^X", names(df), value = TRUE)   # the mediator columns
#' # The mediator columns are named, so the footer names the active ones.
#' pe_mediate(y ~ x, data = df, mediators = meds, outcome = "continuous")
#'
#' @export
pe_mediate <- function(formula, data,
                       mediators, confounders = NULL,
                       outcome = c("continuous", "binary", "count"), ...) {
  outcome <- match.arg(outcome)
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("`formula` must be a two-sided formula, outcome ~ exposure.",
         call. = FALSE)
  if (!is.data.frame(data))
    stop("`data` must be a data.frame.", call. = FALSE)

  # Outcome and exposures come from the formula (the exposure side may name
  # several columns, e.g. y ~ x1 + x2). model.frame drops incomplete rows.
  yname <- all.vars(formula)[1L]
  mf <- stats::model.frame(formula, data = data)
  Y <- stats::model.response(mf)
  X <- stats::model.matrix(formula, data = mf)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]

  miss <- setdiff(mediators, names(data))
  if (length(miss))
    stop("mediator column(s) not in `data`: ", paste(miss, collapse = ", "),
         ".", call. = FALSE)
  M <- as.matrix(data[rownames(mf), mediators, drop = FALSE])

  # Confounders: numeric columns pass through, factor/character columns become
  # indicator variables (no intercept column).
  Z <- NULL
  if (!is.null(confounders)) {
    miss <- setdiff(confounders, names(data))
    if (length(miss))
      stop("confounder column(s) not in `data`: ", paste(miss, collapse = ", "),
           ".", call. = FALSE)
    zf <- stats::reformulate(confounders)
    Z <- stats::model.matrix(zf, data = data[rownames(mf), , drop = FALSE])
    Z <- Z[, colnames(Z) != "(Intercept)", drop = FALSE]
  }

  pe_mediation(X = X, Y = Y, M = M, Z = Z, outcome = outcome, ...)
}
