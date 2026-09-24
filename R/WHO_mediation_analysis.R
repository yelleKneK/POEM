# Fit one subgroup with the article's exact preprocessing: standardize the
# exposure and mediators, center the outcome, and leave the confounder
# indicators UNSCALED (so the benchmark Wald test follows the article, which
# does not scale the region/income/year covariates). Returns one summary row.
# Not exported.
#' @keywords internal
#' @noRd
.WHO_fit_group <- function(grouping, group, des, lambda_grid, error_level,
                           full_table = FALSE) {
  # Some indicators are constant within a small region/income cell; drop them
  # (they carry no signal and break standardization) and keep the surviving
  # indicator codes so the reported active mediators map back correctly.
  keep <- which(apply(des$M, 2L, stats::sd) > 0)
  M <- des$M[, keep, drop = FALSE]
  ind <- des$indicators[keep]
  fit <- pe_mediation_linear(
    scale(des$X), des$Y - mean(des$Y), scale(M), Z = des$Z,
    scale = FALSE, error_level = error_level,
    report_all_methods = full_table,
    lambda_grid = lambda_grid, lambda_grid_reduced = lambda_grid)
  pval_hdmm <- fit$value[fit$term == "pval_hdmm"]
  codes <- function(idx) if (length(idx))
    paste(ind[idx], collapse = ", ") else "none"

  if (!full_table) {
    active_idx <- attr(fit, "active_mediators")
    return(data.frame(
      grouping = grouping, group = group, n_countries = des$n_countries,
      pval_hdmm = pval_hdmm,
      pval_pe   = fit$value[fit$term == "pval_pe"],
      n_active  = length(active_idx),
      active_mediators = codes(active_idx),
      stringsAsFactors = FALSE))
  }

  # full_table: report the PE p-value and the active mediator codes under each
  # multiplicity method (Bonferroni / BH / BY), as the article's extended
  # data-analysis tables do.
  pbm <- attr(fit, "pe_by_method")
  sbm <- attr(fit, "selection_by_method")
  meth_key <- c(Bonferroni = "bonferroni", BH = "bh", BY = "by")
  row <- data.frame(grouping = grouping, group = group,
                    n_countries = des$n_countries,
                    pval_hdmm = pval_hdmm, stringsAsFactors = FALSE)
  for (m in names(meth_key)) {
    k <- meth_key[[m]]
    row[[paste0("pval_pe_", k)]] <-
      if (!is.null(pbm)) pbm$pval_pe[match(m, pbm$method)] else NA_real_
    idx <- if (!is.null(sbm)) sbm[[m]] else integer(0)
    row[[paste0("active_", k)]] <- codes(idx)
  }
  row
}

#' Reproduce the article's empirical mediation analysis
#'
#' Runs the article's real-data analysis for one health outcome across the
#' groupings it reports: a global model using all WHO members, then
#' separate models within each WHO region, within each World Bank income
#' group, and (optionally) within each region-by-income cell. For every
#' grouping it reports the benchmark Wald test and the power-enhanced test
#' on the total indirect effect of health expenditure, together with the
#' individual mediators the PE test identifies as active. This reproduces
#' the structure of the article's data-analysis
#' tables (for example its Table for infant mortality) from the shipped
#' [WHO_health_mediation] data.
#'
#' The output is the package's \strong{benchmark fit}: applied to the
#' benchmark [WHO_health_mediation] data with the default settings, it
#' returns the active mediators reported in the article, and the
#' package's test suite checks those values, so the fit
#' serves as a fixed reference point for the method.
#'
#' @details
#' Each subgroup is fit with the article's preprocessing: the exposure and
#' the 57 mediators are standardized, the outcome is centered, and the
#' confounder indicators (region, income, year) are left unscaled. Groups
#' with too few complete observations (some region-by-income cells are
#' empty or nearly so) are skipped and reported as `NA` rows, named in a
#' message; a group whose fit fails is likewise reported as `NA`, with a
#' warning that quotes the error. Both are recorded in the `"notes"`
#' attribute (a `data.frame` with `group`, `stage`, and `message`, or
#' `NULL` when every group was fit). Because every cell fits a penalized
#' model over a grid of tuning parameters, a full run with all groupings
#' fits dozens of models and takes about ten seconds on a current laptop; start with the
#' default groupings, or a single outcome, before scaling up.
#'
#' Exact p-values depend on modeling choices that the article fix
#' (the tuning-parameter grid, the unscaled confounders); the active
#' mediators the PE test identifies are the stable, substantive output and
#' reproduce the article's findings (for example, general government
#' expenditure as a percent of GDP, `gge_gdp`, for the global
#' infant-mortality model).
#'
#' @param outcome Which health outcome to analyze: `"imr"`, `"u5mr"`,
#'   `"leb"`, `"lbw"`, or `"pou"`. See [WHO_health_mediation].
#' @param groupings Character vector of groupings to run, any of
#'   `"global"`, `"region"`, `"income"`, and `"region_income"`. Default
#'   `c("global", "region", "income")` (the region-by-income cells are the
#'   slowest and are opt-in).
#' @param error_level Target familywise error rate for the
#'   mediator-screening step. Default 0.05, as in Yu and Kelley (in press).
#' @param lambda_grid Tuning-parameter grid for the subgroup penalized fits
#'   (regional, income, and region-by-income models). Default
#'   `seq(0.1, 10, length.out = 100)`, as in Yu and Kelley (in press).
#' @param lambda_grid_global Tuning-parameter grid for the global (all-country)
#'   model only. Default `seq(0.1, 5, length.out = 100)`, the narrower
#'   range Yu and Kelley (in press) use for the global fit.
#' @param data The source data. Defaults to [WHO_health_mediation].
#' @param cores Number of CPU cores; values above 1 fit the subgroups in
#'   parallel via the base \pkg{parallel} package (Unix only). Default 1.
#' @param verbose Logical; if `TRUE`, print each group as it is fit.
#'   Default `FALSE`.
#' @param full_table Logical; if `TRUE`, return the article's
#'   \emph{extended} table layout, with the power-enhanced p-value and the
#'   active mediator codes reported separately under each multiplicity
#'   method (Bonferroni for FWER, BH and BY for FDR) instead of a single
#'   primary method. Default `FALSE`.
#'
#' @return A `data.frame` with one row per fitted subgroup. With
#'   `full_table = FALSE` (the default) the columns are `grouping`
#'   (global / region / income / region_income), `group` (the specific
#'   group, for example `AFR` or `Low`), `n_countries` (number of countries),
#'   `pval_hdmm` and `pval_pe` (the benchmark and power-enhanced p-values
#'   for the total indirect effect), `n_active` (number of active mediators
#'   identified), and `active_mediators` (their codes, or `"none"`). With
#'   `full_table = TRUE` the single `pval_pe`/`n_active`/`active_mediators`
#'   columns are replaced by `pval_pe_bonferroni`, `pval_pe_bh`,
#'   `pval_pe_by` and the matching `active_bonferroni`, `active_bh`,
#'   `active_by` code columns. The result is a `poem_tbl` with the
#'   attributes `"outcome"`, `"full_table"`, and `"notes"` (see Details);
#'   call [summary()][summary.poem_tbl] on it for a cross-group digest of
#'   where the power-enhanced test detects mediation the benchmark misses.
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @references
#' Yu, X., and Kelley, K. (in press). Power Enhancement in
#' High-Dimensional Heterogeneous Mediation Analysis. \emph{Journal of the
#' American Statistical Association}.
#'
#' @seealso [WHO_mediation_design()] for a single design,
#'   [WHO_health_mediation] for the data, [pe_mediation()] for the test.
#'
#' @family WHO mediation data
#'
#' @examples
#' # The global infant-mortality model, the article's headline row (under a
#' # second). Add "region" and "income" to groupings for the full table.
#' WHO_mediation_analysis("imr", groupings = "global")
#'
#' @export
WHO_mediation_analysis <- function(outcome = c("imr", "u5mr", "leb",
                                               "lbw", "pou"),
                                   groupings = c("global", "region", "income"),
                                   error_level = 0.05,
                                   lambda_grid = seq(0.1, 10, length.out = 100),
                                   lambda_grid_global = seq(0.1, 5,
                                                            length.out = 100),
                                   data = WHO_health_mediation,
                                   cores = 1L, verbose = FALSE,
                                   full_table = FALSE) {
  outcome <- match.arg(outcome)
  groupings <- match.arg(groupings,
                         c("global", "region", "income", "region_income"),
                         several.ok = TRUE)
  .validate_error_level(error_level)
  .validate_lambda_grid(lambda_grid, "lambda_grid")
  .validate_lambda_grid(lambda_grid_global, "lambda_grid_global")
  if (!is.numeric(cores) || length(cores) != 1L || !is.finite(cores) ||
      cores < 1 || cores != round(cores))
    stop("`cores` must be a single positive whole number.", call. = FALSE)
  for (nm in c("verbose", "full_table")) {
    v <- get(nm)
    if (!is.logical(v) || length(v) != 1L || is.na(v))
      stop(sprintf("`%s` must be TRUE or FALSE.", nm), call. = FALSE)
  }
  .WHO_validate_data(data, outcome)

  regions <- c("AFR", "AMR", "EMR", "EUR", "SEAR", "WPR")
  incomes <- c("Low", "Lower-middle", "Upper-middle", "High")

  # Build the (grouping, group, region-filter, income-filter) work list from
  # the requested groupings.
  jobs <- list()
  if ("global" %in% groupings)
    jobs <- c(jobs, list(list(g = "global", lab = "ALL", reg = NULL, inc = NULL)))
  if ("region" %in% groupings)
    jobs <- c(jobs, lapply(regions,
              function(r) list(g = "region", lab = r, reg = r, inc = NULL)))
  if ("income" %in% groupings)
    jobs <- c(jobs, lapply(incomes,
              function(i) list(g = "income", lab = i, reg = NULL, inc = i)))
  if ("region_income" %in% groupings)
    jobs <- c(jobs, unlist(lapply(regions, function(r) lapply(incomes,
              function(i) list(g = "region_income", lab = paste(r, i, sep = "-"),
                               reg = r, inc = i))), recursive = FALSE))

  # Fit one grouping cell, returning a one-row summary. Cells with too few
  # observations (WHO_mediation_design() errors) or a failed fit are recorded as
  # NA rather than aborting the whole run; the reason travels as a "note"
  # attribute on the row and is reported once, after the loop.
  fit_job <- function(job) {
    if (verbose) message("Fitting ", outcome, " / ", job$lab, " ...")
    na_row <- function(n_countries) {
      if (!full_table)
        return(data.frame(
          grouping = job$g, group = job$lab, n_countries = n_countries,
          pval_hdmm = NA_real_, pval_pe = NA_real_, n_active = NA_integer_,
          active_mediators = NA_character_, stringsAsFactors = FALSE))
      data.frame(
        grouping = job$g, group = job$lab, n_countries = n_countries,
        pval_hdmm = NA_real_,
        pval_pe_bonferroni = NA_real_, active_bonferroni = NA_character_,
        pval_pe_bh = NA_real_, active_bh = NA_character_,
        pval_pe_by = NA_real_, active_by = NA_character_,
        stringsAsFactors = FALSE)
    }
    noted <- function(row, stage, cond) {
      attr(row, "note") <- c(stage = stage, message = conditionMessage(cond))
      row
    }
    des <- tryCatch(
      WHO_mediation_design(outcome, region = job$reg, income = job$inc,
                           data = data),
      error = function(e) e)
    if (inherits(des, "error"))
      return(noted(na_row(NA_integer_), "design", des))
    grid <- if (job$g == "global") lambda_grid_global else lambda_grid
    res <- tryCatch(
      .WHO_fit_group(job$g, job$lab, des, grid, error_level,
                     full_table = full_table),
      error = function(e) e)
    if (inherits(res, "error"))
      return(noted(na_row(des$n_countries), "fit", res))
    res
  }
  rows <- .pe_lapply(jobs, fit_job, cores = cores)

  # Collect the per-cell notes (skipped designs, failed fits) and report
  # them once: skips as a message (they are expected in sparse cells), fit
  # failures as a warning quoting the first error.
  notes <- do.call(rbind, lapply(seq_along(rows), function(i) {
    nt <- attr(rows[[i]], "note")
    if (is.null(nt)) return(NULL)
    data.frame(group = jobs[[i]]$lab, stage = nt[["stage"]],
               message = nt[["message"]], stringsAsFactors = FALSE)
  }))
  .WHO_report_notes(notes)

  out <- do.call(rbind, lapply(rows, function(r) { attr(r, "note") <- NULL; r }))
  rownames(out) <- NULL
  # Return the package's general poem_tbl so summary() gives the cross-group
  # benchmark-vs-PE digest; it still behaves as an ordinary data.frame.
  out <- .as_poem_tbl(out)
  attr(out, "outcome") <- outcome
  attr(out, "full_table") <- full_table
  attr(out, "notes") <- notes
  out
}

# The data must be a data.frame in the shape of WHO_health_mediation: the key
# and confounder columns, the exposure, the requested outcome, and at least
# two indicator columns. Checked once up front so a wrong `data` stops with a
# clear message instead of every cell silently coming back NA. Not exported.
#' @keywords internal
#' @noRd
.WHO_validate_data <- function(data, outcome) {
  if (!is.data.frame(data))
    stop("`data` must be a data.frame in the shape of WHO_health_mediation.",
         call. = FALSE)
  need <- c("code", "region", "income", "year", "gdp_growth", outcome)
  miss <- setdiff(need, names(data))
  if (length(miss))
    stop("`data` lacks the column", if (length(miss) > 1L) "s" else "", " ",
         paste0("`", miss, "`", collapse = ", "),
         " that WHO_health_mediation carries.", call. = FALSE)
  if (length(.WHO_indicator_names(data)) < 2L)
    stop("`data` carries fewer than two indicator (mediator) columns.",
         call. = FALSE)
  invisible(TRUE)
}

# Report the per-cell notes collected by WHO_mediation_analysis(). Not
# exported.
#' @keywords internal
#' @noRd
.WHO_report_notes <- function(notes) {
  if (is.null(notes) || !nrow(notes)) return(invisible(NULL))
  label <- function(df) paste(sprintf("%s (%s)", df$group,
                                      sub("[.]$", "", df$message)),
                              collapse = "; ")
  skipped <- notes[notes$stage == "design", , drop = FALSE]
  failed  <- notes[notes$stage == "fit", , drop = FALSE]
  if (nrow(skipped))
    message(sprintf("%d group%s skipped and reported as NA: %s",
                    nrow(skipped), if (nrow(skipped) == 1L) "" else "s",
                    label(skipped)))
  if (nrow(failed))
    warning(sprintf("The fit failed in %d group%s, reported as NA: %s",
                    nrow(failed), if (nrow(failed) == 1L) "" else "s",
                    label(failed)), call. = FALSE)
  invisible(NULL)
}
