# The non-key, non-outcome columns of WHO_health_mediation are exactly the 57
# candidate mediators. Computed once here so the design builder and the data
# documentation agree on the set.
#' @keywords internal
#' @noRd
.WHO_indicator_names <- function(data) {
  setdiff(names(data),
          c("code", "country", "region", "income", "year",
            "gdp_growth", "imr", "u5mr", "leb", "lbw", "pou"))
}

#' Build a mediation design from the WHO health-expenditure data
#'
#' Assembles the exposure, outcome, mediator, and confounder matrices for
#' one of the analyses in the article from [WHO_health_mediation], ready
#' to pass to [pe_mediation()]. The exposure is GDP per capita growth, the
#' mediators are the 57 health-expenditure indicators, and the outcome is
#' the chosen health indicator. Optionally restrict to a WHO region and/or
#' a World Bank income group, reproducing the article's region-specific,
#' income-specific, and region-by-income models.
#'
#' @details
#' Rows with a missing outcome or missing exposure are dropped (the
#' article handles each outcome's coverage separately, which is why the
#' low-birthweight and undernourishment outcomes have fewer observations).
#' The confounder matrix \eqn{Z} contains the WHO region and income group
#' (as indicator variables) and the year (numeric), exactly as in the
#' article, but only for covariates that still vary after subsetting: a
#' region-specific model drops the (now constant) region, and unused
#' factor levels are dropped so no all-zero indicator columns remain.
#'
#' To reproduce the article's exact test statistics, fit with the
#' article's preprocessing, standardizing the exposure and mediators and
#' centering the outcome while leaving the confounder indicators
#' unscaled:
#' \preformatted{
#'   des <- WHO_mediation_design("imr")
#'   pe_mediation(scale(des$X), des$Y - mean(des$Y), scale(des$M),
#'                Z = des$Z, outcome = "continuous", scale = FALSE,
#'                lambda_grid = seq(0.1, 5, length.out = 100))
#' }
#' The simpler call `pe_mediation(des$X, des$Y, des$M, des$Z, outcome =
#' "continuous")` standardizes every block, including \eqn{Z}, and so
#' fits a different model: scaling the indicator columns changes both the
#' benchmark Wald test and the penalized selection, so its p-values and
#' its active set need not match the article's (on the default tuning
#' grid the global IMR active set happens to agree; on the article's
#' 100-point grid it does not). [WHO_mediation_analysis()] applies the
#' article's preprocessing for you.
#'
#' @param outcome Which health outcome to use as \eqn{Y}: `"imr"` (infant
#'   mortality rate), `"u5mr"` (under-five mortality rate), `"leb"` (life
#'   expectancy at birth), `"lbw"` (prevalence of low birthweight), or
#'   `"pou"` (prevalence of undernourishment).
#' @param region Optional character vector of WHO regions to keep
#'   (`"AFR"`, `"AMR"`, `"EMR"`, `"EUR"`, `"SEAR"`, `"WPR"`). Default
#'   `NULL` keeps all regions.
#' @param income Optional character vector of income groups to keep
#'   (`"Low"`, `"Lower-middle"`, `"Upper-middle"`, `"High"`). Default
#'   `NULL` keeps all groups.
#' @param data The source data frame. Defaults to [WHO_health_mediation].
#'
#' @return A list with `X` (exposure, \eqn{n \times 1}), `Y` (outcome
#'   vector), `M` (the \eqn{n \times 57} mediator matrix), `Z` (the
#'   confounder matrix, or `NULL` if no covariate varies), and the
#'   metadata `outcome`, `n` (country-year observations), `n_countries`
#'   (unique countries), `indicators`, `region`, and `income`.
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @seealso [pe_mediation()], [WHO_health_mediation],
#'   [WHO_indicator_codebook].
#'
#' @family WHO mediation data
#'
#' @examples
#' # Global infant-mortality model: GDP growth -> health spending -> IMR.
#' des <- WHO_mediation_design("imr")
#' c(n = des$n, mediators = ncol(des$M))
#' fit <- pe_mediation(scale(des$X), des$Y - mean(des$Y), scale(des$M),
#'   Z = des$Z, outcome = "continuous", scale = FALSE,
#'   lambda_grid = seq(0.1, 5, length.out = 100))
#' fit
#' # The active mediator's full name:
#' WHO_indicator_codebook$description[
#'   match(des$indicators[attr(fit, "active_mediators")],
#'         WHO_indicator_codebook$indicator)]
#'
#' # A region-specific model (Europe).
#' des_eur <- WHO_mediation_design("imr", region = "EUR")
#'
#' @export
WHO_mediation_design <- function(outcome = c("imr", "u5mr", "leb",
                                             "lbw", "pou"),
                                 region = NULL, income = NULL,
                                 data = WHO_health_mediation) {
  outcome <- match.arg(outcome)
  df <- data
  if (!is.null(region)) {
    bad <- setdiff(region, levels(factor(df$region)))
    if (length(bad)) stop("Unknown region(s): ", paste(bad, collapse = ", "),
                          ".", call. = FALSE)
    df <- df[df$region %in% region, , drop = FALSE]
  }
  if (!is.null(income)) {
    bad <- setdiff(income, levels(df$income))
    if (length(bad)) stop("Unknown income group(s): ", paste(bad, collapse = ", "),
                          ".", call. = FALSE)
    df <- df[as.character(df$income) %in% income, , drop = FALSE]
  }
  # Drop observations with no outcome or no exposure value.
  df <- df[!is.na(df[[outcome]]) & !is.na(df$gdp_growth), , drop = FALSE]
  if (nrow(df) < 10L)
    stop("Too few complete observations (", nrow(df),
         ") for this outcome/region/income combination.", call. = FALSE)

  indicators <- .WHO_indicator_names(data)
  X <- matrix(df$gdp_growth, ncol = 1L, dimnames = list(NULL, "gdp_growth"))
  Y <- df[[outcome]]
  M <- as.matrix(df[, indicators])

  # Confounders: region and income indicators (only if they still vary after
  # subsetting; unused factor levels dropped so no zero columns appear) plus
  # year as a numeric trend, the article's covariate set.
  terms <- character(0)
  reg <- droplevels(factor(as.character(df$region)))
  inc <- droplevels(factor(as.character(df$income)))
  if (nlevels(reg) > 1L) { df$.region <- reg; terms <- c(terms, ".region") }
  if (nlevels(inc) > 1L) { df$.income <- inc; terms <- c(terms, ".income") }
  terms <- c(terms, "year")
  Z <- stats::model.matrix(stats::reformulate(terms), data = df)[, -1,
                                                                 drop = FALSE]
  if (ncol(Z) == 0L) Z <- NULL

  list(X = X, Y = Y, M = M, Z = Z, outcome = outcome, n = nrow(df),
       n_countries = length(unique(df$code)),
       indicators = indicators, region = region, income = income)
}
