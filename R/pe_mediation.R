#' Power-enhanced test for high-dimensional mediation
#'
#' Tests whether any mediator is active among a large, possibly
#' intercorrelated set of candidate mediators, for a continuous, binary,
#' or count outcome. This is the main entry point of \pkg{POEM}: name the
#' outcome type and it dispatches to the appropriate model, returning one
#' tidy table that reports the benchmark Wald test on the total indirect
#' effect side by side with the power-enhanced (PE) test, plus the set of
#' individual mediators identified as active.
#'
#' @details
#' Standard high-dimensional mediation tests target the \emph{total}
#' indirect effect \eqn{\beta = \Gamma_x \alpha_m}, the sum of every
#' mediator's individual indirect effect. When some indirect effects are
#' positive and others negative they can cancel, making \eqn{\beta = 0}
#' even though active mediators exist, and a test built on \eqn{\beta} is
#' then powerless. The power-enhanced test adds a component
#' \deqn{J_m = \sqrt{p} \sum_{i=1}^{q} \sum_{j \in \hat{S}}
#'   \left| \frac{\hat\alpha_{m,j}}{\hat\sigma_{m,j}} \right|
#'   \left| \frac{\hat\Gamma_{x,i,j}}{\hat\sigma_{\Gamma,i,j}} \right|
#'   \mathbf{1}\!\left\{ \max(p_{1,j}, p_{2,i,j}) <
#'   \frac{\alpha_{\mathrm{lvl}}}{s q \log\log n} \right\},}
#' a sum over the selected mediators \eqn{\hat S} of the product of the
#' (standardized) mediator-on-outcome and exposure-on-mediator
#' statistics, kept only for pairs where both paths are individually
#' significant. Because \eqn{J_m} accumulates magnitudes, opposite-signed
#' indirect effects reinforce rather than cancel, so the test
#' \eqn{M_{PE} = S_n + J_m} stays powerful under heterogeneous and
#' contrasting mediation while keeping the chi-square reference
#' distribution (and hence the Type I error rate) of the benchmark test
#' \eqn{S_n} under the global null. The indicator also identifies which
#' individual mediators are active, with familywise error rate control
#' under `method = "Bonferroni"` or false discovery rate control under
#' `"BH"` / `"BY"`.
#'
#' @param X Numeric matrix of exposures with \eqn{n} rows and \eqn{q}
#'   columns (\eqn{q} is usually 1). A numeric vector is treated as a
#'   single exposure.
#' @param Y Numeric outcome vector of length \eqn{n}. Continuous for
#'   `outcome = "continuous"`, coded 0/1 for `"binary"`, and nonnegative
#'   integer counts for `"count"`.
#' @param M Numeric matrix of candidate mediators with \eqn{n} rows and
#'   \eqn{p} columns; \eqn{p} may exceed \eqn{n}. Column names, when
#'   present, label the active mediators in the print footer and the
#'   mediator table; without them mediators are reported by column
#'   position.
#' @param Z Optional numeric matrix of confounders with \eqn{n} rows.
#'   Default `NULL` (no confounders).
#' @param outcome Outcome type: `"continuous"` (linear model),
#'   `"binary"` (logistic model), or `"count"` (Poisson log-link model).
#' @param method Multiplicity adjustment used to screen individual
#'   mediators in the PE component: `"Bonferroni"` (familywise error
#'   rate control, the default), `"BH"`, or `"BY"` (false discovery rate
#'   control; `"BH"` assumes independence or positive dependence among
#'   mediators, `"BY"` is valid under arbitrary dependence).
#' @param scale Logical; if `TRUE` (default) the exposures, mediators,
#'   and confounders are standardized (and a continuous outcome is
#'   centered) before fitting, as the method assumes.
#' @param error_level Target error rate for the mediator-screening step,
#'   in \eqn{(0, 1)}. Interpreted as the familywise error rate when
#'   `method = "Bonferroni"` and as the false discovery rate when
#'   `method = "BH"` or `"BY"`. Default 0.05. This is distinct from the
#'   significance level used to test the global null (which is the user's
#'   choice when reading `pval_pe`); the test is insensitive to
#'   `error_level` and remains valid for any fixed value in \eqn{(0, 1)}.
#' @param conf_level Confidence level for the interval on the total
#'   indirect effect. Default 0.95.
#' @param report_all_methods Logical; if `TRUE`, the active-mediator set is
#'   computed for all three multiplicity methods (Bonferroni, BH, BY) off the
#'   same fit and recorded for comparison, retrievable with [pe_selection()].
#'   The `method` argument still drives the primary result. Default `FALSE`.
#' @param drop_constant Logical; how to handle constant (zero-variance)
#'   mediator columns, which cannot be standardized or carry signal. If
#'   `FALSE` (the default) they are an error; if `TRUE` they are dropped
#'   with a warning and the remaining mediators are renumbered back to their
#'   original column positions in the reported active set. Useful for small
#'   subgroups where some mediators happen to be constant.
#' @param lambda_grid Numeric vector of candidate SCAD tuning parameters
#'   for the penalized mediator fit; the fit is repeated at each value and
#'   the one minimizing the high-dimensional BIC (HBIC) is kept. Default
#'   `NULL` uses, for a continuous outcome, [pe_lambda_grid()]: the grid
#'   the article's simulations used at their design, rescaled to this `n`
#'   and `p` by the rate \eqn{\sqrt{\log p / n}} the theory requires of the
#'   tuning parameter; for a binary or count outcome it uses the review-era
#'   grid `seq(0.05, 1, length.out = 20)`. The grid's lower end matters
#'   most, because the HBIC minimum often sits there and because it decides
#'   whether the identified set keeps its error-rate guarantee; see
#'   [pe_lambda_grid()] for the trade-off with power, measured. The value
#'   chosen is reported in the `"tuning"` attribute and the print footer.
#' @param lambda_grid_reduced Numeric vector of candidate tuning
#'   parameters for the reduced-model fit that the continuous-outcome Wald
#'   test refits (unused for binary and count outcomes). Default `NULL`
#'   uses [pe_lambda_grid()] with `model = "reduced"` when `lambda_grid` is
#'   also `NULL` (the article's reduced-model grid, rescaled), and
#'   otherwise `lambda_grid` itself.
#'
#' @return A tidy `data.frame` of class `poem_tbl`. See
#'   [pe_mediation_linear()] for the row schema and the attributes. The
#'   identified active mediators are in the `"active_mediators"` attribute
#'   and the print footer; the tuning parameter HBIC chose is in the
#'   `"tuning"` attribute and the footer.
#'
#' @section When to use POEM:
#' POEM is for the global question, \dQuote{is there any active mediator
#' among many candidates?}, followed by the identification of the active
#' ones. It assumes a linear or generalized linear outcome model, a sparse
#' set of active mediators, standardized inputs, and complete data. Some
#' practical limits, seen on simulated data during the package's release
#' audit:
#' \itemize{
#'   \item Sample size. The penalized selection needs \eqn{n} large
#'     relative to \eqn{\log p} and to the signal. With \eqn{n} near 100
#'     and \eqn{p} in the thousands under the contrasting pattern the
#'     selection keeps a single mediator and the power-enhanced test has
#'     nothing to add. The default grid's floor rises as \eqn{n} falls
#'     (see [pe_lambda_grid()]), which protects the identified set but
#'     costs global power at small designs.
#'   \item Strongly correlated mediators. With an autoregressive
#'     correlation of 0.9 among neighboring mediators the selection keeps
#'     one representative of a correlated block (an active mediator was
#'     recovered in 3 of 10 seeds, against 10 of 10 at a correlation of
#'     0.5), so the identified set is unstable from sample to sample even
#'     when the global test rejects.
#'   \item Rare binary outcomes. With about 2 to 3 percent events at
#'     \eqn{n = 300} the penalized logistic fit selects nothing at any grid
#'     value; the fit is empty and the function says so.
#'   \item Overdispersed counts. The count model is Poisson; under negative
#'     binomial overdispersion its size was not inflated at \eqn{n = 200},
#'     \eqn{p = 100}, but about a tenth of the fits were empty. The count
#'     path is also the slowest, 15 to 30 times the continuous one.
#'   \item Missing values are not handled; supply complete cases.
#' }
#' For a single mediator, or a few mediators fit as one structural model
#' (indirect effects with confidence intervals, likelihood ratio tests of
#' arbitrary indirect effects, moderated mediation), the \pkg{DMAR}
#' package is the tool; POEM's contribution is the high-dimensional global
#' test. The other high-dimensional tests the article benchmarks against
#' (HILMA, GlobalTest, HDMT, DACT) live in their own packages; the vignette
#' `poem-vs-competitors` shows how to run them beside a POEM fit.
#'
#' @references
#' Yu, X., and Kelley, K. (in press). Power Enhancement in
#' High-Dimensional Heterogeneous Mediation Analysis. \emph{Journal of the
#' American Statistical Association}. (The article these methods
#' implement.)
#'
#' Fan, J., Liao, Y., and Yao, J. (2015). Power enhancement in
#' high-dimensional cross-sectional tests. \emph{Econometrica, 83}(4),
#' 1497--1541. \doi{10.3982/ECTA12749}
#'
#' Guo, X., Li, R., Liu, J., and Zeng, M. (2022). High-dimensional
#' mediation analysis for selecting DNA methylation loci mediating
#' childhood trauma and cortisol stress reactivity. \emph{Journal of the
#' American Statistical Association, 117}(539), 1110--1121.
#' \doi{10.1080/01621459.2022.2053136}
#'
#' Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery
#' rate: A practical and powerful approach to multiple testing.
#' \emph{Journal of the Royal Statistical Society, Series B, 57}(1),
#' 289--300. \doi{10.1111/j.2517-6161.1995.tb02031.x}
#'
#' Benjamini, Y., and Yekutieli, D. (2001). The control of the false
#' discovery rate in multiple testing under dependency. \emph{The Annals
#' of Statistics, 29}(4), 1165--1188. \doi{10.1214/aos/1013699998}
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @seealso The outcome-specific workers [pe_mediation_linear()],
#'   [pe_mediation_logistic()], [pe_mediation_poisson()];
#'   [pe_lambda_grid()] for the default tuning grid;
#'   [simulate_mediation_data()] to generate data; [pe_power_curve()] to
#'   reproduce the article's size and power studies.
#'
#' @family mediation tests
#'
#' @examples
#' set.seed(113)
#' # Contrasting mediation: active mediators whose effects cancel, so the
#' # total indirect effect is zero. Compare the two p-values: the benchmark
#' # Wald test (pval_hdmm) is large while the PE test (pval_pe) is tiny.
#' d <- simulate_mediation_data(n = 200, p = 60, pattern = "contrasting",
#'                              outcome = "continuous", c1 = 1)
#' pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
#'
#' @export
pe_mediation <- function(X, Y, M, Z = NULL,
                         outcome = c("continuous", "binary", "count"),
                         method = c("Bonferroni", "BH", "BY"),
                         scale = TRUE, error_level = 0.05,
                         conf_level = 0.95, report_all_methods = FALSE,
                         drop_constant = FALSE,
                         lambda_grid = NULL,
                         lambda_grid_reduced = NULL) {
  outcome <- match.arg(outcome)
  method <- match.arg(method)
  worker <- switch(outcome,
                   continuous = pe_mediation_linear,
                   binary     = pe_mediation_logistic,
                   count      = pe_mediation_poisson)
  worker(X = X, Y = Y, M = M, Z = Z, method = method, scale = scale,
         error_level = error_level, conf_level = conf_level,
         report_all_methods = report_all_methods, drop_constant = drop_constant,
         lambda_grid = lambda_grid, lambda_grid_reduced = lambda_grid_reduced)
}
