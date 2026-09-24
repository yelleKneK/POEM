#' Power Enhancement for High-Dimensional Mediation Analysis (POEM)
#'
#' Powerful global tests for high-dimensional mediation analysis, for
#' continuous, binary, and count outcomes. The package answers the
#' question \dQuote{is there any active mediator among a large set of
#' candidate mediators?} and, when the answer is yes, reports which
#' individual mediators are active. The name comes from POwer-Enhanced
#' Mediation.
#'
#' @section The problem POEM solves:
#' A mediation analysis asks how an exposure \eqn{X} affects an outcome
#' \eqn{Y} through intermediate variables (mediators) \eqn{M}. With many
#' candidate mediators, the standard approach tests the \emph{total}
#' indirect effect \eqn{\beta = \Gamma_x \alpha_m}, the sum of every
#' mediator's individual indirect effect. That test has a blind spot:
#' when some mediators carry positive indirect effects and others carry
#' negative ones, the individual effects can cancel and the total
#' indirect effect is zero even though active mediators plainly exist.
#' A test built on the total indirect effect is then powerless, by
#' construction, against exactly the heterogeneous mediation that is
#' common in practice.
#'
#' The power enhancement (PE) tests in this package add a component
#' \eqn{J_m} that accumulates the marginal signal from each individual
#' mediator. Because \eqn{J_m} is a sum of magnitudes, opposite-signed
#' indirect effects reinforce rather than cancel, so the PE test detects
#' active mediators whether their effects are homogeneous,
#' heterogeneous, or exactly contrasting, while keeping the Type I error
#' rate of the original test under the global null of no mediation.
#'
#' @section Main functions:
#' \describe{
#'   \item{[pe_mediation()]}{The front end. Give it the exposure,
#'     outcome, mediators, and optional confounders, name the outcome
#'     type (continuous, binary, or count), and it returns a tidy table
#'     with the benchmark Wald test, the power-enhanced test, the
#'     estimated total indirect effect, and the set of active mediators
#'     it identified.}
#'   \item{[pe_mediation_linear()], [pe_mediation_logistic()],
#'     [pe_mediation_poisson()]}{The outcome-specific workers that
#'     [pe_mediation()] dispatches to. Call them directly when the
#'     outcome type is known.}
#'   \item{[simulate_mediation_data()]}{Generates mediation data under
#'     the homogeneous, heterogeneous, and contrasting patterns studied
#'     in the article, for any of the three outcome types.}
#'   \item{[pe_power_curve()]}{Runs the article's Monte Carlo size and
#'     power study: it sweeps a signal-strength grid and returns the
#'     empirical rejection rate of the benchmark and power-enhanced
#'     tests at each point.}
#'   \item{[pe_mediate()]}{A data-frame / formula front end to
#'     [pe_mediation()] for users who prefer naming columns to matrices.}
#'   \item{[pe_mediators()] and [pe_selection()]}{The per-mediator screen
#'     detail behind a fit, and the active set under each multiplicity
#'     method.}
#'   \item{[pe_simulation_study()] and [ss_power_pe_mediation()]}{Run the
#'     size and power study across patterns, and plan the sample size for
#'     a target power.}
#'   \item{[plot.poem_tbl()], [tidy.poem_tbl()], [glance.poem_tbl()]}{A
#'     base-graphics plot and broom verbs for POEM results.}
#'   \item{[mediation_ar1_cov()]}{The autoregressive mediator covariance
#'     used throughout the simulations, exported as a reusable utility.}
#'   \item{[WHO_health_mediation], [WHO_mediation_design()], and
#'     [WHO_mediation_analysis()]}{The article's empirical data (economic
#'     growth, 57 health-expenditure mediators, five health outcomes for
#'     91 WHO members), a helper that assembles an analysis design from
#'     it, and a function that reproduces the article's full set of
#'     subgroup analyses.}
#' }
#'
#' @section Reading the output:
#' Every estimation and testing function returns a tidy `data.frame`
#' with a `term` column and a numeric `value` column, carrying the
#' `poem_tbl` class so it prints with whole numbers shown without a
#' decimal part, other quantities to a few significant figures, and
#' p-values to four decimal places. The stored numbers keep full
#' precision; only the display rounds. See [poem_tbl] for the details.
#'
#' @references
#' Yu, X., and Kelley, K. (in press). Power Enhancement in
#' High-Dimensional Heterogeneous Mediation Analysis. \emph{Journal of the
#' American Statistical Association}.
#'
#' Fan, J., Liao, Y., and Yao, J. (2015). Power enhancement in
#' high-dimensional cross-sectional tests. \emph{Econometrica, 83}(4),
#' 1497--1541. \doi{10.3982/ECTA12749}
#'
#' Guo, X., Li, R., Liu, J., and Zeng, M. (2022). High-dimensional
#' mediation analysis for selecting DNA methylation loci mediating
#' childhood trauma and cortisol stress reactivity. \emph{Journal of
#' the American Statistical Association, 117}(539), 1110--1121.
#' \doi{10.1080/01621459.2022.2053136}
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @keywords internal
"_PACKAGE"
