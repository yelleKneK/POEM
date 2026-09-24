#' Example high-dimensional mediation data sets
#'
#' Three small simulated data sets, one per outcome type, for trying the
#' POEM tests and reproducing the help-page examples. Each was generated
#' by [simulate_mediation_data()] under a homogeneous mediation pattern
#' (a few same-sign active mediators among many null ones) with three
#' confounders, at seed 113, so the generating truth is known and travels
#' with the data.
#'
#' @format Each is a list, the return value of
#'   [simulate_mediation_data()], with components:
#' \describe{
#'   \item{X}{Exposure matrix, \eqn{n \times 1}.}
#'   \item{M}{Candidate mediator matrix, \eqn{n \times 100}.}
#'   \item{Y}{Outcome vector of length \eqn{n}: continuous for
#'     `example_continuous`, 0/1 for `example_binary`, nonnegative counts
#'     for `example_count`.}
#'   \item{Z}{Confounder matrix, \eqn{n \times 3}.}
#'   \item{alpha_m, Gamma_x, alpha_x}{The mediator-on-outcome,
#'     exposure-on-mediator, and direct-effect coefficients used to
#'     generate the data.}
#'   \item{beta}{The total indirect effect \eqn{\Gamma_x \alpha_m}.}
#'   \item{active_mediators}{Indices of the truly active mediators.}
#'   \item{n, p, q, outcome, pattern}{The generating settings.}
#' }
#'   The sample sizes are \eqn{n = 200} (`example_continuous`) and
#'   \eqn{n = 300} (`example_binary`, `example_count`).
#'
#' @source Simulated with [simulate_mediation_data()]; see
#'   `data-raw/example_data.R` in the package sources.
#'
#' @seealso [pe_mediation()], [simulate_mediation_data()].
#'
#' @examples
#' data(example_continuous)
#' str(example_continuous$active_mediators)
#' fit <- pe_mediation(example_continuous$X, example_continuous$Y,
#'                     example_continuous$M, Z = example_continuous$Z,
#'                     outcome = "continuous")
#' fit
#'
#' @name example_data
#' @aliases example_continuous example_binary example_count
NULL
