#' Autoregressive (AR(1)) covariance matrix for mediator simulation
#'
#' Builds the \eqn{p \times p} first-order autoregressive covariance
#' matrix \eqn{\Sigma = (\rho^{|i-j|})_{i,j}}, the intercorrelation
#' structure imposed on the mediator noise throughout the article's
#' Monte Carlo studies. It is exported as a small reusable utility so a
#' user can inspect or reuse the exact covariance the simulations use.
#'
#' @param p Number of mediators (matrix dimension); a positive integer.
#' @param rho Autocorrelation parameter in \eqn{(-1, 1)}. The covariance
#'   between mediators \eqn{i} and \eqn{j} is \eqn{\rho^{|i-j|}}, so
#'   nearby mediators are more strongly correlated and the correlation
#'   decays geometrically with separation. Default 0.5, the value used in
#'   the article.
#'
#' @return A \eqn{p \times p} numeric matrix with ones on the diagonal
#'   and \eqn{\rho^{|i-j|}} off the diagonal.
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @seealso [simulate_mediation_data()], which uses this covariance to
#'   draw correlated mediator noise.
#'
#' @family mediation simulation
#'
#' @examples
#' mediation_ar1_cov(4, rho = 0.5)
#'
#' @export
mediation_ar1_cov <- function(p, rho = 0.5) {
  if (!is.numeric(p) || length(p) != 1L || p < 1 || p != round(p))
    stop("`p` must be a single positive integer.", call. = FALSE)
  if (!is.numeric(rho) || length(rho) != 1L || abs(rho) >= 1)
    stop("`rho` must be a single number in (-1, 1).", call. = FALSE)
  # |i - j| as a matrix of lag distances, then rho raised to it. The outer()
  # form builds the whole banded matrix at once.
  lag <- abs(outer(seq_len(p), seq_len(p), "-"))
  rho^lag
}
