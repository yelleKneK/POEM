#' Calibration constants for the real-data-motivated heterogeneous setting
#'
#' The fixed constants for the article's real-data-motivated heterogeneous
#' mediation simulation, calibrated to the DNA-methylation case study of
#' Guo et al. (2022): a \eqn{p = 1008}-mediator linear model with eleven active
#' loci whose effects mix positive and negative signs. [simulate_guo_mediation()]
#' generates data from these constants; the object is exported so the exact
#' calibration is inspectable.
#'
#' Pairing the exposure-mediator coefficients `Gamma_x` with the
#' outcome-mediator coefficients `alpha_m` gives a total indirect effect of
#' `beta_per_c1` \eqn{\approx -1.597} per unit of the signal scale \eqn{c_1}.
#' (The article reports \eqn{-1.5977}, computed from the full-precision Guo
#' coefficients; the values shipped here are those coefficients rounded to the
#' three decimals printed in the supplement, which give \eqn{-1.597}.)
#'
#' @format A list with components:
#' \describe{
#'   \item{p, s, n}{The mediator count (1008), confounder count (9, an
#'     intercept plus eight covariates), and case-study sample size (85).}
#'   \item{locations}{Indices of the eleven active mediators.}
#'   \item{alpha_m}{Outcome-mediator coefficients at the active loci (the
#'     designed simulation values).}
#'   \item{Gamma_x}{Exposure-mediator coefficients at the active loci.}
#'   \item{alpha_z}{Confounder-outcome coefficients (length 9).}
#'   \item{Gamma_z}{Confounder-mediator coefficients at the active loci (an
#'     11 by 9 matrix).}
#'   \item{beta_per_c1}{The total indirect effect per unit `c1`, equal to
#'     `sum(Gamma_x * alpha_m)` \eqn{\approx -1.597} (the article reports
#'     -1.5977 from the unrounded coefficients).}
#'   \item{alpha_m_estimated}{The raw Guo et al. (2022) estimate of the
#'     outcome-mediator coefficients (for reference; the simulation uses the
#'     designed `alpha_m`).}
#'   \item{alpha_m_variants}{Two alternative outcome-mediator coefficient sets
#'     the article also studies, `homogeneous_like` and `contrasting_like`,
#'     which remain heterogeneous when paired with `Gamma_x`.}
#' }
#'
#' @source Transcribed from the supplement of the article (its
#'   parameter-configuration section), which in turn calibrates to the
#'   `simulation_allS.Rdata` of Guo et al. (2022).
#'
#' @references
#' Guo, X., Li, R., Liu, J., & Zeng, M. (2022). High-dimensional mediation
#' analysis for selecting DNA methylation loci mediating childhood trauma and
#' cortisol stress reactivity. \emph{Journal of the American Statistical
#' Association, 117}(539), 1110--1121.
#' \doi{10.1080/01621459.2022.2053136}
#'
#' @seealso [simulate_guo_mediation()]
#'
#' @family mediation simulation
"guo_calibration"
