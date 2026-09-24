# Canonical mediator-coefficient (alpha_m) presets per outcome and pattern,
# taken directly from the article's Monte Carlo designs. "homogeneous" means
# all active mediation effects share a sign; "contrasting" means the active
# effects cancel so the total indirect effect beta = Gamma_x alpha_m is exactly
# zero (the case where a total-indirect-effect test is powerless). The presets
# place the nonzero coefficients in the first few mediators and pad the rest
# with zeros to length p. Not exported.
#' @keywords internal
#' @noRd
.pe_alpha_m_preset <- function(outcome, pattern, p) {
  nz <- switch(paste(outcome, pattern, sep = "/"),
    "continuous/homogeneous" = c(1, 0.8, 0.6, 0.4, 0.2),
    "continuous/contrasting" = c(1, -0.5, 0.4, -0.3),
    "binary/homogeneous"     = c(3, 1.5, 0, 0, 2),
    "binary/contrasting"     = c(3, -1.5),
    "count/homogeneous"      = c(0.9, 0.8, 0, 0, 0.7),
    "count/contrasting"      = c(0, 0, 0, 0.8, -0.7),
    stop("No preset for outcome/pattern combination.", call. = FALSE))
  if (length(nz) > p)
    stop("`p` is too small for this pattern (need at least ", length(nz),
         " mediators).", call. = FALSE)
  c(nz, rep(0, p - length(nz)))
}

#' Simulate high-dimensional mediation data
#'
#' Generates a data set from the mediation models studied in the article,
#' for a continuous, binary, or count outcome, under the homogeneous or
#' contrasting (heterogeneous) mediation patterns. The mediators carry an
#' AR(1) correlation structure, a handful are truly active, and the rest
#' are null, so the data exercise exactly the situation the power-enhanced
#' tests are designed for. This is the generator behind the package's
#' example data sets and behind [pe_power_curve()].
#'
#' @details
#' One exposure column (or \eqn{q} of them) is drawn standard normal. The
#' mediators follow \eqn{M = X \Gamma_x + Z \Gamma_z + \varepsilon_m}
#' with \eqn{\varepsilon_m \sim N(0, \Sigma)}, \eqn{\Sigma} the AR(1)
#' covariance [mediation_ar1_cov()] with parameter `rho`. The
#' exposure-on-mediator coefficient is \eqn{\Gamma_x = c_1 \tau}, where
#' the base pattern \eqn{\tau} has small increasing loadings on the first
#' five mediators and random noise loadings on the rest, so `c1` scales
#' the overall \eqn{X \to M} signal (and `c1 = 0` gives the global null
#' of no mediation). The outcome is generated from its model with
#' mediator coefficients \eqn{\alpha_m} (the preset for the chosen
#' `pattern`, or a user-supplied vector) and direct effect
#' \eqn{\alpha_x = c_2}:
#' \itemize{
#'   \item continuous: \eqn{Y = M\alpha_m + X\alpha_x + Z\alpha_z +
#'     \varepsilon_y}, \eqn{\varepsilon_y \sim N(0, \sigma_y^2)};
#'   \item binary: \eqn{Y \sim \mathrm{Bernoulli}(\mathrm{logit}^{-1}(
#'     M\alpha_m + X\alpha_x + Z\alpha_z))};
#'   \item count: \eqn{Y \sim \mathrm{Poisson}(\exp(\eta))}, with the
#'     log-mean \eqn{\eta} clamped to \eqn{[-5, 5]} to guard against
#'     overflow, as in the article.
#' }
#'
#' @param n Number of observations.
#' @param p Number of candidate mediators.
#' @param outcome Outcome type: `"continuous"`, `"binary"`, or `"count"`.
#' @param pattern Mediation pattern: `"homogeneous"` (active effects share
#'   a sign) or `"contrasting"` (active effects cancel; total indirect
#'   effect zero). `"heterogeneous"` is accepted as a synonym for
#'   `"contrasting"`. Ignored when `alpha_m` is supplied directly.
#' @param c1 Scale on the exposure-on-mediator coefficients \eqn{\Gamma_x
#'   = c_1 \tau}. `c1 = 0` makes every mediator inactive (the global
#'   null). Default 1.
#' @param c2 The direct effect \eqn{\alpha_x} of each exposure on the
#'   outcome. Default 0.5.
#' @param q Number of exposures. Default 1. For `q > 1` you must supply
#'   `tau` as a \eqn{q \times p} matrix.
#' @param d Number of confounders. Default 0 (none). When positive,
#'   confounders are drawn standard normal and enter both models through
#'   `alpha_z` and `Gamma_z`.
#' @param alpha_m Optional length-\eqn{p} vector of mediator-on-outcome
#'   coefficients, overriding the `pattern` preset.
#' @param tau Optional base pattern for \eqn{\Gamma_x}. Default builds the
#'   article's \eqn{\tau = (0.1, 0.2, 0.3, 0.4, 0.5, \text{noise})} for
#'   `q = 1`; required (as a \eqn{q \times p} matrix) when `q > 1`.
#' @param rho AR(1) correlation parameter for the mediator noise. Default
#'   0.5.
#' @param sigma_y Outcome noise standard deviation for a continuous
#'   outcome. Default 0.5.
#' @param alpha_z,Gamma_z Optional confounder coefficients (length \eqn{d}
#'   vector and \eqn{d \times p} matrix). Default zero when `d > 0`.
#' @param seed Optional integer seed. When supplied it is set locally and
#'   the caller's random number generator state is restored on exit;
#'   `NULL` (default) leaves the random number generator alone.
#'
#' @return A list with the simulated data and the truth used to generate
#'   it: `X` (\eqn{n \times q}), `M` (\eqn{n \times p}), `Y` (length
#'   \eqn{n}), `Z` (\eqn{n \times d} or `NULL`), the coefficient values
#'   `alpha_m`, `Gamma_x`, `alpha_x`, the total indirect effect `beta =
#'   Gamma_x alpha_m`, and `active_mediators` (the indices that are truly
#'   active, i.e. nonzero in both paths).
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @seealso [pe_mediation()] to analyze the data, [pe_power_curve()] to
#'   sweep `c1`, [mediation_ar1_cov()] for the mediator covariance.
#'
#' @family mediation simulation
#'
#' @examples
#' set.seed(113)
#' d <- simulate_mediation_data(n = 100, p = 50, outcome = "continuous",
#'                              pattern = "contrasting")
#' dim(d$M)
#' d$active_mediators       # truly active mediators
#' d$beta                   # total indirect effect (zero under "contrasting")
#'
#' @export
simulate_mediation_data <- function(n, p,
                                     outcome = c("continuous", "binary", "count"),
                                     pattern = c("homogeneous", "contrasting",
                                                 "heterogeneous"),
                                     c1 = 1, c2 = 0.5, q = 1L, d = 0L,
                                     alpha_m = NULL, tau = NULL, rho = 0.5,
                                     sigma_y = 0.5, alpha_z = NULL,
                                     Gamma_z = NULL, seed = NULL) {
  outcome <- match.arg(outcome)
  pattern <- match.arg(pattern)
  if (pattern == "heterogeneous") pattern <- "contrasting"
  .validate_sim_design(n, p, q, d, c1, c2, rho, sigma_y)

  # Local seeding that restores the caller's RNG state on exit (see poem_seed.R).
  .poem_local_seed(seed)

  # Mediator-on-outcome coefficients: a user vector or the pattern preset.
  if (is.null(alpha_m)) {
    alpha_m <- .pe_alpha_m_preset(outcome, pattern, p)
  } else {
    if (length(alpha_m) != p)
      stop("`alpha_m` must have length `p`.", call. = FALSE)
  }

  # Base exposure-on-mediator pattern tau (q-by-p), scaled by c1 to form
  # Gamma_x. The default single-exposure tau loads 0.1..0.5 on the first five
  # mediators and N(0, 0.5^2) noise on the rest, following Yu and Kelley (in press).
  if (is.null(tau)) {
    if (q != 1L)
      stop("For `q > 1` you must supply `tau` as a q-by-p matrix.",
           call. = FALSE)
    if (p < 5L) stop("`p` must be at least 5 for the default `tau`.",
                     call. = FALSE)
    tau <- matrix(c(0.1, 0.2, 0.3, 0.4, 0.5,
                    stats::rnorm(p - 5, 0, 0.5)), nrow = 1L)
  } else {
    tau <- matrix(tau, nrow = q)
    if (ncol(tau) != p)
      stop("`tau` must have `p` columns.", call. = FALSE)
  }
  Gamma_x <- c1 * tau                          # q-by-p
  alpha_x <- rep(c2, q)                         # direct effect per exposure

  # Exposures, correlated mediator noise, and mediators.
  X <- matrix(stats::rnorm(n * q), nrow = n, ncol = q)
  Sigma_m <- mediation_ar1_cov(p, rho = rho)
  eps_m <- matrix(stats::rnorm(n * p), nrow = n, ncol = p) %*% chol(Sigma_m)
  M <- X %*% Gamma_x + eps_m

  # Optional confounders enter both the mediator and outcome models.
  Z <- NULL
  if (d > 0L) {
    Z <- matrix(stats::rnorm(n * d), nrow = n, ncol = d)
    if (is.null(Gamma_z)) Gamma_z <- matrix(0, nrow = d, ncol = p)
    if (is.null(alpha_z)) alpha_z <- rep(0, d)
    M <- M + Z %*% Gamma_z
  }

  # Linear predictor for the outcome model.
  eta <- as.numeric(M %*% alpha_m + X %*% alpha_x)
  if (d > 0L) eta <- eta + as.numeric(Z %*% alpha_z)

  Y <- switch(outcome,
    continuous = eta + stats::rnorm(n, 0, sigma_y),
    binary     = stats::rbinom(n, 1, stats::plogis(eta)),
    # Clamp the log-mean to [-5, 5] before exponentiating, as in the article,
    # so a few large linear predictors do not produce overflowing counts.
    count      = stats::rpois(n, exp(pmax(pmin(eta, 5), -5))))

  beta <- as.numeric(Gamma_x %*% alpha_m)       # total indirect effect
  # Truly active mediators: nonzero on BOTH paths (outcome and exposure).
  active <- which(alpha_m != 0 & apply(Gamma_x != 0, 2L, any))

  list(X = X, M = M, Y = Y, Z = Z, alpha_m = alpha_m, Gamma_x = Gamma_x,
       alpha_x = alpha_x, beta = beta, active_mediators = active,
       n = n, p = p, q = q, outcome = outcome, pattern = pattern)
}
