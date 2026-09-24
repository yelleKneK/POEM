# Build the `guo_calibration` data object: the constants for the article's
# real-data-motivated heterogeneous simulation setting, transcribed from the
# supplement's "Parameter Configurations of the Real-Data Motivated Simulation
# Studies" section. These calibrate a p = 1008, s = 9 linear mediation model to
# the DNA-methylation case study of Guo et al. (2022): 11 active mediators with
# a mix of positive and negative effects. The exposure-mediator coefficients
# Gamma_x paired with the outcome-mediator coefficients alpha_m give a total
# indirect effect of beta = Gamma_x' alpha_m = -1.597 (per unit c1; the
# article reports -1.5977, computed from the unrounded Guo coefficients).
#
# Run from the package root:  Rscript data-raw/guo_calibration.R

# Outcome-mediator coefficients used in the SIMULATION at the 11 active loci
# (the designed values; the article's "setting 0"). Two alternative sets
# ("setting 12") swap these for homogeneous- and contrasting-looking values
# that, paired with this Gamma_x, are still heterogeneous.
alpha_m <- c(1.0, 0.9, 0.8, -0.9, -0.8, -0.7, 0.6, 0.5, 0.4, 0.3, 0.2)

# Exposure-mediator coefficients (Gamma_x^G) at the 11 active loci.
Gamma_x <- c(-0.251, -0.221, -0.233, 0.251, 0.295, 0.282,
             -0.332, -0.359, 0.335, -0.345, 0.234)

# Confounder-outcome coefficients (alpha_z^G), length s = 9 (an intercept term
# plus eight clinical covariates).
alpha_z <- c(-0.336, -0.070, 0.665, 0.278, 0.315, 0.201, 0.173, 0.510, 0.315)

# Confounder-mediator coefficients (Gamma_z^G) at the 11 active loci: an
# 11 (locus) by 9 (confounder) matrix.
Gamma_z <- matrix(c(
  -0.045,  0.076,  0.089,  0.127, -0.408, -0.233, -0.104, -0.442, -0.242,
  -0.197,  0.100,  0.390, -0.357, -0.310, -0.195, -0.270, -0.497, -0.419,
  -0.076,  0.149,  0.151, -0.590, -0.813, -0.242, -0.273, -1.217, -0.614,
  -0.052,  0.048,  0.103,  0.115,  0.065, -0.065,  0.114,  0.199,  0.003,
  -0.033, -0.184,  0.065,  0.156,  0.008,  0.116,  0.070,  0.365,  0.287,
  -0.012, -0.095,  0.023,  0.112,  0.200, -0.147,  0.105,  0.240, -0.080,
   0.203,  0.396, -0.400,  0.377,  0.674,  0.543,  0.452,  1.180,  0.671,
   0.017,  0.234, -0.034,  0.167, -0.411, -0.063,  0.039, -0.260, -0.091,
  -0.364,  0.082,  0.719, -0.252, -0.517,  0.007, -0.066, -0.793, -0.349,
  -0.001,  0.098,  0.001,  0.040, -0.052,  0.064, -0.040, -0.076, -0.070,
  -0.099,  0.022,  0.196, -0.103,  0.468,  0.138,  0.213,  0.491,  0.122),
  nrow = 11, ncol = 9, byrow = TRUE)

# The raw Guo et al. (2022) estimate of alpha_m (kept for reference; the
# simulation uses the designed `alpha_m` above, not these).
alpha_m_estimated <- c(0.166, 0.243, 0.248, -0.049, -0.294, -0.187,
                       0.148, 0.087, 0.112, 0.223, 0.165)

# Sanity check: the documented total indirect effect per unit c1.
beta_per_c1 <- sum(Gamma_x * alpha_m)
stopifnot(abs(beta_per_c1 - (-1.5977)) < 5e-3)
cat(sprintf("beta = Gamma_x' alpha_m = %.4f (article: -1.5977)\n", beta_per_c1))

guo_calibration <- list(
  p = 1008L, s = 9L, n = 85L,
  locations = 1:11,
  alpha_m = alpha_m,
  Gamma_x = Gamma_x,
  alpha_z = alpha_z,
  Gamma_z = Gamma_z,
  beta_per_c1 = beta_per_c1,
  alpha_m_estimated = alpha_m_estimated,
  alpha_m_variants = list(
    homogeneous_like = c(1, 0.8, 0.6, 0.4, 0.2, 0, 0, 0, 0, 0, 0),
    contrasting_like = c(1, -0.5, 0.4, -0.3, 0, 0, 0, 0, 0, 0, 0)))

dir.create("data", showWarnings = FALSE)
save(guo_calibration, file = "data/guo_calibration.rda",
     compress = "xz", version = 2)
cat("wrote data/guo_calibration.rda\n")
