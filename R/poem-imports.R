# Namespace imports for POEM. The estimation internals (ported verbatim from
# the reference implementation) call stats and glmnet functions unqualified, so
# those namespaces are imported wholesale here; ncvreg contributes only its
# column standardizer. Roxygen turns these tags into the NAMESPACE import
# directives.
#' @import stats
#' @import glmnet
#' @importFrom ncvreg std
NULL

# WHO_mediation_design() takes the shipped data set as a default argument, and
# simulate_guo_mediation() reads the shipped guo_calibration object, so their
# names are "global variables" from R CMD check's point of view. Declaring them
# here keeps the check clean without changing behavior.
utils::globalVariables(c("WHO_health_mediation", "guo_calibration"))
