#' Estimate average treatment effect under additive hazard model
#'
#' @param pscoreAdjustment Takes a value of either 'stratification' or 'spline'.
#'   'match' may be supported in the future.
#' @param timeVariedPscoreEffect Time-varying covariates effects collapse into
#'   time-varying propensity score effects. In other words, set this parameter
#'   to TRUE if you want the effects of the original confounders to be time-varying.
#' @export
estimateAdditiveHazard <- function(eventTime, censorTime, treatment, pscore,
                                   pscoreAdjustment = 'spline',
                                   adjustmentSpec = list(splineDf = 5, numStrata = 5),
                                   timeVariedPscoreEffect = FALSE,
                                   timeVariedTreatmentEffect = FALSE) {
  stop("The function is not yet implemented")
  # Return average treatment effect
}

estimateDoublyRobustAdditiveHazard <- function() {
  stop("The function is not yet implemented")
}
