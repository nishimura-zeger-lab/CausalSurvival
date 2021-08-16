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

  ## Observed time
  time <- eventTime
  time[which(is.na(time))] <- censorTime[which(is.na(time))]
  ## subject ID
  rowId <- 1:length(eventTime)

  ## stratify by PS score
  data <- data.frame(rowId = rowId, treatment = treatment, propensityScore = pscore)
  strata <- CohortMethod::stratifyByPs(data, 5)

  ## dataset
  strata$time <- time
  strata$y <- ifelse(is.na(eventTime), 0, 1)
  strata$treatment <- as.factor(treatment)
  strata$stratumId <- as.factor(strata$stratumId)

  ## model
  if(timeVariedTreatmentEffect == FALSE & timeVariedPscoreEffect == FALSE & pscoreAdjustment == 'stratification'){
    fit <- timereg::aalen(formula = Surv(time, y) ~ const(treatment)+const(stratumId), data=strata, n.sim=1000)
  }else if(timeVariedTreatmentEffect == FALSE & timeVariedPscoreEffect == TRUE & pscoreAdjustment == 'stratification'){
    fit <- timereg::aalen(formula = Surv(time, y) ~ const(treatment)+stratumId, data=strata, n.sim=1000)
  }else if(timeVariedTreatmentEffect == TRUE & timeVariedPscoreEffect == TRUE & pscoreAdjustment == 'stratification'){
    fit <- timereg::aalen(formula = Surv(time, y) ~ treatment+stratumId, data=strata, n.sim=1000)
  }else if(timeVariedTreatmentEffect == TRUE & timeVariedPscoreEffect == FALSE & pscoreAdjustment == 'stratification'){
    fit <- timereg::aalen(formula = Surv(time, y) ~ treatment+const(stratumId), data=strata, n.sim=1000)
  }

  # Return average treatment effect
  if(timeVariedTreatmentEffect == FALSE){
    out <- list(ATE=fit$gamma, sd=fit$robvar.gamma)
  }else{
    out <- list(ATE=fit$cum, sd=fit$robvar.cum)
  }
  return(out)
}

estimateDoublyRobustAdditiveHazard <- function() {
  stop("The function is not yet implemented")
}
