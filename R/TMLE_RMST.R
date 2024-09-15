#' Estimate cross-fitted TMLE of RMST at time tau
#'
#' @param covariates Design matrix in triplet format (row index, col index, and value)
#' @param covariates.names Corresponding covariates name in ''covariates''
#' @param crossFitnum For cross-fitting: random partition of subjects into XXX prediction sets of approximately the same size.
#' @param SurvHaz.estimate Model for estimating survival hazards. Options currently include logistic LASSO, glm
#' @param CenHaz.estimate Model for estimating censoring hazards. Options currently include logistic LASSO, glm
#' @param TreatProb.estimate Model for estimating treatment probability. Options currently include logistic LASSO
#' @param maxCohortSizeForFitting If the target or comparator cohort are larger than this number, they
#'                                 will be downsampled before fitting the propensity model. The model
#'                                 will be used to compute propensity scores for all subjects. The
#'                                 purpose of the sampling is to gain speed.
#' @param covID.SurvHaz Covariates id to include in modeling discrete survival hazards. If NULL, then include all covariates in the model
#' @param covID.CenHaz Covariates id to include in modeling discrete censoring hazards. If NULL, then include all covariates in the model
#' @param covID.TreatProb Covariates id to include in modeling treatment probability. If NULL, then include all covariates in the model
#' @param freq.time Coarsen observed time to XXX days intervals
#' @param tau Time of interest. Can be a vector (multiple time of interest)
#' @export A list SurvRMST1: RMST(1, tau); SurvRMST0: RMST(0, tau); std.error.diff: standard error of SurvRMST1-SurvRMST0

estimateTMLErmst <- function(eventTime, censorTime, treatment, covariates, covariates.names,
                             crossFitnum=5, SurvHaz.estimate="LASSO", CenHaz.estimate="LASSO",
                             TreatProb.estimate="LASSO", maxCohortSizeForFitting=250000,
                             covID.SurvHaz=NULL, covID.CenHaz=NULL, covID.TreatProb=NULL,
                             freq.time=90, tau){

  ## Indicator for event
  eventObserved <- ifelse(is.na(eventTime), 0, 1)
  ## Observed time
  censored <- is.na(eventTime)
  time <- eventTime
  time[censored] <- censorTime[censored]
  ## subject id
  id <- 1:length(time)


  ## Get index for cross-fitting
  index_ls <- crossFit(eventObserved=eventObserved, id=id, crossFitnum=crossFitnum)


  ## Wide-form dataset
  dwide <- data.frame(id=id, time=time, eventObserved=eventObserved, treatment=treatment)


  ## transform data into long format
  dlong <- transformData(dwide=dwide, freq.time=freq.time)



  ## Estimate cross-fitted nuisance parameter
  SurvHaz <- estimateSurvHaz(dlong=dlong, covariates=covariates, covID.SurvHaz=covID.SurvHaz, crossFitnum=crossFitnum, SurvHaz.estimate=SurvHaz.estimate)
  CenHaz <- estimateCenHaz(dlong=dlong, covariates=covariates, covID.CenHaz=covID.CenHaz, crossFitnum=crossFitnum, index_ls=index_ls, CenHaz.estimate=CenHaz.estimate)
  TreatProb <- estimateTreatProb(id=id, treatment=treatment,
                                 covariates=covariates, covID.TreatProb=covID.TreatProb,
                                 TreatProb.estimate=TreatProb.estimate,
                                 maxCohortSizeForFitting=maxCohortSizeForFitting,
                                 index_ls=index_ls, crossFitnum=crossFitnum)


  ## into the same order
  SurvHaz <- SurvHaz[order(SurvHaz$ID), ]
  CenHaz <- CenHaz[order(CenHaz$ID), ]
  TreatProb <- TreatProb[order(TreatProb$id), ]


  ## container
  SurvRMST1_result <- SurvRMST0_result <- std.error.diff <- rep(0, length=length(tau))


  ## calculate censoring probability that doesn't need iterative updates
  CenProb1List <- tapply(1 - CenHaz$CenHaz1, dlong$id, cumprod, simplify = FALSE)
  CenProb0List <- tapply(1 - CenHaz$CenHaz0, dlong$id, cumprod, simplify = FALSE)

  CenProb1 <- unlist(CenProb1List, use.names = FALSE)
  CenProb0 <- unlist(CenProb0List, use.names = FALSE)

  rm(list=c("CenProb1List", "CenProb0List", "CenHaz"))


  for (TimePoint in tau){

    if(TimePoint == 1){next}

    ## initial SurvHaz
    SurvHaz1 <- SurvHaz$SurvHaz1
    SurvHaz0 <- SurvHaz$SurvHaz0
    SurvHaz_obs  <- dlong$treatment*SurvHaz1 + (1-dlong$treatment)*SurvHaz0

    ## parameter
    ind <- (dlong$t <= TimePoint-1)
    crit <- TRUE
    iter <- 1

    ## iterate
    while(crit && iter <= 20){

      SurvProb1List <- tapply(1 - SurvHaz1, dlong$id, cumprod, simplify = FALSE)
      SurvProb0List <- tapply(1 - SurvHaz0, dlong$id, cumprod, simplify = FALSE)

      SurvProb1 <- unlist(SurvProb1List, use.names = FALSE)
      SurvProb0 <- unlist(SurvProb0List, use.names = FALSE)

      rm(list=c("SurvProb1List", "SurvProb0List"))

      ## clever covariate for survival hazard
      H1 <- - unlist(tapply(ind * SurvProb1, dlong$id, function(x){rev(cumsum(rev(x)))}), use.names = FALSE) / bound(SurvProb1 * TreatProb$TreatProb[dlong$id] * CenProb1)
      H0 <- - unlist(tapply(ind * SurvProb0, dlong$id, function(x){rev(cumsum(rev(x)))}), use.names = FALSE) / bound(SurvProb0 * (1-TreatProb$TreatProb[dlong$id]) * CenProb0)
      H <- dlong$treatment * H1 + (1-dlong$treatment) * H0

      ## update for survival hazard
      eps   <- coef(glm2::glm2(Lt ~ 0 + offset(qlogis(SurvHaz_obs)) + H,
                               family = binomial(), subset = It == 1, data = dlong))

      ## NA as 0 for the new values
      eps[is.na(eps)] <- 0

      ## update values
      SurvHaz1 <- bound01(plogis(qlogis(SurvHaz1) + eps * H1))
      SurvHaz0 <- bound01(plogis(qlogis(SurvHaz0) + eps * H0))
      SurvHaz_obs  <- dlong$treatment * SurvHaz1  + (1 - dlong$treatment) * SurvHaz0

      iter <-  iter + 1
      crit <- abs(eps) > 1e-3/length(id)^(0.6)

      ## clear
      rm(list=c("H1", "H0", "H","SurvProb1", "SurvProb0"))
    }

    ## result for time tau
    SurvProb1List <- tapply(1 - SurvHaz1, dlong$id, cumprod, simplify = FALSE)
    SurvProb0List <- tapply(1 - SurvHaz0, dlong$id, cumprod, simplify = FALSE)

    SurvProb1 <- unlist(SurvProb1List, use.names = FALSE)
    SurvProb0 <- unlist(SurvProb0List, use.names = FALSE)

    rm(list=c("SurvProb1List", "SurvProb0List", "SurvHaz1", "SurvHaz0"))

    H1 <- - unlist(tapply(ind * SurvProb1, dlong$id, function(x){rev(cumsum(rev(x)))}), use.names = FALSE) / bound(SurvProb1 * TreatProb$TreatProb[dlong$id] * CenProb1)
    H0 <- - unlist(tapply(ind * SurvProb0, dlong$id, function(x){rev(cumsum(rev(x)))}), use.names = FALSE) / bound(SurvProb0 * (1-TreatProb$TreatProb[dlong$id]) * CenProb0)

    DT <- with(dlong, tapply(It * (treatment * H1 - (1 - treatment) * H0) * (Lt - SurvHaz_obs), id, sum))

    DW1 <- tapply(ind * SurvProb1, dlong$id, sum)
    DW0 <- tapply(ind * SurvProb0, dlong$id, sum)

    ## RMST(1, tau)
    SurvRMST1_mean <- mean(DW1)
    ## RMST(0, tau)
    SurvRMST0_mean <- mean(DW0)
    ## standard error of RMST(1, tau)-RMST(0, tau)
    D <- DT + DW1 - DW0
    sdn <- sqrt(var(D) / length(id))

    ## store
    SurvRMST1_result[TimePoint] <- SurvRMST1_mean
    SurvRMST0_result[TimePoint] <- SurvRMST0_mean
    std.error.diff[TimePoint] <- sdn

    ## clear workspace
    rm(list=c("SurvProb1", "SurvProb0", "H1", "H0", "DT", "DW1", "DW0", "D"))

  }

  ## result
  out <- list(SurvRMST1=SurvRMST1, SurvRMST0=SurvRMST0, std.error.diff=std.error.diff)
  return(out)
}





