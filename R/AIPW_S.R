#' Estimate cross-fitted Augmented IPW of survival probability at time tau
#'
#' @param covariants Design matrix in triplet format (row index, col index, and value)
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
#' @return A list SurvProb1: S(1, tau); SurvProb0: S(0, tau); std.error.diff: standard error of SurvProb1-SurvProb0

estimateAIPWprob <- function(eventTime, censorTime, treatment, covariates, covariates.names,
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
  SurvProb1_result <- SurvProb0_result <- std.error.diff <- rep(0, length=length(tau))

  ## calculate survival and censoring probability
  SurvProb1List <- tapply(1 - SurvHaz$SurvHaz1, dlong$id, cumprod, simplify = FALSE)
  SurvProb0List <- tapply(1 - SurvHaz$SurvHaz0, dlong$id, cumprod, simplify = FALSE)

  SurvProb1 <- unlist(SurvProb1List, use.names = FALSE)
  SurvProb0 <- unlist(SurvProb0List, use.names = FALSE)

  rm(list=c("SurvProb1List", "SurvProb0List", "SurvHaz"))

  CenProb1List <- tapply(1 - CenHaz$CenHaz1, dlong$id, cumprod, simplify = FALSE)
  CenProb0List <- tapply(1 - CenHaz$CenHaz0, dlong$id, cumprod, simplify = FALSE)

  CenProb1 <- unlist(CenProb1List, use.names = FALSE)
  CenProb0 <- unlist(CenProb0List, use.names = FALSE)

  rm(list=c("CenProb1List", "CenProb0List", "CenHaz"))

  SurvHaz_obs <- dlong$treatment*SurvHaz$SurvHaz1 + (1-dlong$treatment)*SurvHaz$SurvHaz0

  for (TimePoint in tau){

    ## parameter
    ind <- (dlong$t <= TimePoint)

    ## solve estimating equation
    H1 <- - (ind * rep(SurvProb1[which(dlong$t == TimePoint)], each=max(dlong$t))) / bound(SurvProb1 * TreatProb$TreatProb[dlong$id] * CenProb1)
    H0 <- - (ind * rep(SurvProb0[which(dlong$t == TimePoint)], each=max(dlong$t))) / bound(SurvProb0 * (1-TreatProb$TreatProb[dlong$id]) * CenProb0)

    DT1 <- with(dlong, tapply(It * treatment * H1 * (Lt - SurvHaz_obs), id, sum))
    DT0 <- with(dlong, tapply(It * (1 - treatment) * H0 * (Lt - SurvHaz_obs), id, sum))

    DW1 <- SurvProb1[which(dlong$t == TimePoint)]
    DW0 <- SurvProb0[which(dlong$t == TimePoint)]

    aipw <- c(mean(DT0 + DW0), mean(DT1 + DW1))

    D <- DT1 - DT0 + DW1 - DW0
    sdn <- sqrt(var(D) / n)

    ## store
    SurvProb1_result[TimePoint] <- aipw[2]
    SurvProb0_result[TimePoint] <- aipw[1]
    std.error.diff[TimePoint] <- sdn

    rm(list=c("DW1", "DW0", "DT1", "DT0", "H1", "H0", "ind"))

  }

  ## result
  out <- list(SurvProb1=SurvProb1_result, SurvProb0=SurvProb0_result, std.error.diff=std.error.diff)
  return(out)
}





