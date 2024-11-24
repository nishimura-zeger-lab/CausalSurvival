#' Estimate cross-fitted TMLE of survival probability at time 1, 2, ...tau
#'
#' @param survHaz Data frame with two columns: SurvHaz1, SurvHaz0
#'                Estimated survival hazards for each person at each time points if receive treatment 1 (SurvHaz1) and if receive treatment 0 (SurvHaz0)
#' @param cenHaz Data frame with two columns: CenHaz1, CenHaz0
#'               Estimated censoring hazards for each person at each time points if receive treatment 1 (CenHaz1) and if receive treatment 0 (CenHaz0)
#' @param treatProb Estimated probability for treatment for each person if receive treatment 1
#' @param tau Max time of interest. A scalar
#' @return A data frame with three columns: SurvProb1, SurvProb0, SEprobDiff

estimateTMLEprob_monotone <- function(treatment, eventObserved, time, survHaz, cenHaz, treatProb, tau){

  ## check if ID is in the same order
  identical(survHaz$ID, cenHaz$ID)

  ## container
  SurvProb1_result <- SurvProb0_result <- SEprobDiff_result <- rep(0, length=length(tau))

  ## calculate censoring probability that doesn't need iterative updates
  CenProb1List <- tapply(1 - cenHaz$Haz1, cenHaz$ID, cumprod, simplify = FALSE)
  CenProb0List <- tapply(1 - cenHaz$Haz0, cenHaz$ID, cumprod, simplify = FALSE)

  CenProb1 <- unlist(CenProb1List, use.names = FALSE)
  CenProb0 <- unlist(CenProb0List, use.names = FALSE)

  rm(list=c("CenProb1List", "CenProb0List", "cenHaz"))

  ## initial SurvHaz
  SurvHaz1 <- survHaz$Haz1
  SurvHaz0 <- survHaz$Haz0
  SurvHaz_obs  <- treatment[survHaz$ID]*SurvHaz1 + (1-treatment[survHaz$ID])*SurvHaz0

  ## parameter
  n <- length(unique(survHaz$ID))
  maxTime <- dim(survHaz)[1]/n
  converged <- FALSE
  iter <- 1

  ## dlong
  dlong <- transformData(dwide=data.frame(eventObserved=eventObserved, time=time), freqTime=1)
  dlong <- dlong[which(dlong$t <= maxTime),]

  ## iterate
  while((!converged) && iter <= 20){

    if(iter == 1){
      initial_coef <- NULL
    }else{initial_coef <- eps}

    SurvProb1List <- tapply(1 - SurvHaz1, SurvHaz$ID, cumprod, simplify = FALSE)
    SurvProb0List <- tapply(1 - SurvHaz0, SurvHaz$ID, cumprod, simplify = FALSE)

    SurvProb1 <- unlist(SurvProb1List, use.names = FALSE)
    SurvProb0 <- unlist(SurvProb0List, use.names = FALSE)

    rm(list=c("SurvProb1List", "SurvProb0List"))

    weightH1 <- 1/bound(SurvProb1 * treatProb[SurvHaz$ID] * CenProb1)
    weightH0 <- 1/bound(SurvProb0 * (1-treatProb[SurvHaz$ID]) * CenProb0)

    H1 <- H0 <- c()

    for (TimePoint in tau){

      ind <- (dlong$t <= TimePoint)
      H1_temp <- as(matrix(- (ind * rep(SurvProb1[which(dlong$t == TimePoint)], each=max(dlong$t))), nrow = 1), "sparseMatrix")
      H0_temp <- as(matrix(- (ind * rep(SurvProb0[which(dlong$t == TimePoint)], each=max(dlong$t))), nrow = 1), "sparseMatrix")

      H1 <- cbind(H1, H1_temp)
      H0 <- cbind(H0, H0_temp)

      rm(list=c("ind", "H1_temp", "H0_temp"))

    }

    H1 <- H1 * weightH1
    H0 <- H0 * weightH0
    H <-  H1 * treatment[survHaz$ID] + H0 * (1-treatment[survHaz$ID])


    ## update for survival hazard
    eps   <- glmSparse(outcome=dlong$Lt, offset=SurvHaz_obs, H=H, subset=which(dlong$It == 1),
                       maxiter=40, threshold=1e-14, initial_coef=initial_coef, printIter=FALSE)

    ## NA as 0 for the new values
    eps[is.na(eps)] <- 0

    ## update values
    SurvHaz1 <- plogis(qlogis(SurvHaz1) + H1 %*% eps)
    SurvHaz0 <- plogis(qlogis(SurvHaz0) + H0 %*% eps)
    SurvHaz_obs  <- treatment[survHaz$ID] * SurvHaz1  + (1 - treatment[survHaz$ID]) * SurvHaz0

    iter <-  iter + 1
    converged <- (abs(eps) <= 1e-3/(n^(0.6)))

    rm(list=c("H1", "H0", "H","SurvProb1", "SurvProb0"))

  }

  ## result for time tau
  SurvProb1List <- tapply(1 - SurvHaz1, survHaz$ID, cumprod, simplify = FALSE)
  SurvProb0List <- tapply(1 - SurvHaz0, survHaz$ID, cumprod, simplify = FALSE)

  SurvProb1 <- unlist(SurvProb1List, use.names = FALSE)
  SurvProb0 <- unlist(SurvProb0List, use.names = FALSE)

  rm(list=c("SurvProb1List", "SurvProb0List", "SurvHaz1", "SurvHaz0"))


  weightH1 <- 1/bound(SurvProb1 * treatProb[SurvHaz$ID] * CenProb1)
  weightH0 <- 1/bound(SurvProb0 * (1-treatProb[SurvHaz$ID]) * CenProb0)

  DT <- c()

  for (TimePoint in tau){

    ind <- (dlong$t <= TimePoint)
    H1_temp <- - (ind * rep(SurvProb1[which(dlong$t == TimePoint)], each=max(dlong$t)))
    H0_temp <- - (ind * rep(SurvProb0[which(dlong$t == TimePoint)], each=max(dlong$t)))

    DT_temp <- tapply(dlong$It * (treatment[survHaz$ID] * H1_temp - (1 - treatment[survHaz$ID]) * H0_temp) * (dlong$Lt - SurvHaz_obs), survHaz$ID, sum)
    DT_temp <- DT_temp + SurvProb1 - SurvProb0

    DT <- cbind(DT, DT_temp)

    rm(list=c("ind", "H1_temp", "H0_temp", "DT_temp"))

  }

  rm(list=c("SurvProb1", "SurvProb0"))


  ## standard error of S(1, tau)-S(0, tau)
  sdn <- sqrt(diag(var(D)) / n)

  rm(list=c("DT", "D"))

  ## store
  SurvProb1_result <- tapply(SurvProb1, dlong$t, mean)
  SurvProb0_result <- tapply(SurvProb0, dlong$t, mean)
  SEprobDiff_result <- sdn


  ## result
  out <- data.frame(SurvProb1=SurvProb1_result, SurvProb0=SurvProb0_result, SEprobDiff=SEprobDiff_result)
  return(out)
}
