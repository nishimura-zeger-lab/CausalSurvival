#' Estimate (cross-fitted) TMLE of survival probability at time 1, 2, ...tau
#'
#' @param survHaz Data frame with two columns: Haz1, Haz0
#'                Estimated survival hazards for each person at each time points if receive treatment 1 (SurvHaz1) and if receive treatment 0 (SurvHaz0)
#' @param cenHaz Data frame with two columns: Haz1, Haz0
#'               Estimated censoring hazards for each person at each time points if receive treatment 1 (CenHaz1) and if receive treatment 0 (CenHaz0)
#' @param treatProb Estimated probability for treatment for each person if receive treatment 1
#' @param tau Max time of interest. A scalar
#' @return A data frame with three columns: estimand1, estimand0, SE

estimateTMLEprob_monotone <- function(treatment, eventObserved, time,
                                      survHaz, cenHaz, treatProb, tau,
                                      timeIntMidPoint, printIter){

  ## container
  estimand1_result <- estimand0_result <- SE_result <- rep(0, length=length(tau))

  ## calculate censoring probability that doesn't need iterative updates
  CenProb1List <- tapply(1 - cenHaz$Haz1, survHaz$ID, cumprod, simplify = FALSE)
  CenProb0List <- tapply(1 - cenHaz$Haz0, survHaz$ID, cumprod, simplify = FALSE)

  CenProb1 <- unlist(CenProb1List, use.names = FALSE)
  CenProb0 <- unlist(CenProb0List, use.names = FALSE)

  rm(list=c("CenProb1List", "CenProb0List"))

  ## initial SurvHaz
  SurvHaz1 <- survHaz$Haz1
  SurvHaz0 <- survHaz$Haz0
  SurvHaz_obs  <- treatment[survHaz$ID]*SurvHaz1 + (1-treatment[survHaz$ID])*SurvHaz0

  ## parameter
  n <- length(unique(survHaz$ID))
  maxTime <- length(timeIntMidPoint)
  converged <- FALSE
  iter <- 1

  ## dlong
  dlong <- transformData(dwide=data.frame(eventObserved=eventObserved, time=time), timeIntMidPoint=timeIntMidPoint, type="survival")
  rownames(dlong) <- NULL
  dlong <- dlong[, c("Lt", "It", "t")]

  ## iterate
  while((!converged) && iter <= 50){

    SurvProb1List <- tapply(1 - SurvHaz1, survHaz$ID, cumprod, simplify = FALSE)
    SurvProb0List <- tapply(1 - SurvHaz0, survHaz$ID, cumprod, simplify = FALSE)

    SurvProb1 <- unlist(SurvProb1List, use.names = FALSE)
    SurvProb0 <- unlist(SurvProb0List, use.names = FALSE)

    rm(list=c("SurvProb1List", "SurvProb0List"))

    weightH1 <- 1/(SurvProb1 * treatProb[survHaz$ID] * CenProb1)
    weightH0 <- 1/(SurvProb0 * (1-treatProb[survHaz$ID]) * CenProb0)

    weightH1[which(weightH1 >= quantile(weightH1, probs = 0.95))] <- quantile(weightH1, probs = 0.95)
    weightH0[which(weightH0 >= quantile(weightH0, probs = 0.95))] <- quantile(weightH0, probs = 0.95)

    H1 <- H0 <- c()

    for (TimePoint in timeIntMidPoint){

      ind <- (dlong$t <= TimePoint)
      H1_temp <- as(matrix(- (ind * rep(SurvProb1[which(dlong$t == TimePoint)], each=maxTime)), ncol = 1), "sparseMatrix")
      H0_temp <- as(matrix(- (ind * rep(SurvProb0[which(dlong$t == TimePoint)], each=maxTime)), ncol = 1), "sparseMatrix")

      H1 <- cbind(H1, H1_temp)
      H0 <- cbind(H0, H0_temp)

      rm(list=c("ind", "H1_temp", "H0_temp"))

    }

    H1 <- H1 * weightH1
    H0 <- H0 * weightH0
    H <-  H1[which(dlong$It == 1), ] * treatment[survHaz$ID[which(dlong$It == 1)]] + H0[which(dlong$It == 1), ] * (1-treatment[survHaz$ID[which(dlong$It == 1)]])


    ## update for survival hazard
    if(iter == 1){
      eps   <- coefSparse(outcome=dlong$Lt[which(dlong$It == 1)], offset=qlogis(SurvHaz_obs[which(dlong$It == 1)]), H=H,
                                              maxiter=40, threshold=1e-8, initial_coef=NULL, printIter=TRUE)
    }else{
      eps   <- coefSparse(outcome=dlong$Lt[which(dlong$It == 1)], offset=qlogis(SurvHaz_obs[which(dlong$It == 1)]), H=H,
                          maxiter=40, threshold=1e-8, initial_coef=eps, printIter=TRUE)
    }

    # if(iter == 1){
    #   eps   <- coef(glm2::glm2(dlong$Lt[which(dlong$It == 1)] ~ 0 + offset(qlogis(SurvHaz_obs[which(dlong$It == 1)])) + as.matrix(H),
    #                            family = binomial()))
    # }else{
    #   eps   <- coef(glm2::glm2(dlong$Lt[which(dlong$It == 1)] ~ 0 + offset(qlogis(SurvHaz_obs[which(dlong$It == 1)])) + as.matrix(H),
    #                          family = binomial(), start=eps))
    # }

    ## NA as 0 for the new values
    eps[is.na(eps)] <- 0


    ## update values
    SurvHaz1 <- plogis(qlogis(SurvHaz1) + (H1 %*% eps)[, 1])
    SurvHaz0 <- plogis(qlogis(SurvHaz0) + (H0 %*% eps)[, 1])
    SurvHaz1 <- bound01(SurvHaz1)
    SurvHaz0 <- bound01(SurvHaz0)
    SurvHaz_obs  <- treatment[survHaz$ID] * SurvHaz1  + (1 - treatment[survHaz$ID]) * SurvHaz0

    iter <-  iter + 1
    converged <- (max(abs(eps)) <= 1e-3/(n^(0.6)))

    rm(list=c("H1", "H0", "H","SurvProb1", "SurvProb0", "weightH1", "weightH0"))

    if(printIter){print(iter - 1)}

  }

  ## result for time tau
  SurvProb1List <- tapply(1 - SurvHaz1, survHaz$ID, cumprod, simplify = FALSE)
  SurvProb0List <- tapply(1 - SurvHaz0, survHaz$ID, cumprod, simplify = FALSE)

  SurvProb1 <- unlist(SurvProb1List, use.names = FALSE)
  SurvProb0 <- unlist(SurvProb0List, use.names = FALSE)

  rm(list=c("SurvProb1List", "SurvProb0List", "SurvHaz1", "SurvHaz0"))


  weightH1 <- 1/(SurvProb1 * treatProb[survHaz$ID] * CenProb1)
  weightH0 <- 1/(SurvProb0 * (1-treatProb[survHaz$ID]) * CenProb0)

  weightH1[which(weightH1 >= quantile(weightH1, probs = 0.95))] <- quantile(weightH1, probs = 0.95)
  weightH0[which(weightH0 >= quantile(weightH0, probs = 0.95))] <- quantile(weightH0, probs = 0.95)

  DT <- c()

  for (TimePoint in tau){

    ind <- (dlong$t <= TimePoint)

    H1_temp <- as(matrix(- (ind * rep(SurvProb1[which(dlong$t == TimePoint)], each=max(dlong$t))), ncol = 1), "sparseMatrix")
    H0_temp <- as(matrix(- (ind * rep(SurvProb0[which(dlong$t == TimePoint)], each=max(dlong$t))), ncol = 1), "sparseMatrix")

    DT_temp <- tapply(dlong$It * (treatment[survHaz$ID] * H1_temp - (1 - treatment[survHaz$ID]) * H0_temp) * (dlong$Lt - SurvHaz_obs), survHaz$ID, sum)
    DT_temp <- DT_temp + SurvProb1[which(dlong$t == TimePoint)] - SurvProb0[which(dlong$t == TimePoint)]

    DT <- cbind(DT, DT_temp)

    rm(list=c("ind", "H1_temp", "H0_temp", "DT_temp"))

  }


  ## standard error
  sdn <- sqrt(diag(var(DT)) / n)
  rm(list=c("DT"))

  ## store
  estimand1_result <- tapply(SurvProb1, dlong$t, mean)
  estimand0_result <- tapply(SurvProb0, dlong$t, mean)
  SE_result <- sdn

  rm(list=c("SurvProb1", "SurvProb0", "sdn", "dlong"))

  ## result
  out <- data.frame(estimand1=estimand1_result, estimand0=estimand0_result, SE=SE_result)
  return(out)
}



