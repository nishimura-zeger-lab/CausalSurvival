#' Estimate (cross-fitted) TMLE of survival probability / rmst at time tau
#'
#' @param survHaz Data frame with two columns: Haz1, Haz0
#'                Estimated survival hazards for each person at each time points if receive treatment 1 (SurvHaz1) and if receive treatment 0 (SurvHaz0)
#' @param cenHaz Data frame with two columns: Haz1, Haz0
#'               Estimated censoring hazards for each person at each time points if receive treatment 1 (CenHaz1) and if receive treatment 0 (CenHaz0)
#' @param treatProb Estimated probability for treatment for each person if receive treatment 1
#' @param tau Time of interest. Can be a vector (multiple time of interest)
#' @param estimand "risk", "rmst"
#' @return A data frame with three columns: estimand1, estimand0, SE

estimateTMLE <- function(treatment, eventObserved, time, survHaz, cenHaz, treatProb, tau, estimand, printIter, printTau){

  ## container
  estimand1_result <- estimand0_result <- SE_result <- rep(0, length=length(tau))

  ## calculate censoring probability that doesn't need iterative updates
  CenProb1List <- tapply(1 - cenHaz$Haz1, survHaz$ID, cumprod, simplify = FALSE)
  CenProb0List <- tapply(1 - cenHaz$Haz0, survHaz$ID, cumprod, simplify = FALSE)

  CenProb1 <- unlist(CenProb1List, use.names = FALSE)
  CenProb0 <- unlist(CenProb0List, use.names = FALSE)

  rm(list=c("CenProb1List", "CenProb0List", "cenHaz"))

  ## parameter
  n <- length(unique(survHaz$ID))
  maxTime <- dim(survHaz)[1]/n

  ## dlong
  dlong <- transformData(dwide=data.frame(eventObserved=eventObserved, time=time), freqTime=1, type="survival")
  rownames(dlong) <- NULL
  dlong <- dlong[which(dlong$t <= maxTime), c("Lt", "It", "t")]

  for (TimePoint in tau){

    if(estimand=="rmst"){
      if(TimePoint == 1){next}
    }

    ## initial SurvHaz
    SurvHaz1 <- survHaz$Haz1
    SurvHaz0 <- survHaz$Haz0
    SurvHaz_obs  <- treatment[survHaz$ID]*SurvHaz1 + (1-treatment[survHaz$ID])*SurvHaz0

    ## parameter
    if(estimand=="rmst"){
      ind <- (dlong$t <= TimePoint -1)
    }else if(estimand=="risk"){
      ind <- (dlong$t <= TimePoint)
    }
    converged <- FALSE
    iter <- 1

    ## iterate
    while((!converged) && iter <= 20){

      SurvProb1List <- tapply(1 - SurvHaz1, survHaz$ID, cumprod, simplify = FALSE)
      SurvProb0List <- tapply(1 - SurvHaz0, survHaz$ID, cumprod, simplify = FALSE)

      SurvProb1 <- unlist(SurvProb1List, use.names = FALSE)
      SurvProb0 <- unlist(SurvProb0List, use.names = FALSE)

      rm(list=c("SurvProb1List", "SurvProb0List"))

      weightH1 <- 1/(SurvProb1 * treatProb[survHaz$ID] * CenProb1)
      weightH0 <- 1/(SurvProb0 * (1-treatProb[survHaz$ID]) * CenProb0)

      weightH1[which(weightH1 >= quantile(weightH1, probs = 0.95))] <- quantile(weightH1, probs = 0.95)
      weightH0[which(weightH0 >= quantile(weightH0, probs = 0.95))] <- quantile(weightH0, probs = 0.95)

      if(estimand=="risk"){

      ## clever covariate for updating survival hazards
      H1 <- - (ind * rep(SurvProb1[which(dlong$t == TimePoint)], each=max(dlong$t))) * weightH1
      H0 <- - (ind * rep(SurvProb0[which(dlong$t == TimePoint)], each=max(dlong$t))) * weightH0

      rm(list=c("weightH1", "weightH0"))

      }else if(estimand=="rmst"){

      cumProb1TillTimePoint <- unlist(tapply(ind * SurvProb1, survHaz$ID, function(x){rev(cumsum(rev(x)))}), use.names = FALSE)
      H1 <- - cumProb1TillTimePoint * weightH1

      cumProb0TillTimePoint <- unlist(tapply(ind * SurvProb0, survHaz$ID, function(x){rev(cumsum(rev(x)))}), use.names = FALSE)
      H0 <- - cumProb0TillTimePoint * weightH0

      rm(list=c("cumProb1TillTimePoint", "cumProb0TillTimePoint", "weightH1", "weightH0"))

      }

      H <- treatment[survHaz$ID[which(dlong$It == 1)]] * H1[which(dlong$It == 1)] + (1-treatment[survHaz$ID[which(dlong$It == 1)]]) * H0[which(dlong$It == 1)]

      ## update for survival hazard
      eps   <- coef(glm2::glm2(dlong$Lt[which(dlong$It == 1)] ~ 0 + offset(qlogis(SurvHaz_obs[which(dlong$It == 1)])) + H,
                               family = binomial()))

      ## NA as 0 for the new values
      eps[is.na(eps)] <- 0

      ## update values
      SurvHaz1 <- plogis(qlogis(SurvHaz1) + eps * H1)
      SurvHaz0 <- plogis(qlogis(SurvHaz0) + eps * H0)
      SurvHaz1 <- bound01(SurvHaz1)
      SurvHaz0 <- bound01(SurvHaz0)
      SurvHaz_obs  <- treatment[survHaz$ID] * SurvHaz1  + (1 - treatment[survHaz$ID]) * SurvHaz0

      iter <-  iter + 1
      converged <- (abs(eps) <= 1e-3/n^(0.6))

      if(printIter){print(iter - 1)}

      ## clear
      rm(list=c("H1", "H0", "H","SurvProb1", "SurvProb0"))

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


    if(estimand=="risk"){

      ## clever covariate for updating survival hazards
      H1 <- - (ind * rep(SurvProb1[which(dlong$t == TimePoint)], each=max(dlong$t))) * weightH1
      H0 <- - (ind * rep(SurvProb0[which(dlong$t == TimePoint)], each=max(dlong$t))) * weightH0

      rm(list=c("weightH1", "weightH0"))

    }else if(estimand=="rmst"){

      cumProb1TillTimePoint <- unlist(tapply(ind * SurvProb1, survHaz$ID, function(x){rev(cumsum(rev(x)))}), use.names = FALSE)
      H1 <- - cumProb1TillTimePoint * weightH1

      cumProb0TillTimePoint <- unlist(tapply(ind * SurvProb0, survHaz$ID, function(x){rev(cumsum(rev(x)))}), use.names = FALSE)
      H0 <- - cumProb0TillTimePoint * weightH0

      rm(list=c("cumProb1TillTimePoint", "cumProb0TillTimePoint", "weightH1", "weightH0"))

    }

    DT <- tapply(dlong$It * (treatment[survHaz$ID] * H1 - (1 - treatment[survHaz$ID]) * H0) * (dlong$Lt - SurvHaz_obs), survHaz$ID, sum)

    rm(list=c("H1", "H0"))

    if(estimand=="risk"){

      DW1 <- SurvProb1[which(dlong$t == TimePoint)]
      DW0 <- SurvProb0[which(dlong$t == TimePoint)]

    }else if(estimand=="rmst"){

      DW1 <- tapply(ind * SurvProb1, survHaz$ID, sum)
      DW0 <- tapply(ind * SurvProb0, survHaz$ID, sum)

    }

    rm(list=c("SurvProb1", "SurvProb0", "ind"))

    ## estimand1 and estimand0
    estimand1 <- mean(DW1)
    estimand0 <- mean(DW0)

    ## standard error of estimand1-estimand0
    D <- DT + DW1 - DW0
    sdn <- sqrt(var(D) / n)

    rm(list=c("DT", "DW1", "DW0", "D"))

    ## store
    estimand1_result[TimePoint] <- estimand1
    estimand0_result[TimePoint] <- estimand0
    SE_result[TimePoint] <- sdn

    rm(list=c("estimand1", "estimand0", "sdn"))

    if(printTau){print(paste("Time point", TimePoint, "finished"))}

  }

  ## result
  out <- data.frame(estimand1=estimand1_result, estimand0=estimand0_result, SE=SE_result)
  return(out)
}


