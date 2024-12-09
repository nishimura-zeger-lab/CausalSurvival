#' Estimate (cross-fitted) Augmented IPW of survival probability / rmst at time tau
#'
#' @param survHaz Data frame with two columns: Haz1, Haz0
#'                Estimated survival hazards for each person at each time points if receive treatment 1 (Haz1) and if receive treatment 0 (Haz0)
#' @param cenHaz Data frame with two columns: Haz1, Haz0
#'               Estimated censoring hazards for each person at each time points if receive treatment 1 (Haz1) and if receive treatment 0 (Haz0)
#' @param treatProb Estimated probability for treatment for each person if receive treatment 1
#' @param tau Time of interest. Can be a vector (multiple time of interest)
#' @return A data frame with three columns: estimand1, estimand0, SE

estimateAIPW <- function(treatment, eventObserved, time, survHaz, cenHaz, treatProb, tau, estimand, printTau){


  ## container
  estimand1_result <- estimand0_result <- SE_result <- rep(0, length=length(tau))

  ## parameter
  n <- length(unique(survHaz$ID))
  maxTime <- dim(survHaz)[1]/n
  ID <- survHaz$ID

  ## calculate survival probability
  SurvProb1List <- tapply(1 - survHaz$Haz1, ID, cumprod, simplify = FALSE)
  SurvProb0List <- tapply(1 - survHaz$Haz0, ID, cumprod, simplify = FALSE)

  SurvProb1 <- unlist(SurvProb1List, use.names = FALSE)
  SurvProb0 <- unlist(SurvProb0List, use.names = FALSE)

  rm(list=c("SurvProb1List", "SurvProb0List"))

  ## Observed estimated survival hazards
  SurvHaz_obs <- treatProb[ID] * survHaz$Haz1 + (1-treatProb[ID]) * survHaz$Haz0
  rm(list=c("survHaz"))

  ## calculate censoring probability
  CenProb1List <- tapply(1 - cenHaz$Haz1, ID, cumprod, simplify = FALSE)
  CenProb0List <- tapply(1 - cenHaz$Haz0, ID, cumprod, simplify = FALSE)

  CenProb1 <- unlist(CenProb1List, use.names = FALSE)
  CenProb0 <- unlist(CenProb0List, use.names = FALSE)

  rm(list=c("CenProb1List", "CenProb0List", "cenHaz"))

  ## dlong
  dlong <- transformData(dwide=data.frame(eventObserved=eventObserved, time=time), freqTime=1)
  rownames(dlong) <- NULL
  dlong <- dlong[which(dlong$t <= maxTime), c("Lt", "It", "t")]

  ## denominator
  weightH1 <- 1/(SurvProb1 * treatProb[ID] * CenProb1)
  weightH0 <- 1/(SurvProb0 * (1-treatProb[ID]) * CenProb0)

  weightH1[which(weightH1 >= quantile(weightH1, probs = 0.95))] <- quantile(weightH1, probs = 0.95)
  weightH0[which(weightH0 >= quantile(weightH0, probs = 0.95))] <- quantile(weightH0, probs = 0.95)


  for (TimePoint in tau){

    if(estimand=="rmst"){
      if(TimePoint == 1){next}
    }

    ## parameter
    if(estimand=="rmst"){
      ind <- (dlong$t <= TimePoint -1)
    }else if(estimand=="risk"){
      ind <- (dlong$t <= TimePoint)
    }

    ## solve estimating equation
    if(estimand=="risk"){

    H1 <- - (ind * rep(SurvProb1[which(dlong$t == TimePoint)], each=max(dlong$t))) * weightH1
    H0 <- - (ind * rep(SurvProb0[which(dlong$t == TimePoint)], each=max(dlong$t))) * weightH0

    }else if(estimand=="rmst"){

      ## solve estimating equation
      cumProb1TillTimePoint <- unlist(tapply(ind * SurvProb1, ID, function(x){rev(cumsum(rev(x)))}), use.names = FALSE)
      H1 <- - cumProb1TillTimePoint * weightH1
      cumProb0TillTimePoint <- unlist(tapply(ind * SurvProb0, ID, function(x){rev(cumsum(rev(x)))}), use.names = FALSE)
      H0 <- - cumProb0TillTimePoint * weightH0

      rm(list=c("cumProb1TillTimePoint", "cumProb0TillTimePoint"))

    }

    DT1 <- tapply(dlong$It * treatment[ID] * H1 * (dlong$Lt - SurvHaz_obs), ID, sum)
    DT0 <- tapply(dlong$It * (1 - treatment[ID]) * H0 * (dlong$Lt - SurvHaz_obs), ID, sum)

    rm(list=c("H1", "H0"))

    if(estimand=="risk"){

      DW1 <- SurvProb1[which(dlong$t == TimePoint)]
      DW0 <- SurvProb0[which(dlong$t == TimePoint)]

    }else if(estimand=="rmst"){

      DW1 <- tapply(ind * SurvProb1, ID, sum)
      DW0 <- tapply(ind * SurvProb0, ID, sum)

    }

    ## AIPW
    aipw <- c(mean(DT0 + DW0), mean(DT1 + DW1))
    ## SE
    D <- DT1 - DT0 + DW1 - DW0
    sdn <- sqrt(var(D) / n)

    rm(list=c("D", "DW1", "DW0", "DT1", "DT0"))

    ## store
    estimand1_result[TimePoint] <- aipw[2]
    estimand0_result[TimePoint] <- aipw[1]
    SE_result[TimePoint] <- sdn

    if(printTau){print(paste("Time point", TimePoint, "finished"))}

  }

  ## result
  out <- data.frame(estimand1=estimand1_result, estimand0=estimand0_result, SE=SE_result)
  return(out)
}





