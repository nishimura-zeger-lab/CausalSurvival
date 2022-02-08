#' Estimate IPW of survival probability / rmst at time tau
#'
#' @param cenHaz Data frame with two columns: Haz1, Haz0
#'               Estimated censoring hazards for each person at each time points if receive treatment 1 (Haz1) and if receive treatment 0 (Haz0)
#' @param treatProb Estimated probability for treatment for each person if receive treatment 1
#' @param tau Time of interest. Can be a vector (multiple time of interest)
#' @return A data frame with three columns: estimand1, estimand0, SE

estimateIPW <- function(treatment, eventObserved, time, cenHaz, treatProb, tau, estimand, printTau){

  ## container
  estimand1_result <- estimand0_result <- SE_result <- rep(0, length=length(tau))

  ## parameter
  n <- length(treatment)
  maxTime <- dim(cenHaz)[1]/n
  ID <- rep(1:n, each=maxTime)

  ## dlong
  dlong <- transformData(dwide=data.frame(eventObserved=eventObserved, time=time, treat=treatment), freqTime=1, type="survival")
  rownames(dlong) <- NULL
  dlong <- dlong[which(dlong$t <= maxTime), c("Lt", "It", "t", "treat")]

  ## Estimate survival hazards
  survHazFit <- lm(Lt ~ as.factor(t) * as.factor(treat), subset = It == 1, data = dlong)
  survHaz1 <- predict(survHazFit, newdata = data.frame(t=1:maxTime, treat=rep(1, maxTime)), type = 'response')
  survHaz0 <- predict(survHazFit, newdata = data.frame(t=1:maxTime, treat=rep(0, maxTime)), type = 'response')

  rm(list=c("survHazFit"))

  survHaz1 <- rep(survHaz1, n)
  survHaz0 <- rep(survHaz0, n)

  ## calculate survival and censoring probability
  SurvProb1List <- tapply(1 - survHaz1, ID, cumprod, simplify = FALSE)
  SurvProb0List <- tapply(1 - survHaz0, ID, cumprod, simplify = FALSE)

  SurvProb1 <- unlist(SurvProb1List, use.names = FALSE)
  SurvProb0 <- unlist(SurvProb0List, use.names = FALSE)

  rm(list=c("SurvProb1List", "SurvProb0List"))

  CenProb1List <- tapply(1 - cenHaz$Haz1, ID, cumprod, simplify = FALSE)
  CenProb0List <- tapply(1 - cenHaz$Haz0, ID, cumprod, simplify = FALSE)

  CenProb1 <- unlist(CenProb1List, use.names = FALSE)
  CenProb0 <- unlist(CenProb0List, use.names = FALSE)

  rm(list=c("CenProb1List", "CenProb0List", "cenHaz"))

  ## Observed estimated survival hazards
  SurvHaz_obs <- treatment[ID] * survHaz1 + (1-treatment[ID]) * survHaz0

  rm(list=c("survHaz1", "survHaz0"))

  weightH1 <- 1/(SurvProb1 * treatProb[ID] * CenProb1)
  weightH0 <- 1/(SurvProb0 * (1-treatProb[ID]) * CenProb0)

  weightH1[which(weightH1 >= quantile(weightH1, probs = 0.95))] <- quantile(weightH1, probs = 0.95)
  weightH0[which(weightH0 >= quantile(weightH0, probs = 0.95))] <- quantile(weightH0, probs = 0.95)

  Z1weight <-  1 / (treatProb[ID] * CenProb1)
  Z0weight <-  1 / ((1-treatProb[ID]) * CenProb0)
  Z1weight[which(Z1weight >= quantile(Z1weight, probs = 0.95))] <- quantile(Z1weight, probs = 0.95)
  Z0weight[which(Z0weight >= quantile(Z0weight, probs = 0.95))] <- quantile(Z0weight, probs = 0.95)


  for (TimePoint in tau){

    if(estimand=="rmst"){
      if(TimePoint == 1){next}
    }

    ## parameter
    if(estimand=="rmst"){
      ind <- (dlong$t <= (TimePoint -1))
    }else if(estimand=="risk"){
      ind <- (dlong$t <= TimePoint)
    }

    ## IPW

    if(estimand=="risk"){
      Z1ipw <- - ind * Z1weight
      Z0ipw <- - ind * Z0weight
    }else if(estimand=="rmst"){
      indTillTimePoint <- unlist(tapply(ind, ID, function(x){rev(cumsum(rev(x)))}), use.names = FALSE)
      Z1ipw <- - indTillTimePoint * Z1weight
      Z0ipw <- - indTillTimePoint * Z0weight
      rm(list=c("indTillTimePoint"))
    }

    DT1ipw <- tapply(dlong$It * treatment[ID] * Z1ipw * dlong$Lt, ID, sum)
    DT0ipw <- tapply(dlong$It * (1-treatment[ID]) * Z0ipw * dlong$Lt , ID, sum)

    if(estimand=="risk"){
      ipw  <- 1 + c(mean(DT0ipw), mean(DT1ipw))
    }else if(estimand=="rmst"){
      ipw  <- TimePoint + c(mean(DT0ipw), mean(DT1ipw))
    }

    rm(list=c("DT1ipw", "DT0ipw", "Z1ipw", "Z0ipw"))

    ## variance
    if(estimand=="risk"){
      H1 <- - (ind * rep(SurvProb1[which(dlong$t == TimePoint)], each=max(dlong$t))) * weightH1
      H0 <- - (ind * rep(SurvProb0[which(dlong$t == TimePoint)], each=max(dlong$t))) * weightH0
    }else if(estimand=="rmst"){
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

    D <- DT1 - DT0 + DW1 - DW0
    sdn <- sqrt(var(D) / n)

    rm(list=c("D", "DW1", "DW0", "DT1", "DT0"))

    ## store
    estimand1_result[TimePoint] <- ipw[2]
    estimand0_result[TimePoint] <- ipw[1]
    SE_result[TimePoint] <- sdn

    if(printTau){print(paste("Time point", TimePoint, "finished"))}

  }

  ## result
  out <- data.frame(estimand1=estimand1_result, estimand0=estimand0_result, SE=SE_result)
  return(out)
}









