#' Estimate IPW of survival probability / rmst at time tau
#'
#' @param cenHaz Data frame with two columns: Haz1, Haz0
#'               Estimated censoring hazards for each person at each time points if receive treatment 1 (Haz1) and if receive treatment 0 (Haz0)
#' @param treatProb Estimated probability for treatment for each person if receive treatment 1
#' @param tau Time of interest. Can be a vector (multiple time of interest)
#' @return A data frame with three columns: estimand1, estimand0, SE

estimateIPW <- function(treatment, eventObserved, time,
                        cenHaz, treatProb,
                        tau, timeIntMidPoint, timeIntLength,
                        estimand, printTau){

  ## container
  estimand1_result <- estimand0_result <- SE_result <- rep(0, length=length(tau))
  if(estimand == "both"){
    Dac <- 0
    SE_transform_result <- rep(0, length=length(tau))
  }

  ## parameter
  n <- length(treatment)
  maxTime <- dim(cenHaz)[1]/n
  ID <- rep(1:n, each=maxTime)

  ## dlong
  dlong <- transformData(dwide=data.frame(eventObserved=eventObserved, time=time, treat=treatment), timeIntMidPoint=timeIntMidPoint, type="survival")
  rownames(dlong) <- NULL

  ## Estimate survival hazards
  survHaz <- coef_lm(dlong=dlong, timeIntMidPoint=timeIntMidPoint, n=n, printIter=FALSE)
  survHaz1 <- survHaz$Haz1
  survHaz0 <- survHaz$Haz0
  rm(list=c("survHaz"))

  ## calculate survival and censoring probability
  SurvProb1List <- rep(cumprod(1 - survHaz1), n)
  SurvProb0List <- rep(cumprod(1 - survHaz0), n)

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


  for (TimePoint in 1:length(tau)){

    ## parameter
    ind <- (dlong$t <= tau[TimePoint])

    ## IPW
    if(estimand %in% c("risk", "both")){
      Z1ipw <- - ind * Z1weight
      Z0ipw <- - ind * Z0weight
    }else if(estimand == "rmst"){
      indTillTimePoint <- unlist(tapply(ind, ID, function(x){rev(cumsum(rev(x)))}), use.names = FALSE)
      Z1ipw <- - indTillTimePoint * Z1weight * rep(timeIntLength, n)
      Z0ipw <- - indTillTimePoint * Z0weight * rep(timeIntLength, n)
      rm(list=c("indTillTimePoint"))
    }

    DT1ipw <- tapply(dlong$It * treatment[ID] * Z1ipw * dlong$Lt, ID, sum)
    DT0ipw <- tapply(dlong$It * (1-treatment[ID]) * Z0ipw * dlong$Lt , ID, sum)

    if(estimand %in% c("risk", "both")){
      ipw  <- 1 + c(mean(DT0ipw), mean(DT1ipw))
    }else if(estimand=="rmst"){
      ipw  <- TimePoint + c(mean(DT0ipw), mean(DT1ipw))
    }

    rm(list=c("DT1ipw", "DT0ipw", "Z1ipw", "Z0ipw"))

    ## variance
    if(estimand %in% c("risk", "both")){

      H1 <- as(matrix( - (ind * rep(SurvProb1[which(dlong$t == tau[TimePoint])], each=maxTime)) * weightH1, ncol = 1), "sparseMatrix")
      H0 <- as(matrix( - (ind * rep(SurvProb0[which(dlong$t == tau[TimePoint])], each=maxTime)) * weightH0, ncol = 1), "sparseMatrix")

    }else if(estimand=="rmst"){

      cumProb1TillTimePoint <- unlist(tapply(ind * SurvProb1, ID, function(x){rev(cumsum(rev(x)))}), use.names = FALSE)
      H1 <- as(matrix(- cumProb1TillTimePoint * weightH1 * rep(timeIntLength, n), ncol = 1), "sparseMatrix")

      cumProb0TillTimePoint <- unlist(tapply(ind * SurvProb0, ID, function(x){rev(cumsum(rev(x)))}), use.names = FALSE)
      H0 <- as(matrix(- cumProb0TillTimePoint * weightH0 * rep(timeIntLength, n), ncol = 1), "sparseMatrix")

      rm(list=c("cumProb1TillTimePoint", "cumProb0TillTimePoint"))
    }

    DT1 <- tapply(dlong$It * treatment[ID] * H1[, 1] * (dlong$Lt - SurvHaz_obs), ID, sum)
    DT0 <- tapply(dlong$It * (1 - treatment[ID]) * H0[, 1] * (dlong$Lt - SurvHaz_obs), ID, sum)

    rm(list=c("H1", "H0"))

    if(estimand %in% c("risk", "both")){

      DW1 <- SurvProb1[which(dlong$t == tau[TimePoint])]
      DW0 <- SurvProb0[which(dlong$t == tau[TimePoint])]

    }else if(estimand=="rmst"){

      DW1 <- tapply(ind * SurvProb1 * rep(timeIntLength, n), ID, sum)
      DW0 <- tapply(ind * SurvProb0 * rep(timeIntLength, n), ID, sum)

    }

    D <- DT1 - DT0 + DW1 - DW0
    sdn <- sqrt(var(D) / n)
    if(estimand == "both"){
      Dac <- Dac + D * timeIntLength[TimePoint]
      SE_transform_result[TimePoint] <- sqrt(var(Dac) / n)
    }

    rm(list=c("D", "DW1", "DW0", "DT1", "DT0"))

    ## store
    estimand1_result[TimePoint] <- ipw[2]
    estimand0_result[TimePoint] <- ipw[1]
    SE_result[TimePoint] <- sdn

    if(printTau & (TimePoint == floor(length(tau)/2))){print("Halfway finished")}

  }

  ## result
  if(estimand == "both"){
    out <- data.frame(S1=estimand1_result, S0=estimand0_result,
                      rmst1=cumsum(estimand1_result), rmst0=cumsum(estimand0_result),
                      SE_S=SE_result, SE_rmst=SE_transform_result)
  }else{
    out <- data.frame(estimand1=estimand1_result, estimand0=estimand0_result, SE=SE_result)
  }
   return(out)
}









