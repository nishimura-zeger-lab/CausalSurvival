#' Estimate TMLE of survival probability / rmst at time tau
#'
#' @param survHaz Data frame with two columns: Haz1, Haz0
#'                Estimated survival hazards for each person at each time points if receive treatment 1 (SurvHaz1) and if receive treatment 0 (SurvHaz0)
#' @param cenHaz Data frame with two columns: Haz1, Haz0
#'               Estimated censoring hazards for each person at each time points if receive treatment 1 (CenHaz1) and if receive treatment 0 (CenHaz0)
#' @param treatProb Estimated probability for treatment for each person if receive treatment 1
#' @param tau Time of interest. Can be a vector (multiple time of interest)
#' @param estimand "risk", "rmst", "both"
#' @return A data frame with three columns: estimand1, estimand0, SE

estimateTMLE <- function(treatment, eventObserved, time,
                         survHaz, cenHaz, treatProb, tau,
                         timeIntMidPoint, timeIntLength,
                         estimand, printIter, printTau){

  ## container
  estimand1_result <- estimand0_result <- SE_result <- rep(0, length=length(tau))
  if(estimand == "both"){
    Dac <- 0
    SE_transform_result <- rep(0, length=length(tau))
  }

  ## calculate censoring probability that doesn't need iterative updates
  prodlag <- function(x) cumprod(c(1, x[-length(x)]))
  CenProb1List <- tapply(1 - cenHaz$Haz1, survHaz$ID, prodlag, simplify = FALSE)
  CenProb0List <- tapply(1 - cenHaz$Haz0, survHaz$ID, prodlag, simplify = FALSE)

  CenProb1 <- unlist(CenProb1List, use.names = FALSE)
  CenProb0 <- unlist(CenProb0List, use.names = FALSE)

  rm(list=c("CenProb1List", "CenProb0List", "cenHaz"))

  ## parameter
  n <- length(unique(survHaz$ID))
  maxTime <- dim(survHaz)[1]/n

  ## dlong
  dlong <- transformData(time=time, eventObserved=eventObserved,
                         timeIntMidPoint=timeIntMidPoint, type="survival")
  rownames(dlong) <- NULL

  for (TimePoint in 1:length(tau)){

    ## initial SurvHaz
    SurvHaz1 <- survHaz$Haz1
    SurvHaz0 <- survHaz$Haz0
    SurvHaz_obs  <- treatment[survHaz$ID]*SurvHaz1 + (1-treatment[survHaz$ID])*SurvHaz0

    ## parameter
    ind <- (dlong$time <= tau[TimePoint])
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

      if(estimand %in% c("risk", "both")){

      ## clever covariate for updating survival hazards
      H1 <- as(matrix(- (ind * rep(SurvProb1[which(dlong$time == tau[TimePoint])], each=maxTime)) * weightH1, ncol = 1), "sparseMatrix")
      H0 <- as(matrix(- (ind * rep(SurvProb1[which(dlong$time == tau[TimePoint])], each=maxTime)) * weightH0, ncol = 1), "sparseMatrix")

      rm(list=c("weightH1", "weightH0"))

      }else if(estimand == "rmst"){

      cumProb1TillTimePoint <- unlist(tapply(ind * SurvProb1, survHaz$ID, function(x){rev(cumsum(rev(x)))}), use.names = FALSE)
      H1 <- as(matrix(- cumProb1TillTimePoint * weightH1 * rep(timeIntLength, n), ncol = 1), "sparseMatrix")

      cumProb0TillTimePoint <- unlist(tapply(ind * SurvProb0, survHaz$ID, function(x){rev(cumsum(rev(x)))}), use.names = FALSE)
      H0 <- as(matrix(- cumProb0TillTimePoint * weightH0 * rep(timeIntLength, n), ncol = 1), "sparseMatrix")

      rm(list=c("cumProb1TillTimePoint", "cumProb0TillTimePoint", "weightH1", "weightH0"))

      }

      H <- treatment[survHaz$ID[which(dlong$valid == 1 & ind)]] * H1[which(dlong$valid == 1 & ind)] + (1-treatment[survHaz$ID[which(dlong$valid == 1 & ind)]]) * H0[which(dlong$valid == 1 & ind)]
      H <- as(matrix(H, ncol = 1), "sparseMatrix")

      ## update for survival hazard
      # eps   <- coef(glm2::glm2(dlong$Lt[which(dlong$It == 1 & ind)] ~ 0 + offset(qlogis(SurvHaz_obs[which(dlong$It == 1 & ind)])) + H,
      #                          family = binomial()))

      if(iter == 1){
        eps   <- coefSparse(outcome=dlong$y[which(dlong$valid == 1 & ind)], offset=qlogis(SurvHaz_obs[which(dlong$valid == 1 & ind)]), H=H,
                            maxiter=40, threshold=1e-8, initial_coef=NULL, printIter=FALSE)
      }else{
        eps   <- coefSparse(outcome=dlong$y[which(dlong$valid == 1 & ind)], offset=qlogis(SurvHaz_obs[which(dlong$valid == 1 & ind)]), H=H,
                            maxiter=40, threshold=1e-8, initial_coef=eps, printIter=FALSE)
      }

      ## NA as 0 for the new values
      eps[is.na(eps)] <- 0

      ## update values
      SurvHaz1 <- plogis(qlogis(SurvHaz1) + eps * H1[, 1])
      SurvHaz0 <- plogis(qlogis(SurvHaz0) + eps * H0[, 1])
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


    if(estimand %in% c("risk", "both")){

      ## clever covariate for updating survival hazards
      H1 <- as(matrix(- (ind * rep(SurvProb1[which(dlong$time == tau[TimePoint])], each=maxTime)) * weightH1, ncol = 1), "sparseMatrix")
      H0 <- as(matrix(- (ind * rep(SurvProb1[which(dlong$time == tau[TimePoint])], each=maxTime)) * weightH0, ncol = 1), "sparseMatrix")

      rm(list=c("weightH1", "weightH0"))

    }else if(estimand == "rmst"){

      cumProb1TillTimePoint <- unlist(tapply(ind * SurvProb1, survHaz$ID, function(x){rev(cumsum(rev(x)))}), use.names = FALSE)
      H1 <- as(matrix(- cumProb1TillTimePoint * weightH1 * rep(timeIntLength, n), ncol = 1), "sparseMatrix")

      cumProb0TillTimePoint <- unlist(tapply(ind * SurvProb0, survHaz$ID, function(x){rev(cumsum(rev(x)))}), use.names = FALSE)
      H0 <- as(matrix(- cumProb0TillTimePoint * weightH0 * rep(timeIntLength, n), ncol = 1), "sparseMatrix")

      rm(list=c("cumProb1TillTimePoint", "cumProb0TillTimePoint", "weightH1", "weightH0"))

    }

    DT <- tapply(dlong$valid * (treatment[survHaz$ID] * H1[, 1] - (1 - treatment[survHaz$ID]) * H0[, 1]) * (dlong$y - SurvHaz_obs), survHaz$ID, sum)

    rm(list=c("H1", "H0"))

    if(estimand %in% c("risk", "both")){

      DW1 <- SurvProb1[which(dlong$time == tau[TimePoint])]
      DW0 <- SurvProb0[which(dlong$time == tau[TimePoint])]

    }else if(estimand == "rmst"){

      DW1 <- tapply(ind * SurvProb1 * rep(timeIntLength, n), survHaz$ID, sum)
      DW0 <- tapply(ind * SurvProb0 * rep(timeIntLength, n), survHaz$ID, sum)

    }

    rm(list=c("SurvProb1", "SurvProb0", "ind"))

    ## estimand1 and estimand0
    estimand1 <- mean(DW1)
    estimand0 <- mean(DW0)

    ## standard error of estimand1-estimand0
    D <- DT + DW1 - DW0
    sdn <- sqrt(var(D) / n)
    if(estimand == "both"){
      Dac <- Dac + D * timeIntLength[TimePoint]
      SE_transform_result[TimePoint] <- sqrt(var(Dac) / n)
    }

    rm(list=c("DT", "DW1", "DW0", "D"))

    ## store
    estimand1_result[TimePoint] <- estimand1
    estimand0_result[TimePoint] <- estimand0
    SE_result[TimePoint] <- sdn

    rm(list=c("estimand1", "estimand0", "sdn"))

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


