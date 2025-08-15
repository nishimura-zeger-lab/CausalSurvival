#' Estimate TMLE of survival probability at time tau, RMST at time tau, log-ratio of log-survival probability
#'
#' @param treatment
#' @param eventObserved
#' @param time
#' @param survHaz Data frame with two columns: Haz1, Haz0
#'                Estimated survival hazards for each person at each time points if receive treatment 1 (SurvHaz1) and if receive treatment 0 (SurvHaz0)
#' @param cenHaz Data frame with two columns: Haz1, Haz0
#'               Estimated censoring hazards for each person at each time points if receive treatment 1 (CenHaz1) and if receive treatment 0 (CenHaz0)
#' @param treatProb Estimated probability for treatment for each person if receive treatment 1
#' @param tau Time of interest. Can be a vector (multiple time of interest)
#' @param weights
#' @param timeIntMidPoint
#' @param timeIntLength
#' @param printIter
#' @param printTau
#' @return A data frame with three columns: estimand1, estimand0, SE

estimateTMLE <- function(treatment, eventObserved, time,
                         survHaz, cenHaz, treatProb, tau,
                         timeIntMidPoint, timeIntLength,
                         printIter, printTau){

  ## container
  S1_result <- S0_result <- SE_S_diff_result <- rep(0, length=length(tau))
  RMST1_result <- RMST0_result <- SE_RMST_diff_result <- rep(0, length=length(tau))
  rho_result <- SE_rho_result <- 0
  D_rmst <- 0
  DT_rmst <- DW1_rmst <- DW0_rmst <- 0
  D_rho <- W_rho <- rho_2term <- 0

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

      ## clever covariate for updating survival hazards
      H1 <- as(matrix(- (ind * rep(SurvProb1[which(dlong$time == tau[TimePoint])], each=maxTime)) * weightH1, ncol = 1), "sparseMatrix")
      H0 <- as(matrix(- (ind * rep(SurvProb1[which(dlong$time == tau[TimePoint])], each=maxTime)) * weightH0, ncol = 1), "sparseMatrix")
      H <- cbind(H1, H0)
      rm(list=c("weightH1", "weightH0"))

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
      SurvHaz1 <- plogis(qlogis(SurvHaz1) + eps[1] * H1[, 1])
      SurvHaz0 <- plogis(qlogis(SurvHaz0) + eps[2] * H0[, 1])
      SurvHaz1 <- bound01(SurvHaz1)
      SurvHaz0 <- bound01(SurvHaz0)
      SurvHaz_obs  <- treatment[survHaz$ID] * SurvHaz1  + (1 - treatment[survHaz$ID]) * SurvHaz0

      iter <-  iter + 1
      converged <- (max(abs(eps)) <= 1e-3/n^(0.6))

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

    ## clever covariate for updating survival hazards
    H1 <- as(matrix(- (ind * rep(SurvProb1[which(dlong$time == tau[TimePoint])], each=maxTime)) * weightH1, ncol = 1), "sparseMatrix")
    H0 <- as(matrix(- (ind * rep(SurvProb1[which(dlong$time == tau[TimePoint])], each=maxTime)) * weightH0, ncol = 1), "sparseMatrix")

    rm(list=c("weightH1", "weightH0"))

    DT1 <- tapply(dlong$valid * treatment[survHaz$ID] * H1[, 1] * (dlong$y - SurvHaz_obs), survHaz$ID, sum)
    DT0 <- tapply(dlong$valid * (1 - treatment[survHaz$ID]) * H0[, 1] * (dlong$y - SurvHaz_obs), survHaz$ID, sum)

    rm(list=c("H1", "H0"))

    DW1 <- SurvProb1[which(dlong$time == tau[TimePoint])]
    DW0 <- SurvProb0[which(dlong$time == tau[TimePoint])]

    DW1_rmst <- DW1_rmst + DW1 * timeIntLength[TimePoint]
    DW0_rmst <- DW0_rmst + DW0 * timeIntLength[TimePoint]

    rm(list=c("SurvProb1", "SurvProb0", "ind"))

    ## standard error of S1-S0
    D <- DT1 - DT0 + DW1 - DW0 + mean(DW1) - mean(DW0)
    sdn_S <- sqrt(var(D) / n)
    ## standard error of rmst1-rmst0
    D_rmst <- D_rmst + D * timeIntLength[TimePoint]
    sdn_rmst <- sqrt(var(D_rmst) / n)
    ## rho
    W_rho <- W_rho + 1/(sdn_S^2)
    D_rho <- D_rho + 1/(sdn_S^2)*((DT1+DW1+mean(DW1))/(mean(DW1)*log(mean(DW1)))-
                                    (DT0+DW0+mean(DW0))/(mean(DW0)*log(mean(DW0))))
    rho_2term <- rho_2term + 1/(sdn_S^2)*log(log(mean(DW1))/log(mean(DW0)))

    rm(list=c("DT1", "DT0", "DW1", "DW0", "D"))

    ## store results
    S1_result[TimePoint] <- mean(DW1)
    S0_result[TimePoint] <- mean(DW0)
    SE_S_diff_result[TimePoint] <- sdn_S

    RMST1_result[TimePoint] <- mean(DW1_rmst)
    RMST0_result[TimePoint] <- mean(DW0_rmst)
    SE_RMST_diff_result[TimePoint] <- sdn_rmst

    rm(list=c("DW1", "DW0", "sdn_S", "DW1_rmst", "DW0_rmst", "sdn_rmst"))

    if(printTau & (TimePoint == floor(length(tau)/2))){print("Halfway finished")}

  }

  rho_result <- 1/W_rho * rho_2term
  D_rho <- 1/W_rho*D_rho
  SE_rho_result <- sqrt(var(D_rho)/n)

  ## result
  out <- list(S1=S1_result, S0=S0_result, SE_S_diff=SE_S_diff_result,
              rmst1=RMST1_result, rmst0=RMST0_result, SE_rmst_diff=SE_RMST_diff_result,
              rho=rho_result, SE_rho=SE_rho_result)

  return(out)
}


