#' Estimate Augmented IPW of survival probability at time tau, RMST at time tau, log-ratio of log-survival probability
#'
#' @param treatment
#' @param eventObserved
#' @param time
#' @param survHaz Data frame with two columns: Haz1, Haz0
#'                Estimated survival hazards for each person at each time points if receive treatment 1 (Haz1) and if receive treatment 0 (Haz0)
#' @param cenHaz Data frame with two columns: Haz1, Haz0
#'               Estimated censoring hazards for each person at each time points if receive treatment 1 (Haz1) and if receive treatment 0 (Haz0)
#' @param treatProb Estimated probability for treatment for each person if receive treatment 1
#' @param tau Time of interest. Can be a vector (multiple time of interest)
#' @param timeIntMidPoint
#' @param timeIntLength
#' @param printTau
#' @return A list with estimates and square root variance of survival probability at time tau, RMST at time tau, log-ratio of log-survival probability

estimateAIPW <- function(treatment, eventObserved, time,
                         survHaz, cenHaz, treatProb,
                         tau, timeIntMidPoint, timeIntLength,
                         printTau){


  ## container
  S1_result <- S0_result <- SE_S_diff_result <- rep(0, length=length(tau))
  RMST1_result <- RMST0_result <- SE_RMST_diff_result <- rep(0, length=length(tau))
  rho_result <- SE_rho_result <- 0
  D_rmst <- 0
  D1_rmst <- D0_rmst <- 0
  D_rho <- W_rho <- rho_2term <- 0

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
  SurvHaz_obs  <- treatment[survHaz$ID] * survHaz$Haz1  + (1 - treatment[survHaz$ID]) * survHaz$Haz0
  rm(list=c("survHaz"))

  ## calculate censoring probability
  prodlag <- function(x) cumprod(c(1, x[-length(x)]))
  CenProb1List <- tapply(1 - cenHaz$Haz1, ID, prodlag, simplify = FALSE)
  CenProb0List <- tapply(1 - cenHaz$Haz0, ID, prodlag, simplify = FALSE)

  CenProb1 <- unlist(CenProb1List, use.names = FALSE)
  CenProb0 <- unlist(CenProb0List, use.names = FALSE)

  rm(list=c("CenProb1List", "CenProb0List", "cenHaz"))

  ## dlong
  dlong <- transformData(time=time, eventObserved=eventObserved,
                         timeIntMidPoint=timeIntMidPoint, type="survival")
  rownames(dlong) <- NULL

  ## denominator
  weightH1 <- 1/(SurvProb1 * treatProb[ID] * CenProb1)
  weightH0 <- 1/(SurvProb0 * (1-treatProb[ID]) * CenProb0)

  weightH1[which(weightH1 >= quantile(weightH1, probs = 0.95))] <- quantile(weightH1, probs = 0.95)
  weightH0[which(weightH0 >= quantile(weightH0, probs = 0.95))] <- quantile(weightH0, probs = 0.95)


  for (TimePoint in 1:length(tau)){

    ## parameter
    ind <- (dlong$time <= tau[TimePoint])

    ## solve estimating equation
    H1 <- as(matrix( - (ind * rep(SurvProb1[which(dlong$time == tau[TimePoint])], each=maxTime)) * weightH1, ncol = 1), "sparseMatrix")
    H0 <- as(matrix( - (ind * rep(SurvProb0[which(dlong$time == tau[TimePoint])], each=maxTime)) * weightH0, ncol = 1), "sparseMatrix")

    DT1 <- tapply(dlong$valid * treatment[ID] * H1[, 1] * (dlong$y - SurvHaz_obs), ID, sum)
    DT0 <- tapply(dlong$valid * (1 - treatment[ID]) * H0[, 1] * (dlong$y - SurvHaz_obs), ID, sum)

    rm(list=c("H1", "H0"))

    DW1 <- SurvProb1[which(dlong$time == tau[TimePoint])]
    DW0 <- SurvProb0[which(dlong$time == tau[TimePoint])]

    aipw_S <- c(mean(DT0 + DW0), mean(DT1 + DW1))
    D1_rmst <- D1_rmst + (DT1 + DW1)* timeIntLength[TimePoint]
    D0_rmst <- D0_rmst + (DT0 + DW0)* timeIntLength[TimePoint]

    ## SE
    D <- DT1 - DT0 + DW1 - DW0
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

    ## store
    S1_result[TimePoint] <- aipw_S[2]
    S0_result[TimePoint] <- aipw_S[1]
    SE_S_diff_result[TimePoint] <- sdn_S

    RMST1_result[TimePoint] <- mean(D1_rmst)
    RMST0_result[TimePoint] <- mean(D0_rmst)
    SE_RMST_diff_result[TimePoint] <- sdn_rmst

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





