#' Estimate cross-fitted TMLE of survival probability at time tau
#'
#' @param dlong Long-format survival data from function transformData(dwide, freqTime), must include columns: id, t, treatment, Lt, It
#' @param survHaz Data frame with two columns: SurvHaz1, SurvHaz0
#'                Estimated survival hazards for each person at each time points if receive treatment 1 (SurvHaz1) and if receive treatment 0 (SurvHaz0)
#' @param cenHaz Data frame with two columns: CenHaz1, CenHaz0
#'               Estimated censoring hazards for each person at each time points if receive treatment 1 (CenHaz1) and if receive treatment 0 (CenHaz0)
#' @param treatProb Estimated probability for treatment for each person if receive treatment 1
#' @param tau Time of interest. Can be a vector (multiple time of interest)
#' @return A data frame with three columns: SurvProb1, SurvProb0, SEprobDiff

estimateTMLEprob <- function(dlong, survHaz, cenHaz, treatProb, tau){

  ## container
  SurvProb1_result <- SurvProb0_result <- SEprobDiff_result <- rep(0, length=length(tau))

  ## calculate censoring probability that doesn't need iterative updates
  CenProb1List <- tapply(1 - cenHaz$CenHaz1, dlong$id, cumprod, simplify = FALSE)
  CenProb0List <- tapply(1 - cenHaz$CenHaz0, dlong$id, cumprod, simplify = FALSE)

  CenProb1 <- unlist(CenProb1List, use.names = FALSE)
  CenProb0 <- unlist(CenProb0List, use.names = FALSE)

  rm(list=c("CenProb1List", "CenProb0List", "cenHaz"))

  for (TimePoint in tau){

    ## initial SurvHaz
    SurvHaz1 <- survHaz$SurvHaz1
    SurvHaz0 <- survHaz$SurvHaz0
    SurvHaz_obs  <- dlong$treatment*SurvHaz1 + (1-dlong$treatment)*SurvHaz0

    ## parameter
    ind <- (dlong$t <= TimePoint)
    converged <- FALSE
    iter <- 1

    ## iterate
    while((!converged) && iter <= 20){

      SurvProb1List <- tapply(1 - SurvHaz1, dlong$id, cumprod, simplify = FALSE)
      SurvProb0List <- tapply(1 - SurvHaz0, dlong$id, cumprod, simplify = FALSE)

      SurvProb1 <- unlist(SurvProb1List, use.names = FALSE)
      SurvProb0 <- unlist(SurvProb0List, use.names = FALSE)

      rm(list=c("SurvProb1List", "SurvProb0List"))

      ## clever covariate for updating survival hazard
      H1 <- - (ind * rep(SurvProb1[which(dlong$t == TimePoint)], each=max(dlong$t))) / bound(SurvProb1 * treatProb[dlong$id] * CenProb1)
      H0 <- - (ind * rep(SurvProb0[which(dlong$t == TimePoint)], each=max(dlong$t))) / bound(SurvProb0 * (1-treatProb[dlong$id]) * CenProb0)
      H <- dlong$treatment * H1 + (1-dlong$treatment) * H0

      ## update for survival hazard
      eps   <- coef(glm2::glm2(Lt ~ 0 + offset(qlogis(SurvHaz_obs)) + H,
                         family = binomial(), subset = It == 1, data = dlong))

      ## NA as 0 for the new values
      eps[is.na(eps)] <- 0

      ## update values
      SurvHaz1 <- bound01(plogis(qlogis(SurvHaz1) + eps * H1))
      SurvHaz0 <- bound01(plogis(qlogis(SurvHaz0) + eps * H0))
      SurvHaz_obs  <- dlong$treatment * SurvHaz1  + (1 - dlong$treatment) * SurvHaz0

      iter <-  iter + 1
      converged <- (abs(eps) <= 1e-3/length(unique(dlong$id))^(0.6))

      ## clear
      rm(list=c("H1", "H0", "H","SurvProb1", "SurvProb0"))

    }

    ## result for time tau
    SurvProb1List <- tapply(1 - SurvHaz1, dlong$id, cumprod, simplify = FALSE)
    SurvProb0List <- tapply(1 - SurvHaz0, dlong$id, cumprod, simplify = FALSE)

    SurvProb1 <- unlist(SurvProb1List, use.names = FALSE)
    SurvProb0 <- unlist(SurvProb0List, use.names = FALSE)

    rm(list=c("SurvProb1List", "SurvProb0List", "SurvHaz1", "SurvHaz0"))

    H1 <- - (ind * rep(SurvProb1[which(dlong$t == TimePoint)], each=max(dlong$t))) / bound(SurvProb1 * treatProb[dlong$id] * CenProb1)
    H0 <- - (ind * rep(SurvProb0[which(dlong$t == TimePoint)], each=max(dlong$t))) / bound(SurvProb0 * (1-treatProb[dlong$id]) * CenProb0)
    DT <- with(dlong, tapply(It * (treatment * H1 - (1 - treatment) * H0) * (Lt - SurvHaz_obs), id, sum))

    rm(list=c("H1", "H0"))

    DW1 <- SurvProb1[which(dlong$t == TimePoint)]
    DW0 <- SurvProb0[which(dlong$t == TimePoint)]

    rm(list=c("SurvProb1", "SurvProb0"))

    ## S(1, tau)
    SurvProb1_mean <- mean(DW1)
    ## S(0, tau)
    SurvProb0_mean <- mean(DW0)
    ## standard error of S(1, tau)-S(0, tau)
    D <- DT + DW1 - DW0
    sdn <- sqrt(var(D) / length(unique(dlong$id)))

    rm(list=c("DT", "DW1", "DW0", "D"))

    ## store
    SurvProb1_result[TimePoint] <- SurvProb1_mean
    SurvProb0_result[TimePoint] <- SurvProb0_mean
    SEprobDiff_result[TimePoint] <- sdn

}

  ## result
  out <- data.frame(SurvProb1=SurvProb1_result, SurvProb0=SurvProb0_result, SEprobDiff=SEprobDiff_result)
  return(out)
}





