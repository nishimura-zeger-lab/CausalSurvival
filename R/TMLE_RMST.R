#' Estimate cross-fitted TMLE of restricted mean survival time (RMST) at time tau
#'
#' @param dlong Long-format survival data from function transformData(dwide, freqTime), must include columns: id, t, treatment, Lt, It
#' @param survHaz Data frame with two columns: SurvHaz1, SurvHaz0
#'                Estimated survival hazards for each person at each time points if receive treatment 1 (SurvHaz1) and if receive treatment 0 (SurvHaz0)
#' @param cenHaz Data frame with two columns: CenHaz1, CenHaz0
#'               Estimated censoring hazards for each person at each time points if receive treatment 1 (CenHaz1) and if receive treatment 0 (CenHaz0)
#' @param treatProb Estimated probability for treatment for each person if receive treatment 1
#' @param tau Time of interest. Can be a vector (multiple time of interest)
#' @return A data frame with three columns: RMST1, RMST0, SErmstDiff

estimateTMLErmst <- function(dlong, survHaz, cenHaz, treatProb, tau){

  ## container
  RMST1_result <- RMST0_result <- SErmstDiff_result <- rep(0, length=length(tau))


  ## calculate censoring probability that doesn't need iterative updates
  CenProb1List <- tapply(1 - cenHaz$CenHaz1, dlong$id, cumprod, simplify = FALSE)
  CenProb0List <- tapply(1 - cenHaz$CenHaz0, dlong$id, cumprod, simplify = FALSE)

  CenProb1 <- unlist(CenProb1List, use.names = FALSE)
  CenProb0 <- unlist(CenProb0List, use.names = FALSE)

  rm(list=c("CenProb1List", "CenProb0List", "cenHaz"))


  for (TimePoint in tau){

    if(TimePoint == 1){next}

    ## initial SurvHaz
    SurvHaz1 <- survHaz$SurvHaz1
    SurvHaz0 <- survHaz$SurvHaz0
    SurvHaz_obs  <- dlong$treatment*SurvHaz1 + (1-dlong$treatment)*SurvHaz0

    ## parameter
    ind <- (dlong$t <= TimePoint-1)
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
      cumProb1TillTimePoint <- unlist(tapply(ind * SurvProb1, dlong$id, function(x){rev(cumsum(rev(x)))}), use.names = FALSE)
      H1 <- - cumProb1TillTimePoint / bound(SurvProb1 * treatProb[dlong$id] * CenProb1)

      cumProb0TillTimePoint <- unlist(tapply(ind * SurvProb0, dlong$id, function(x){rev(cumsum(rev(x)))}), use.names = FALSE)
      H0 <- - cumProb0TillTimePoint / bound(SurvProb0 * (1-treatProb[dlong$id]) * CenProb0)

      rm(list=c("cumProb1TillTimePoint", "cumProb0TillTimePoint"))

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

    cumProb1TillTimePoint <- unlist(tapply(ind * SurvProb1, dlong$id, function(x){rev(cumsum(rev(x)))}), use.names = FALSE)
    H1 <- - cumProb1TillTimePoint / bound(SurvProb1 * treatProb[dlong$id] * CenProb1)

    cumProb0TillTimePoint <- unlist(tapply(ind * SurvProb0, dlong$id, function(x){rev(cumsum(rev(x)))}), use.names = FALSE)
    H0 <- - cumProb0TillTimePoint / bound(SurvProb0 * (1-treatProb[dlong$id]) * CenProb0)

    rm(list=c("cumProb1TillTimePoint", "cumProb0TillTimePoint"))

    DT <- with(dlong, tapply(It * (treatment * H1 - (1 - treatment) * H0) * (Lt - SurvHaz_obs), id, sum))

    DW1 <- tapply(ind * SurvProb1, dlong$id, sum)
    DW0 <- tapply(ind * SurvProb0, dlong$id, sum)

    rm(list=c("SurvProb1", "SurvProb0"))

    ## RMST(1, tau)
    RMST1_mean <- mean(DW1)
    ## RMST(0, tau)
    RMST0_mean <- mean(DW0)
    ## standard error of RMST(1, tau)-RMST(0, tau)
    D <- DT + DW1 - DW0
    sdn <- sqrt(var(D) / length(unique(dlong$id)))

    rm(list=c("DT", "DW1", "DW0", "D"))

    ## store
    RMST1_result[TimePoint] <- RMST1_mean
    RMST0_result[TimePoint] <- RMST0_mean
    SErmstDiff_result[TimePoint] <- sdn

  }

  ## result
  out <- data.frame(RMST1=SurvRMST1, SurvRMST0=SurvRMST0, SErmstDiff=SErmstDiff_result)
  return(out)
}





