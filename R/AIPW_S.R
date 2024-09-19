#' Estimate cross-fitted Augmented IPW of survival probability at time tau
#'
#' @param dlong Long-format survival data from function transformData(dwide, freqTime), must include columns: id, t, treatment, Lt, It
#' @param survHaz Data frame with two columns: SurvHaz1, SurvHaz0
#'                Estimated survival hazards for each person at each time points if receive treatment 1 (SurvHaz1) and if receive treatment 0 (SurvHaz0)
#' @param cenHaz Data frame with two columns: CenHaz1, CenHaz0
#'               Estimated censoring hazards for each person at each time points if receive treatment 1 (CenHaz1) and if receive treatment 0 (CenHaz0)
#' @param treatProb Estimated probability for treatment for each person if receive treatment 1
#' @param tau Time of interest. Can be a vector (multiple time of interest)
#' @return A data frame with three columns: SurvProb1, SurvProb0, SEprobDiff

estimateAIPWprob <- function(dlong, survHaz, cenHaz, treatProb, tau){


  ## container
  SurvProb1_result <- SurvProb0_result <- SEprobDiff_result <- rep(0, length=length(tau))

  ## calculate survival and censoring probability
  SurvProb1List <- tapply(1 - survHaz$SurvHaz1, dlong$id, cumprod, simplify = FALSE)
  SurvProb0List <- tapply(1 - survHaz$SurvHaz0, dlong$id, cumprod, simplify = FALSE)

  SurvProb1 <- unlist(SurvProb1List, use.names = FALSE)
  SurvProb0 <- unlist(SurvProb0List, use.names = FALSE)

  rm(list=c("SurvProb1List", "SurvProb0List"))

  CenProb1List <- tapply(1 - cenHaz$CenHaz1, dlong$id, cumprod, simplify = FALSE)
  CenProb0List <- tapply(1 - cenHaz$CenHaz0, dlong$id, cumprod, simplify = FALSE)

  CenProb1 <- unlist(CenProb1List, use.names = FALSE)
  CenProb0 <- unlist(CenProb0List, use.names = FALSE)

  rm(list=c("CenProb1List", "CenProb0List", "cenHaz"))

  ## Observed estimated survival hazards
  SurvHaz_obs <- dlong$treatment * survHaz$SurvHaz1 + (1-dlong$treatment) * survHaz$SurvHaz0

  rm(list=c("survHaz"))


  for (TimePoint in tau){

    ## parameter
    ind <- (dlong$t <= TimePoint)

    ## solve estimating equation
    H1 <- - (ind * rep(SurvProb1[which(dlong$t == TimePoint)], each=max(dlong$t))) / bound(SurvProb1 * treatProb[dlong$id] * CenProb1)
    H0 <- - (ind * rep(SurvProb0[which(dlong$t == TimePoint)], each=max(dlong$t))) / bound(SurvProb0 * (1-treatProb[dlong$id]) * CenProb0)
    DT1 <- with(dlong, tapply(It * treatment * H1 * (Lt - SurvHaz_obs), id, sum))
    DT0 <- with(dlong, tapply(It * (1 - treatment) * H0 * (Lt - SurvHaz_obs), id, sum))

    rm(list=c("ind", "H1", "H0"))

    DW1 <- SurvProb1[which(dlong$t == TimePoint)]
    DW0 <- SurvProb0[which(dlong$t == TimePoint)]

    ## AIPW
    aipw <- c(mean(DT0 + DW0), mean(DT1 + DW1))
    ## SE
    D <- DT1 - DT0 + DW1 - DW0
    sdn <- sqrt(var(D) / length(unique(dlong$id)))

    rm(list=c("D", "DW1", "DW0", "DT1", "DT0"))

    ## store
    SurvProb1_result[TimePoint] <- aipw[2]
    SurvProb0_result[TimePoint] <- aipw[1]
    SEprobDiff_result[TimePoint] <- sdn

  }

  ## result
  out <- data.frame(SurvProb1=SurvProb1_result, SurvProb0=SurvProb0_result, SEprobDiff=SEprobDiff_result)
  return(out)
}





