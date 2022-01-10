#' Estimate cross-fitted Augmented IPW of restricted mean survival time (RMST)
#'
#' @param dlong Long-format survival data from function transformData(dwide, freqTime), must include columns: id, t, treatment, Lt, It
#' @param survHaz Data frame with two columns: SurvHaz1, SurvHaz0
#'                Estimated survival hazards for each person at each time points if receive treatment 1 (SurvHaz1) and if receive treatment 0 (SurvHaz0)
#' @param cenHaz Data frame with two columns: CenHaz1, CenHaz0
#'               Estimated censoring hazards for each person at each time points if receive treatment 1 (CenHaz1) and if receive treatment 0 (CenHaz0)
#' @param treatProb Estimated probability for treatment for each person if receive treatment 1
#' @param tau Time of interest. Can be a vector (multiple time of interest)
#' @return A data frame with three columns: RMST1, RMST0, SErmstDiff

estimateAIPWrmst <- function(dlong, survHaz, cenHaz, treatProb, tau){

  ## container
  RMST1_result <- RMST0_result <- SErmstDiff_result <- rep(0, length=length(tau))

  ## calculate survival and censoring probability
  SurvProb1List <- tapply(1 - survHaz$Haz1, dlong$id, cumprod, simplify = FALSE)
  SurvProb0List <- tapply(1 - survHaz$Haz0, dlong$id, cumprod, simplify = FALSE)

  SurvProb1 <- unlist(SurvProb1List, use.names = FALSE)
  SurvProb0 <- unlist(SurvProb0List, use.names = FALSE)

  rm(list=c("SurvProb1List", "SurvProb0List"))

  CenProb1List <- tapply(1 - cenHaz$Haz1, dlong$id, cumprod, simplify = FALSE)
  CenProb0List <- tapply(1 - cenHaz$Haz0, dlong$id, cumprod, simplify = FALSE)

  CenProb1 <- unlist(CenProb1List, use.names = FALSE)
  CenProb0 <- unlist(CenProb0List, use.names = FALSE)

  rm(list=c("CenProb1List", "CenProb0List", "cenHaz"))

  ## Observed estimated survival hazards
  SurvHaz_obs <- dlong$treatment * survHaz$SurvHaz1 + (1-dlong$treatment) * survHaz$SurvHaz0

  rm(list=c("survHaz"))


  for (TimePoint in tau){

    if(TimePoint == 1){next}

    ## parameter
    ind <- (dlong$t <= TimePoint-1)

    ## solve estimating equation
    cumProb1TillTimePoint <- unlist(tapply(ind * SurvProb1, dlong$id, function(x){rev(cumsum(rev(x)))}), use.names = FALSE)
    H1 <- - cumProb1TillTimePoint / bound(SurvProb1 * treatProb[dlong$id] * CenProb1)
    cumProb0TillTimePoint <- unlist(tapply(ind * SurvProb0, dlong$id, function(x){rev(cumsum(rev(x)))}), use.names = FALSE)
    H0 <- - cumProb0TillTimePoint / bound(SurvProb0 * (1-treatProb[dlong$id]) * CenProb0)

    rm(list=c("cumProb1TillTimePoint", "cumProb0TillTimePoint"))

    DT1 <- with(dlong, tapply(It * treatment * H1 * (Lt - SurvHaz_obs), id, sum))
    DT0 <- with(dlong, tapply(It * (1 - treatment) * H0 * (Lt - SurvHaz_obs), id, sum))

    rm(list=c("H1", "H0"))

    DW1 <- tapply(ind * SurvProb1, dlong$id, sum)
    DW0 <- tapply(ind * SurvProb0, dlong$id, sum)
    ## AIPW
    aipw <- c(mean(DT0 + DW0), mean(DT1 + DW1))
    ## SE
    D <- DT1 - DT0 + DW1 - DW0
    sdn <- sqrt(var(D) / length(unique(dlong$id)))

    rm(list=c("D", "DW1", "DW0", "DT1", "DT0", "ind"))

    ## store
    RMST1_result[TimePoint] <- aipw[2]
    RMST0_result[TimePoint] <- aipw[1]
    SErmstDiff_result[TimePoint] <- sdn

  }

  ## result
  out <- data.frame(RMST1=RMST1_result, RMST0=RMST0_result, SErmstDiff=SErmstDiff_result)
  return(out)
}























