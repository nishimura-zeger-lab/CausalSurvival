#' Estimate IPW of restricted mean survival time (RMST)
#'
#' @param dlong Long-format survival data from function transformData(dwide, freqTime), must include columns: id, t, treatment, Lt, It
#' @param cenHaz Data frame with two columns: CenHaz1, CenHaz0
#'               Estimated censoring hazards for each person at each time points if receive treatment 1 (CenHaz1) and if receive treatment 0 (CenHaz0)
#' @param treatProb Estimated probability for treatment for each person if receive treatment 1
#' @param tau Time of interest. Can be a vector (multiple time of interest)
#' @return A data frame with three columns: RMST1, RMST0, SErmstDiff

estimateIPWrmst <- function(dlong, cenHaz, treatProb, tau){

  ## container
  RMST1_result <- RMST0_result <- SErmstDiff_result <- rep(0, length=length(tau))


  ## Estimate survival hazards
  survHazFit <- lm(Lt ~ as.factor(t) * as.factor(treatment), subset = It == 1, data = dlong)
  survHaz1 <- bound01(predict(survHazFit, newdata = data.frame(t=1:max(dlong$t), treatment=rep(1, max(dlong$t))), type = 'response'))
  survHaz0 <- bound01(predict(survHazFit, newdata = data.frame(t=1:max(dlong$t), treatment=rep(0, max(dlong$t))), type = 'response'))

  rm(list=c("survHazFit"))

  survHaz1 <- rep(survHaz1, length(unique(dlong$id)))
  survHaz0 <- rep(survHaz0, length(unique(dlong$id)))


  ## calculate survival and censoring probability
  SurvProb1List <- tapply(1 - survHaz1, dlong$id, cumprod, simplify = FALSE)
  SurvProb0List <- tapply(1 - survHaz0, dlong$id, cumprod, simplify = FALSE)

  SurvProb1 <- unlist(SurvProb1List, use.names = FALSE)
  SurvProb0 <- unlist(SurvProb0List, use.names = FALSE)

  rm(list=c("SurvProb1List", "SurvProb0List"))

  CenProb1List <- tapply(1 - cenHaz$Haz1, dlong$id, cumprod, simplify = FALSE)
  CenProb0List <- tapply(1 - cenHaz$Haz0, dlong$id, cumprod, simplify = FALSE)

  CenProb1 <- unlist(CenProb1List, use.names = FALSE)
  CenProb0 <- unlist(CenProb0List, use.names = FALSE)

  rm(list=c("CenProb1List", "CenProb0List", "cenHaz"))

  ## Observed estimated survival hazards
  SurvHaz_obs <- dlong$treatment * survHaz1 + (1-dlong$treatment) * survHaz0

  rm(list=c("survHaz1", "survHaz0"))

  for (TimePoint in tau){

    if(TimePoint == 1){next}

    ## parameter
    ind <- (dlong$t <= TimePoint-1)

    ## IPW
    indTillTimePoint <- unlist(tapply(ind, dlong$id, function(x){rev(cumsum(rev(x)))}), use.names = FALSE)
    Z1ipw <- - indTillTimePoint / bound(treatProb[dlong$id] * CenProb1)
    Z0ipw <- - indTillTimePoint / bound((1-treatProb[dlong$id]) * CenProb0)

    rm(list=c("indTillTimePoint"))

    DT1ipw <- with(dlong, tapply(It * treatment * Z1ipw * Lt, id, sum))
    DT0ipw <- with(dlong, tapply(It * (1-treatment) * Z0ipw * Lt , id, sum))

    ipw  <- tau + c(mean(DT0ipw), mean(DT1ipw))

    rm(list=c("Z1ipw", "Z0ipw", "DT1ipw", "DT0ipw"))

    ## variance
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
    ## standard error
    D <- DT1 - DT0 + DW1 - DW0
    sdn <- sqrt(var(D) / length(unique(dlong$id)))

    rm(list=c("D", "DW1", "DW0", "DT1", "DT0", "ind"))

    ## store
    RMST1_result[TimePoint] <- ipw[2]
    RMST0_result[TimePoint] <- ipw[1]
    SErmstDiff_result[TimePoint] <- sdn

  }

  ## result
  out <- data.frame(RMST1=RMST1_result, RMST0=RMST0_result, SErmstDiff=SErmstDiff_result)
  return(out)

}



