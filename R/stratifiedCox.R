#' Estimate stratified cox model of survival probability / rmst
#'
#' @param treatProb Estimated probability for treatment for each person if receive treatment 1
#' @param cenHaz
#'

strataCox <- function(treatment, eventObserved, time, treatProb, cenHaz, nsim, printSim){

  ## create strata
  psStrata <- quantile(treatProb, probs=c(0.2, 0.4, 0.6, 0.8))
  breaks <- c(0, psStrata, 1)
  breaks[1] <- -1
  stratumId <- as.integer(as.character(cut(treatProb, breaks=breaks, labels=1:5)))

  ## parameters
  n <- length(treatment)
  maxTime <- dim(cenHaz)[1]/n
  ID <- rep(1:n, each=maxTime)

  ## create dlong
  dlong <- transformData(dwide=data.frame(eventObserved=eventObserved, time=time, stratumId=stratumId, treat=treatment), freqTime=1, type="survival")
  rownames(dlong) <- NULL
  dlong <- dlong[which(dlong$t <= maxTime), c("Lt", "It", "t", "stratumId", "treat")]


  ## fit glm
  if(is.null(cenHaz)){

    fit <- glm(Lt ~ splines::ns(t, df=4) * treat + as.factor(stratumId)*splines::ns(t, df=4), subset = It == 1, data = dlong)

  }else{

    CenProb1List <- tapply(1 - cenHaz$Haz1, ID, cumprod, simplify = FALSE)
    CenProb0List <- tapply(1 - cenHaz$Haz0, ID, cumprod, simplify = FALSE)

    weight1 <- 1/unlist(CenProb1List, use.names = FALSE)
    weight0 <- 1/unlist(CenProb0List, use.names = FALSE)

    rm(list=c("CenProb1List", "CenProb0List", "cenHaz"))

    weight1[which(weight1 >= quantile(weight1, probs = 0.95))] <- quantile(weight1, probs = 0.95)
    weight0[which(weight0 >= quantile(weight0, probs = 0.95))] <- quantile(weight0, probs = 0.95)

    weight <- weight1 * dlong$treat + weight0 * (1 - dlong$treat)

    fit <- glm(Lt ~ splines::ns(t, df=4) * treat + as.factor(stratumId)*splines::ns(t, df=4),
               subset = It == 1, weights=weight, data = dlong)

  }

  ## get estimates
  haz1 <- predict(fit, newdata=dplyr::mutate(dlong, treat=1), type = 'response')
  haz0 <- predict(fit, newdata=dplyr::mutate(dlong, treat=0), type = 'response')

  ## S
  SurvProb1 <- unlist(tapply(1 - haz1, ID, cumprod, simplify = FALSE), use.names = FALSE)
  SurvProb0 <- unlist(tapply(1 - haz0, ID, cumprod, simplify = FALSE), use.names = FALSE)

  S1_result <- tapply(SurvProb1, dlong$t, mean)
  S0_result <- tapply(SurvProb0, dlong$t, mean)

  rm(list = c("haz1", "haz0", "SurvProb1", "SurvProb0"))

  rmst1_result <- cumsum(S1_result)
  rmst0_result <- cumsum(S0_result)

  ## get CI
  S1 <- S0 <- c()
  rmst1 <- rmst0 <- c()
  for (i in 1:nsim){

    ## simulate coef
    coef_temp <- MASS::mvrnorm(1, mu=coef(fit), Sigma=vcov(fit))
    ## get estimates
    fit$coefficients <- coef_temp
    haz1 <- predict(fit, newdata=dplyr::mutate(dlong, treat=1), type = 'response')
    haz0 <- predict(fit, newdata=dplyr::mutate(dlong, treat=0), type = 'response')
    ## S
    SurvProb1 <- unlist(tapply(1 - haz1, ID, cumprod, simplify = FALSE), use.names = FALSE)
    SurvProb0 <- unlist(tapply(1 - haz0, ID, cumprod, simplify = FALSE), use.names = FALSE)

    SProb1 <- tapply(SurvProb1, dlong$t, mean)
    SProb0 <- tapply(SurvProb0, dlong$t, mean)

    S1 <- rbind(S1, SProb1)
    S0 <- rbind(S0, SProb0)

    rmst1 <- rbind(rmst1, cumsum(SProb1))
    rmst0 <- rbind(rmst0, cumsum(SProb0))

    rm(list = c("haz1", "haz0", "SurvProb1", "SurvProb0", "SProb1", "SProb0"))

    if(printSim){print(paste("Simulation", i, "finished"))}

  }

  ## result
  SE_result_S <- apply(S1-S0, 2, sd)
  SE_result_rmst <- apply(rmst1-rmst0, 2, sd)

  out <- data.frame(S1=S1_result, S0=S0_result, SE_S=SE_result_S,
                    rmst1=rmst1_result, rmst0=rmst0_result, SE_rmst=SE_result_rmst)
  return(out)

}

