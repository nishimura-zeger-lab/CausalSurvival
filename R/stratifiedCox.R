#' Estimate stratified cox model of survival probability / rmst
#'
#' @param treatProb Estimated probability for treatment for each person if receive treatment 1
#' @param cenHaz
#'

strataCox <- function(treatment, eventObserved, time,
                      treatProb, cenHaz,
                      timeIntMidPoint, timeIntLength, breaks,
                      nsim, printSim, hazMethod){

  ## create strata
  psStrata <- quantile(treatProb, probs=c(0.2, 0.4, 0.6, 0.8))
  breaks2 <- c(0, psStrata, 1)
  breaks2[1] <- -1
  stratumId <- as.integer(as.character(cut(treatProb, breaks=breaks2, labels=1:5)))

  ## parameters
  n <- length(treatment)
  maxTime <- length(timeIntMidPoint)
  ID <- rep(1:n, each=maxTime)
  id <- rep(1:16, maxTime)

  ## create dlong
  dlong <- transformData(dwide=data.frame(eventObserved=eventObserved, time=time, stratumId=stratumId, treat=treatment), timeIntMidPoint=timeIntMidPoint, type="survival")
  rownames(dlong) <- NULL

  if(hazMethod == "twoStage"){

  ## mgcv
  fit <- mgcv::bam(Lt ~ s(t, bs="ps"), family = binomial, subset = It == 1, data = dlong, method="REML")
  offset_t <- predict(fit, newdata = data.frame(t=timeIntMidPoint))

  ## other parameters
  timeEffect <- "linear"
  evenKnot <- NULL

  }else if(hazMethod == "ns"){

    offset_t <- NULL
    timeEffect <- "ns"
    evenKnot <- FALSE

  }

  ## fit glm
  if(is.null(cenHaz)){

    interactWithTime <- as(as.matrix(model.matrix(~ as.factor(stratumId), data.frame(stratumId=stratumId))[, -1]), "sparseMatrix")
    covariates <- Matrix::summary(interactWithTime)
    colnames(covariates) <- c("rowId", "covariateId", "covariateValue")

    eps <- estimateHaz(id=1:length(treatment), treatment=treatment, eventObserved=eventObserved, time=time,
                       offset_t=offset_t, offset_X=FALSE, breaks=breaks, weight=NULL,
                       covariates=covariates, covIdHaz=NULL, crossFitNum=1, index_ls=NULL,
                       timeEffect=timeEffect, evenKnot=evenKnot, penalizeTimeTreatment=NULL,
                       interactWithTime=as.matrix(interactWithTime), hazEstimate="glm", intercept=TRUE,
                       estimate_hazard="survival", getHaz=FALSE, coef_H=NULL, sigma=NULL,
                       robust=FALSE, threshold=1e-10)

  }else{

    CenProb1List <- tapply(1 - cenHaz$Haz1, ID, cumprod, simplify = FALSE)
    CenProb0List <- tapply(1 - cenHaz$Haz0, ID, cumprod, simplify = FALSE)

    weight1 <- 1/unlist(CenProb1List, use.names = FALSE)
    weight0 <- 1/unlist(CenProb0List, use.names = FALSE)

    rm(list=c("CenProb1List", "CenProb0List"))

    weight1[which(weight1 >= quantile(weight1, probs = 0.95))] <- quantile(weight1, probs = 0.95)
    weight0[which(weight0 >= quantile(weight0, probs = 0.95))] <- quantile(weight0, probs = 0.95)

    weight <- weight1 * dlong$treat + weight0 * (1 - dlong$treat)
    weight <- weight[order(-dlong$t, dlong$time, 1-dlong$eventObserved, decreasing = TRUE)]

    interactWithTime <- as(as.matrix(model.matrix(~ as.factor(stratumId), data.frame(stratumId=stratumId))[, -1]), "sparseMatrix")
    covariates <- Matrix::summary(interactWithTime)
    colnames(covariates) <- c("rowId", "covariateId", "covariateValue")

    eps <- estimateHaz(id=1:length(treatment), treatment=treatment, eventObserved=eventObserved, time=time,
                       offset_t=offset_t, offset_X=FALSE, breaks=breaks, weight=weight,
                       covariates=covariates, covIdHaz=NULL, crossFitNum=1, index_ls=NULL,
                       timeEffect=timeEffect, evenKnot=evenKnot, penalizeTimeTreatment=NULL,
                       interactWithTime=as.matrix(interactWithTime), hazEstimate="glm", intercept=TRUE,
                       estimate_hazard="survival", getHaz=FALSE, coef_H=NULL, sigma=NULL,
                       robust=TRUE, threshold=1e-10)

  }

  ## get estimates
  if(hazMethod == "twoStage"){

  designM <- expand.grid(treat=1, stratumId2=c(0, 1), stratumId3=c(0, 1), stratumId4=c(0, 1), stratumId5=c(0, 1), t=timeIntMidPoint)
  designM$stratumId2t <- designM$stratumId2*designM$t
  designM$stratumId3t <- designM$stratumId3*designM$t
  designM$stratumId4t <- designM$stratumId4*designM$t
  designM$stratumId5t <- designM$stratumId5*designM$t
  haz1 <- plogis(rep(offset_t, each=16) + (as.matrix(designM) %*% eps$coef_fit[-1, 1])[, 1] + eps$coef_fit[1, 1])
  designM$treat <- 0
  haz0 <- plogis(rep(offset_t, each=16) + (as.matrix(designM) %*% eps$coef_fit[-1, 1])[, 1] + eps$coef_fit[1, 1])

  }else if(hazMethod == "ns"){

      indx_subset <- sapply(timeIntMidPoint, function(x) sum(time >= x), USE.NAMES = FALSE)
      nsBase <- splines::ns(timeIntMidPoint, knots=quantile(rep(timeIntMidPoint, times=indx_subset), probs=c(0.2, 0.4, 0.6, 0.8)))

      designM <- expand.grid(treat=1, stratumId2=c(0, 1), stratumId3=c(0, 1), stratumId4=c(0, 1), stratumId5=c(0, 1), t=1:length(timeIntMidPoint))
      designM <- cbind(designM[, -6], nsBase[t, ])
      designM <- cbind(designM, designM[, 2:5] * designM[, 6])
      designM <- cbind(designM, designM[, 2:5] * designM[, 7])
      designM <- cbind(designM, designM[, 2:5] * designM[, 8])
      designM <- cbind(designM, designM[, 2:5] * designM[, 9])
      designM <- cbind(designM, designM[, 2:5] * designM[, 10])
      haz1 <- plogis((as.matrix(designM) %*% eps$coef_fit[-1, 1])[, 1] + eps$coef_fit[1, 1])
      designM$treat <- 0
      haz0 <- plogis((as.matrix(designM) %*% eps$coef_fit[-1, 1])[, 1] + eps$coef_fit[1, 1])

  }


  ## S
  SurvProb1 <- unlist(tapply(1 - haz1, id, cumprod, simplify = FALSE), use.names = FALSE)
  SurvProb0 <- unlist(tapply(1 - haz0, id, cumprod, simplify = FALSE), use.names = FALSE)

  S1_result <- tapply(SurvProb1, rep(timeIntMidPoint, 16), mean)
  S0_result <- tapply(SurvProb0, rep(timeIntMidPoint, 16), mean)

  rm(list = c("haz1", "haz0", "SurvProb1", "SurvProb0"))

  rmst1_result <- cumsum(timeIntLength * S1_result)
  rmst0_result <- cumsum(timeIntLength * S0_result)

  ## get CI
  S1 <- S0 <- c()
  rmst1 <- rmst0 <- c()

  for (i in 1:nsim){

    ## simulate coef
    if(is.null(cenHaz)){
    coef_temp <- MASS::mvrnorm(1, mu=eps$coef_fit, Sigma=eps$coef_var)
    }else{
    coef_temp <- MASS::mvrnorm(1, mu=eps$coef_fit, Sigma=eps$coef_robustVar)
    }
    ## get estimates
    designM$treat <- 1
    if(hazMethod == "twoStage"){
    haz1 <- plogis(rep(offset_t, each=16) + (as.matrix(designM) %*% coef_temp[-1])[, 1] + coef_temp[1])
    }else if(hazMethod == "ns"){
    haz1 <- plogis((as.matrix(designM) %*% coef_temp[-1])[, 1] + coef_temp[1])
    }
    designM$treat <- 0
    if(hazMethod == "twoStage"){
    haz0 <- plogis(rep(offset_t, each=16) + (as.matrix(designM) %*% coef_temp[-1])[, 1] + coef_temp[1])
    }else if(hazMethod == "ns"){
    haz0 <- plogis((as.matrix(designM) %*% coef_temp[-1])[, 1] + coef_temp[1])
    }

    ## S
    SurvProb1 <- unlist(tapply(1 - haz1, id, cumprod, simplify = FALSE), use.names = FALSE)
    SurvProb0 <- unlist(tapply(1 - haz0, id, cumprod, simplify = FALSE), use.names = FALSE)

    SProb1 <- tapply(SurvProb1, rep(timeIntMidPoint, 16), mean)
    SProb0 <- tapply(SurvProb0, rep(timeIntMidPoint, 16), mean)

    S1 <- rbind(S1, SProb1)
    S0 <- rbind(S0, SProb0)

    rmst1 <- rbind(rmst1, cumsum(timeIntLength * SProb1))
    rmst0 <- rbind(rmst0, cumsum(timeIntLength * SProb0))

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

