#' Get initial parameters for simulation
#' @param nInt number of time intervals for coarsening the data
#' @param hazEstimate "survival" or "censoring"
#' @param hazMethod "twoStage" or "ns"
#' @param seed to set.seed()

estimateSimulationParams <- function(outcome, time, treatment, covariates,
                                     simOutcome=NULL, simTime=NULL, covId=NULL,
                                     nInt=NULL, hazEstimate, hazMethod, seed=0){

  ## parameters
  cov <- Matrix::sparseMatrix(i = covariates$rowId, j = covariates$covariateId, x = covariates$covariateValue, repr = "T")
  rowId <- 1:length(outcome)
  if(is.null(nInt)){nInt <- min(50, floor((sum(outcome)/10)))}
  cData <- coarsenData(time=time, outcome=outcome, nInt=nInt)

  if(is.null(covId)){
  ## variable selection from real data
  if(hazEstimate == "censoring"){
    d_outcome <- 1 - outcome
    sigma <- exp(seq(log(0.5), log(0.01), length.out = 20))
  }else if(hazEstimate == "survival"){
    d_outcome <- outcome
    sigma <- exp(seq(log(1), log(0.01), length.out = 20))
  }
  set.seed(seed)
  fit <- glmnet::cv.glmnet(x=cbind(treatment, cov), y=Surv(time=time, event=d_outcome), family = "cox", nfolds = 5, penalty.factor = c(0, rep(1, dim(cov)[2])))
  cf <- coef(fit, s = fit$lambda.1se)
  cov_indx <- setdiff(which(cf != 0)-1, 0)
  rm(list=c("fit", "cf"))
  }else{
  ## directly choose covariates
    if(hazEstimate == "censoring"){
      sigma <- exp(seq(log(0.5), log(0.01), length.out = 20))
    }else if(hazEstimate == "survival"){
      sigma <- exp(seq(log(1), log(0.01), length.out = 20))
    }
    cov <- cov[, covId]
    cov_indx <- covId
  }

  ## hazards from real data
  if(!is.null(simTime)){
    outcome <- simOutcome
    cData$timeInt <- simTime
  }

  if(hazMethod == "twoStage"){

    ## long-format data
    dlong <- transformData(dwide=data.frame(eventObserved=outcome, time=cData$timeInt), timeIntMidPoint=cData$timeIntMidPoint, type="survival")
    rownames(dlong) <- NULL
    dlong <- dlong[, c("Lt", "It", "t")]

    ## offset
    fit <- mgcv::bam(Lt ~ s(t, bs="ps"), family = binomial, subset = It == 1, data = dlong, method="REML")
    offset_t <- predict(fit, newdata = data.frame(t=cData$timeIntMidPoint))
    rm(list=c("dlong", "fit"))

    timeEffect <- "linear"
    evenKnot <- NULL

  }else if(hazMethod == "ns"){

    offset_t <- NULL
    timeEffect <- "ns"
    evenKnot <- FALSE

  }

  haz <- estimateHaz(id=rowId, treatment=treatment, eventObserved=outcome, time=cData$timeInt,
                     offset_t=offset_t, breaks=cData$breaks, covariates=covariates, covIdHaz=cov_indx,
                     timeEffect=timeEffect, evenKnot=evenKnot, penalizeTimeTreatment=FALSE,
                     interactWithTime=treatment, hazEstimate="ridge",
                     sigma=sigma, estimate_hazard=hazEstimate, getHaz=TRUE, coef_H=NULL)

  return(list(haz=haz, cov_indx=cov_indx, coarsenedData=cData))

}



#' Simulate time-to-event data with survival and censoring hazards
#' @param survHaz Output from estimateHaz
#' @param cenHaz Output from estimateHaz
#' @param nInt number of time intervals for coarsening the data
#'

simData <- function(treatment, survHaz, cenHaz, coarsenedTime, seed){

  n <- length(treatment)
  nInt <- length(coarsenedTime$timeIntMidPoint)

  survHaz_all <- survHaz$Haz1 * treatment[survHaz$ID] + survHaz$Haz0 * (1 - treatment[survHaz$ID])
  cenHaz_all <- cenHaz$Haz1 * treatment[cenHaz$ID] + cenHaz$Haz0 * (1 - treatment[cenHaz$ID])

  set.seed(seed)
  rS <- rbinom(length(survHaz_all), 1, survHaz_all)
  rG <- rbinom(length(cenHaz_all), 1, cenHaz_all)

  tS <- tapply(rS, rep(1:n, each=nInt), function(x){which(x == 1)[1]})
  tG <- tapply(rG, rep(1:n, each=nInt), function(x){which(x == 1)[1]})
  tG <- ifelse(is.na(tG), nInt, tG)

  ObservedTime <- apply(cbind(tS, tG), 1, function(x){min(x, na.rm = TRUE)})
  ObservedEvent <- 1*(tS<=tG)
  ObservedEvent <- ifelse(is.na(ObservedEvent), 0, ObservedEvent)

  ObservedEvent <- as.double(ObservedEvent)
  ObservedTime <- as.double(as.character(factor(ObservedTime, labels = cData$timeIntMidPoint)))

  return(list(ObservedEvent=ObservedEvent,
              ObservedTime=ObservedTime))

}



#' Calculate counterfactuals from simulated data
#' @param coarsenedTime Output from coarsenData
#' @param survHaz Output from estimateHaz
#' @param survCurve Output from calculateSurvCurve
#'

calculateSurvCurve <- function(coarsenedTime, survHaz){

  n <- length(coarsenedTime$timeInt)
  nInt <- length(coarsenedTime$timeIntMidPoint)

  Sm0 <- unlist(tapply(1-survHaz$Haz0, rep(1:n, each=nInt), cumprod, simplify = FALSE), use.names = FALSE)
  S0 <- tapply(Sm0, rep(1:nInt, n), mean)
  rm(list=c("Sm0"))
  Sm1 <- unlist(tapply(1-survHaz$Haz1, rep(1:n, each=nInt), cumprod, simplify = FALSE), use.names = FALSE)
  S1 <- tapply(Sm1, rep(1:nInt, n), mean)
  rm(list=c("Sm1"))

  return(data.frame(S0=S0, S1=S1))

}

calculateRMST <- function(coarsenedTime, survCurve){

  rmst0 <- cumsum(coarsenedTime$timeIntLength * survCurve$S0)
  rmst1 <- cumsum(coarsenedTime$timeIntLength * survCurve$S1)

  return(data.frame(RMST0=rmst0, RMST1=rmst1))

}



#' Simulate censoring time, nonproportional
#'
#
simCenTime <- function(treatment, covariates, lambda1, nu1, lambda0, nu0, seed){

  n <- length(treatment)
  nvar <- length(unique(covariates$covariateId))
  cov <- Matrix::sparseMatrix(i = covariates$rowId, j = covariates$covariateId, x = covariates$covariateValue, repr = "T")
  time <- rep(0, n)

  set.seed(seed)
  u <- runif(n)
  beta <- log(runif(nvar, min=1, max=5))

  time[treatment == 0] <- -(u[treatment == 0]/(lambda0 * exp(cov %*% beta)[, 1]))^(1/nu0)
  time[treatment == 1] <- -(u[treatment == 1]/(lambda1 * exp(log(5) + cov %*% beta)[, 1]))^(1/nu1)

  return(time)

}

























