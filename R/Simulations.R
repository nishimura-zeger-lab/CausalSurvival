#' Get initial parameters for simulation
#' @param nInt number of time intervals for coarsening the data
#' @param hazEstimate "survival" or "censoring"
#' @param hazMethod "twoStage" or "ns"
#' @param seed to set.seed()

estimateSimulationParams <- function(outcome, time, treatment, covariates,
                                     simOutcome, simTime, covId,
                                     nInt, hazEstimate, hazMethod, seed){

  ## parameters ##
  ## cov
  cov <- Matrix::sparseMatrix(i = covariates$i, j = covariates$j, x = covariates$val, repr = "T")
  ## id
  rowId <- 1:length(outcome)

  if(is.null(covId)){


  ## variable selection from real data ##
  if(hazEstimate == "censoring"){
    d_outcome <- 1 - outcome
    sigma <- exp(seq(log(0.5), log(0.01), length.out = 20))
  }else if(hazEstimate == "survival"){
    d_outcome <- outcome
    sigma <- exp(seq(log(1), log(0.01), length.out = 20))
  }

  ## variable selection
  set.seed(seed)
  fit <- glmnet::cv.glmnet(x=cbind(treatment, cov), y=Surv(time=time, event=d_outcome), family = "cox", nfolds = 5, penalty.factor = c(0, rep(1, dim(cov)[2])))
  cf <- coef(fit, s = fit$lambda.1se)
  cov_indx <- setdiff(which(cf != 0)-1, 0)
  rm(list=c("fit", "cf"))

  }else{

    if(hazEstimate == "censoring"){
      sigma <- exp(seq(log(0.5), log(0.01), length.out = 20))
    }else if(hazEstimate == "survival"){
      sigma <- exp(seq(log(1), log(0.01), length.out = 20))
    }
    cov <- cov[, covId]
    cov_indx <- covId

  }


  ## hazards from real data ##
  ## coarsen data
  cData <- coarseData(time=time, outcome=outcome, nInt=nInt)
  if(!is.null(simTime)){
    outcome <- simOutcome
    cData$timeInt <- simTime
  }

  if(hazMethod == "twoStage"){

    # dlong
    dlong <- transformData(dwide=data.frame(eventObserved=outcome, time=cData$timeInt), timeIntMidPoint=cData$timeIntMidPoint, type="survival")
    rownames(dlong) <- NULL
    dlong <- dlong[, c("Lt", "It", "t")]

    ## offset
    fit <- mgcv::bam(Lt ~ s(t, bs="ps"), family = binomial, subset = It == 1, data = dlong, method="REML")
    offset_t <- predict(fit, newdata = data.frame(t=cData$timeIntMidPoint))
    rm(list=c("dlong", "fit"))

    ## other parameters
    timeEffect <- "linear"
    evenKnot <- NULL

  }else if(hazMethod == "ns"){

    offset_t <- NULL
    timeEffect <- "ns"
    evenKnot <- FALSE

  }

  haz <- estimateHaz(id=rowId, treatment=treatment, eventObserved=outcome, time=cData$timeInt,
                     offset_t=offset_t, offset_X=FALSE, intercept=TRUE, breaks=cData$breaks,
                     covariates=covariates, covIdHaz=cov_indx, crossFitNum=1, index_ls=NULL,
                     timeEffect=timeEffect, evenKnot=evenKnot, penalizeTimeTreatment=FALSE,
                     interactWithTime=treatment, hazEstimate="ridge", weight=NULL,
                     sigma=sigma, estimate_hazard=hazEstimate, getHaz=TRUE, coef_H=NULL,
                     robust=FALSE, threshold=1e-10)

  ## output ##
  return(list(haz=haz, cov_indx=cov_indx))

}



#' Simulate time-to-event data with survival and censoring hazards
#' @param survHaz Output from estimateHaz
#' @param cenHaz Output from estimateHaz
#' @param nInt number of time intervals for coarsening the data
#'

simData <- function(time, outcome, treatment, survHaz, cenHaz, nInt){

  ## parameters
  n <- length(treatment)

  ## coarsen data
  cData <- coarseData(time=time, outcome=outcome, nInt=nInt)

  ## simulation
  survHaz_all <- survHaz$Haz1 * treatment[survHaz$ID] + survHaz$Haz0 * (1 - treatment[survHaz$ID])
  cenHaz_all <- cenHaz$Haz1 * treatment[cenHaz$ID] + cenHaz$Haz0 * (1 - treatment[cenHaz$ID])

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
#' @param survHaz Output from estimateHaz
#' @param nInt number of time intervals for coarsening the data

counterFactuals <- function(time, outcome, survHaz, nInt){

  ## paramters
  n <- dim(survHaz)[1]/nInt

  ## coarsen data
  cData <- coarseData(time=time, outcome=outcome, nInt=nInt)

  ## surv prob
  Sm0 <- unlist(tapply(1-survHaz$Haz0, rep(1:n, each=nInt), cumprod, simplify = FALSE), use.names = FALSE)
  S0 <- tapply(Sm0, rep(1:nInt, n), mean)
  rm(list=c("Sm0"))
  Sm1 <- unlist(tapply(1-survHaz$Haz1, rep(1:n, each=nInt), cumprod, simplify = FALSE), use.names = FALSE)
  S1 <- tapply(Sm1, rep(1:nInt, n), mean)
  rm(list=c("Sm1"))

  ## rmst
  rmst0 <- cumsum(cData$timeIntLength * S0)
  rmst1 <- cumsum(cData$timeIntLength * S1)

  return(data.frame(S0=S0, S1=S1, rmst0=rmst0, rmst1=rmst1))

}


#' Algorithms for simulated data
#' @param nInt number of time intervals for coarsening the data
#' @param estimand "risk" or "rmst"
#' @param algorithm "TMLE", "AIPW", "IPW", "cox" or "weightedCox"

algorithmSim <- function(treatment, outcome, time,
                         survHaz, cenHaz, treatProb,
                         simOutcome, simTime, nInt,
                         estimand, algorithm){

  ## coarsen data parameters
  cData <- coarseData(time=time, outcome=outcome, nInt=nInt)

  if(estimand == "risk"){
    tau <- cData$timeIntMidPoint
  }else if(estimand == "rmst"){
    tau <- cData$timeIntMidPoint
    tau <- tau[-c(1, length(tau))]
  }

  ## algorithm
  if(algorithm == "TMLE"){

    result <- estimateTMLE(treatment=treatment, eventObserved=outcome, time=time,
                           survHaz=survHaz, cenHaz=cenHaz, treatProb=treatProb, tau=tau,
                           timeIntMidPoint=cData$timeIntMidPoint, timeIntLength=cData$timeIntLength,
                           estimand=estimand, printIter=TRUE, printTau=TRUE, tempCompare=FALSE)

  }else if(algorithm == "AIPW"){

    result <- estimateAIPW(treatment=treatment, eventObserved=outcome, time=time,
                           survHaz=SurvHaz, cenHaz=CenHaz, treatProb=treatProb, tau=tau,
                           timeIntMidPoint=cData$timeIntMidPoint, timeIntLength=cData$timeIntLength,
                           estimand=estimand, printTau=TRUE)

  }else if(algorithm == "IPW"){

    result <- estimateIPW(treatment=treatment, eventObserved=outcome, time=time,
                          cenHaz=CenHaz, treatProb=treatProb, tau=tau,
                          timeIntMidPoint=cData$timeIntMidPoint, timeIntLength=cData$timeIntLength,
                          estimand=estimand, printTau=TRUE)

  }else if(algorithm == "cox"){

    result <- strataCox(treatment=treatment, eventObserved=outcome,
                        time=time, treatProb=treatProb, cenHaz=NULL,
                        timeIntMidPoint=cData$timeIntMidPoint, timeIntLength=cData$timeIntLength, breaks=cData$breaks,
                        nsim=10000, printSim=TRUE)

  }else if(algorithm == "weightedCox"){

    result <- strataCox(treatment=treatment, eventObserved=outcome, time=time,
                        treatProb=treatProb, cenHaz=cenHaz,
                        timeIntMidPoint=cData$timeIntMidPoint, timeIntLength=cData$timeIntLength, breaks=cData$breaks,
                        nsim=10000, printSim=TRUE)

  }

  ## output
  return(result)

}




























