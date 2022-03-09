#' Get initial parameters for simulation
#'

estimateSimulationParams <- function(treatment, covariates, outcome, hazEstimate, hazMethod, seed){

  ################
  ## parameters ##
  ################

  ## cov
  cov <- Matrix::sparseMatrix(i = covariates$i, j = covariates$j, x = covariates$val, repr = "T")
  ## id
  rowId <- 1:length(outcome)

  #######################################
  ## variable selection from real data ##
  #######################################

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

  ############################
  ## hazards from real data ##
  ############################

  if(hazMethod == "twoStage"){

    ## coarsen data
    timeStrata <- floor(quantile(time[outcome == 1], probs = seq(0.02, 0.98, by=0.02)))
    breaks <- unname(c(0, timeStrata, max(time)))
    timeIntMidPoint <- breaks[-length(breaks)] + (diff(breaks)/2)
    timeInt <- as.double(as.character(cut(time, breaks=breaks, labels=breaks[-length(breaks)] + (diff(breaks)/2))))

    # dlong
    dlong <- transformData(dwide=data.frame(eventObserved=outcome, time=timeInt), timeIntMidPoint=timeIntMidPoint, type="survival")
    rownames(dlong) <- NULL
    dlong <- dlong[, c("Lt", "It", "t")]

    ## offset
    fit <- mgcv::bam(Lt ~ s(t, bs="ps"), family = binomial, subset = It == 1, data = dlong, method="REML")
    offset_t <- predict(fit, newdata = data.frame(t=timeIntMidPoint))
    rm(list=c("dlong", "fit"))

    ## other parameters
    timeEffect <- "linear"
    evenKnot <- NULL

  }else if(hazMethod == "ns"){

    offset_t <- NULL
    timeEffect <- "ns"
    evenKnot <- FALSE

  }

  coef_h <- estimateHaz(id=rowId, treatment=treatment, eventObserved=outcome, time=timeInt,
                        offset_t=offset_t, offset_X=FALSE, intercept=TRUE, breaks=breaks,
                        covariates=covariates, covIdHaz=cov_indx, crossFitNum=1, index_ls=NULL,
                        timeEffect=timeEffect, evenKnot=evenKnot, penalizeTimeTreatment=FALSE,
                        interactWithTime=treatment, hazEstimate="ridge", weight=NULL,
                        sigma=sigma, estimate_hazard=hazEstimate, getHaz=FALSE, coef_H=NULL,
                        robust=FALSE, threshold=1e-10)

  ############
  ## output ##
  ############

  return(list(coef_h=coef_h, cov_indx=cov_indx))

}



#' Simulate time-to-event data with survival and censoring hazards
#'
#' @param survHaz Output from estimateHaz
#' @param cenHaz Output from estimateHaz
#' @param maxTimeSurv maximum time for event
#' @param maxTimeCen maximum time for censoring
#'

simData <- function(survHaz, cenHaz, treatment, maxTime, timeIntMidPoint){

  ## parameters
  n <- length(treatment)

  survHaz_all <- survHaz$Haz1 * treatment[survHaz$ID] + survHaz$Haz0 * (1 - treatment[survHaz$ID])
  cenHaz_all <- cenHaz$Haz1 * treatment[cenHaz$ID] + cenHaz$Haz0 * (1 - treatment[cenHaz$ID])

  rS <- rbinom(length(survHaz_all), 1, survHaz_all)
  rG <- rbinom(length(cenHaz_all), 1, cenHaz_all)

  tS <- tapply(rS, rep(1:n, each=maxTime), function(x){which(x == 1)[1]})
  tG <- tapply(rG, rep(1:n, each=maxTime), function(x){which(x == 1)[1]})
  tG <- ifelse(is.na(tG), maxTime, tG)

  ObservedTime <- apply(cbind(tS, tG), 1, function(x){min(x, na.rm = TRUE)})
  ObservedEvent <- 1*(tS<=tG)
  ObservedEvent <- ifelse(is.na(ObservedEvent), 0, ObservedEvent)

  ObservedEvent <- as.double(ObservedEvent)
  ObservedTime <- as.double(as.character(factor(ObservedTime, labels = timeIntMidPoint)))

  return(list(ObservedEvent=ObservedEvent,
              ObservedTime=ObservedTime))

}


simData2 <- function(survHaz, cenHaz, treatment, maxTime){

  n <- length(treatment)

  survHaz_all <- survHaz$Haz1 * treatment[survHaz$ID] + survHaz$Haz0 * (1 - treatment[survHaz$ID])
  cenHaz_all <- cenHaz$Haz1 * treatment[cenHaz$ID] + cenHaz$Haz0 * (1 - treatment[cenHaz$ID])

  S_all <- unlist(tapply(1-survHaz_all, survHaz$ID, cumprod, simplify = FALSE), use.names = FALSE)
  G_all <- unlist(tapply(1-cenHaz_all, cenHaz$ID, cumprod, simplify = FALSE), use.names = FALSE)

  RandS <- rep(runif(n), each=maxTime)
  RandG <- rep(runif(n), each=maxTime)

  Ts <- tapply((S_all>RandS), survHaz$ID, sum, simplify = TRUE)+1
  Tg <- tapply((G_all>RandG), cenHaz$ID, sum, simplify = TRUE)+1
  Ts <- ifelse(Ts == maxTimeSurv + 1, NA, Ts)

  ObservedTime <- apply(cbind(Ts, Tg), 1, function(x){min(x, na.rm = TRUE)})
  ObservedEvent <- 1*(Ts<Tg)
  ObservedTime <- ifelse(ObservedTime == maxTimeCen + 1, maxTimeCen, ObservedTime)
  ObservedEvent <- ifelse(is.na(ObservedEvent), 0, ObservedEvent)

  ObservedEvent <- as.double(ObservedEvent)
  ObservedTime <- as.double(ObservedTime)

  return(list(ObservedEvent=ObservedEvent,
              ObservedTime=ObservedTime))

}

#' Calculate counterfactuals from simulated data
#'
#'

counterFactuals <- function(survHaz, maxTime, timeIntLength){

  ## paramters
  n <- dim(survHaz)[1]/maxTime

  ## surv prob
  Sm0 <- unlist(tapply(1-survHaz$Haz0, rep(1:n, each=maxTime), cumprod, simplify = FALSE), use.names = FALSE)
  S0 <- tapply(Sm0, rep(1:maxTime, n), mean)
  rm(list=c("Sm0"))
  Sm1 <- unlist(tapply(1-survHaz$Haz1, rep(1:n, each=maxTime), cumprod, simplify = FALSE), use.names = FALSE)
  S1 <- tapply(Sm1, rep(1:maxTime, n), mean)
  rm(list=c("Sm1"))

  ## rmst
  rmst0 <- cumsum(timeIntLength * S0)
  rmst1 <- cumsum(timeIntLength * S1)

  return(data.frame(S0=S0, S1=S1, rmst0=rmst0, rmst1=rmst1))

}








