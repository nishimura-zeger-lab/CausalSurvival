## bound number to (0, 1)
bound01 <- function(x, r = 1e-15){
  xx <- x
  xx[x < r] <- r
  xx[x > 1-r] <- 1-r
  return(as.numeric(xx))
}


## bound number to (r, Inf)
bound <- function(x, r = 1e-7){
  xx <- x
  xx[x < r] <- r
  return(as.numeric(xx))
}


#' Transform survival data from wide-format to long-format
#'
#' @param time
#' @param eventObserved Indicator of outcome-of-interest
#' @param timeIntMidPoint Midpoints of time intervals
#' @param type "survival" or "censoring"
#'
#' @return A long-format survival data (with coarsening if freqTime > 1)
#'               with columns: rowId, stratumId (subject id), time and y
#'

transformData <- function(time, eventObserved, timeIntMidPoint, type){

  n <- length(time)
  maxTime <- length(timeIntMidPoint)
  t <- rep(timeIntMidPoint, n)
  longOut <- rep(NA, n*maxTime)
  valid <- 1*(t == timeIntMidPoint[1])
  stratumId <- as.numeric(gl(n, maxTime))

  if(type == "survival"){
  for(i in timeIntMidPoint){
    longOut[t == i] <- eventObserved * (time == i)
    valid[t == i] <- (time >= i)
    }
  }else if(type == "censoring"){
    for(i in timeIntMidPoint){
      longOut[t == i] <- (1 - eventObserved) * (time == i)
      valid[t == i] <- (time > i) * eventObserved + (time >= i) * (1 - eventObserved)
    }
  }

  longOut <- longOut[valid]
  t <- t[valid]
  stratumId <- stratumId[valid]

  dlong <- data.frame(rowId = 1:length(longOut), stratumId = stratumId, time = t, y = longOut)
  return(dlong)

}



#' Index for cross-fitting
#'
#' @param crossFitNum For cross-fitting: random partition of subjects into J prediction sets of approximately the same size.
#' @return

crossFit <- function(eventObserved, id, crossFitNum){

  ## divide data into XXX groups with equal percentage of events
  set.seed(08082021)
  n_folds <- crossFitNum
  ## ID for subjects with or without observed events
  index_event <- id[which(eventObserved==1)]
  index_noevent <- id[which(eventObserved==0)]
  nid_event <- length(index_event)
  nid_noevent <- length(index_noevent)
  ## divide data into J groups
  index_ls_event <- split(sample(index_event, size=nid_event, replace=FALSE), rep(1:n_folds, each=ceiling(nid_event/n_folds))[1:nid_event])
  index_ls_noevent <- split(sample(index_noevent, size=nid_noevent, replace=FALSE), rep(1:n_folds, each=ceiling(nid_noevent/n_folds))[1:nid_noevent])
  ## combine list
  index_ls <- lapply(1:n_folds, function(x) c(index_ls_event[[x]], index_ls_noevent[[x]]))

  ## result
  return(index_ls)
}



#' Coarsen data to non-uniform time intervals
#' @param nInt number of time intervals for coarsening the data
#' @export

coarsenData <- function(time, outcome, nInt=NULL){

  if(is.null(nInt)){nInt <- min(50, floor((sum(outcome)/10)))}

  probSeq <- seq(0, 1, length.out=nInt+1)[-c(1, nInt+1)]
  timeStrata <- floor(quantile(time[outcome == 1], probs = probSeq))

  breaks <- unname(c(0, timeStrata, max(time)))
  timeIntMidPoint <- breaks[-length(breaks)] + (diff(breaks)/2)
  timeIntLength <- diff(breaks)

  timeInt <- as.double(as.character(cut(time, breaks=breaks, labels=timeIntMidPoint)))

  return(list(breaks=breaks, timeIntMidPoint=timeIntMidPoint, timeIntLength=timeIntLength, timeInt=timeInt))

}


#' Run the algorithm
#'
#' @export

computeEstimator <- function(treatment, outcome, time, initial_survHaz, cenHaz, treatProb,
                             coarseningParam, estimand="both", method, hazMethod=NULL){

  if(is.null(coarseningParam)){
    coarsenedTime <- list(timeIntMidPoint=1:max(time), timeIntLength=rep(1, max(time)), breaks=c(0, 1:max(time)))
  }
    tau <- coarseningParam$timeIntMidPoint

    print("Start algorithm")

  ## algorithm
  if(method == "TMLE"){

    result <- estimateTMLE(treatment=treatment, eventObserved=outcome, time=time,
                           survHaz=initial_survHaz, cenHaz=cenHaz, treatProb=treatProb$TreatProb, tau=tau,
                           timeIntMidPoint=coarseningParam$timeIntMidPoint, timeIntLength=coarseningParam$timeIntLength,
                           estimand=estimand, printIter=FALSE, printTau=TRUE)

  }else if(method == "AIPW"){

    result <- estimateAIPW(treatment=treatment, eventObserved=outcome, time=time,
                           survHaz=initial_survHaz, cenHaz=cenHaz, treatProb=treatProb$TreatProb, tau=tau,
                           timeIntMidPoint=coarseningParam$timeIntMidPoint, timeIntLength=coarseningParam$timeIntLength,
                           estimand=estimand, printTau=TRUE)

  }else if(method == "IPW"){

    result <- estimateIPW(treatment=treatment, eventObserved=outcome, time=time,
                          cenHaz=cenHaz, treatProb=treatProb$TreatProb, tau=tau,
                          timeIntMidPoint=coarseningParam$timeIntMidPoint, timeIntLength=coarseningParam$timeIntLength,
                          estimand=estimand, printTau=TRUE)

  }else if(method == "cox"){

    result <- strataCox(treatment=treatment, eventObserved=outcome,
                        time=time, treatProb=treatProb$TreatProb, cenHaz=NULL,
                        timeIntMidPoint=coarseningParam$timeIntMidPoint,
                        timeIntLength=coarseningParam$timeIntLength, breaks=coarseningParam$breaks,
                        nsim=5000, printSim=TRUE, hazMethod=hazMethod)

  }else if(method == "weightedCox"){

    result <- strataCox(treatment=treatment, eventObserved=outcome, time=time,
                        treatProb=treatProb$TreatProb, cenHaz=cenHaz,
                        timeIntMidPoint=coarseningParam$timeIntMidPoint,
                        timeIntLength=coarseningParam$timeIntLength, breaks=coarseningParam$breaks,
                        nsim=5000, printSim=TRUE, hazMethod=hazMethod)

  }

    print("Complete algorithm")

    if(method %in% c("TMLE", "AIPW", "IPW") & estimand != "both"){
      result_final <- data.frame(estimand = result$estimand1 - result$estimand0, SE = result$SE,
                                 estimand1 = result$estimand1, estimand0=result$estimand0)
    }else if(method %in% c("TMLE", "AIPW", "IPW") & estimand == "both"){
      result_final <- data.frame(risk_diff = result$S1 - result$S0, SE_risk_diff = result$SE_S,
                                 S1 = result$S1, S0=result$S0,
                                 rmst_diff = result$rmst1 - result$rmst0, SE_rmst_diff = result$SE_rmst,
                                 rmst1 = result$rmst1, rmst0=result$rmst0)
    }else if(method  %in% c("cox", "weightedCox")){
      result_final <- data.frame(risk_diff = result$S1 - result$S0, SE_risk_diff = result$SE_S,
                                 S1 = result$S1, S0=result$S0,
                                 rmst_diff = result$rmst1 - result$rmst0, SE_rmst_diff = result$SE_rmst,
                                 rmst1 = result$rmst1, rmst0=result$rmst0)
    }

  return(result_final)

}
