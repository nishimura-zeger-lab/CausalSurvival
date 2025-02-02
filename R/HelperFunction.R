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
#' @param dwide Wide-format survival data with columns: time (observed time), eventObserved (Observed event), id
#' @param freqTime Coarsen observed time to XXX days intervals
#' @param type "survival" or "censoring"
#'
#' @return A long-format survival data (with coarsening if freqTime > 1)
#'               with columns: t (time points), It, Jt, Rt, Lt (four indicator functions) and other covariates

transformData <- function(dwide, timeIntMidPoint, type){

  n <- dim(dwide)[1]
  maxtime <- length(timeIntMidPoint)
  t <- rep(timeIntMidPoint, n)

  if(type == "survival"){

  Lt <- rep(NA, n*maxtime)
  It <- 1*(t == timeIntMidPoint[1])

  for(i in timeIntMidPoint){
    Lt[t == i] <- dwide$eventObserved * (dwide$time == i)
    It[t == i] <- (dwide$time >= i)
  }

  ## Long-format dataset
  dlong <- data.frame(dwide[as.numeric(gl(n, maxtime)), ], t = t, It, Lt)

  }else if(type == "censoring"){

    Rt <- rep(NA, n*maxtime)
    Jt <- 1*(t == timeIntMidPoint[1])

    for(i in timeIntMidPoint){
      Rt[t == i] <- (1 - dwide$eventObserved) * (dwide$time == i)
      Jt[t == i] <- (dwide$time > i) * dwide$eventObserved + (dwide$time >= i) * (1 - dwide$eventObserved)
    }

    ## Long-format dataset
    dlong <- data.frame(dwide[as.numeric(gl(n, maxtime)), ], t = t, Rt, Jt)

  }

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
#'

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
#'

estimateCounterfactSurvival <- function(treatment, outcome, time,
                                        survHaz, cenHaz, treatProb,
                                        coarsenedTime, estimand, method, hazMethod=NULL){

  if(is.null(coarsenedTime)){
    coarsenedTime <- list(timeIntMidPoint=1:max(time), timeIntLength=rep(1, max(time)), breaks=c(0, 1:max(time)))
  }
    tau <- coarsenedTime$timeIntMidPoint

  ## algorithm
  if(method == "TMLE"){

    result <- estimateTMLE(treatment=treatment, eventObserved=outcome, time=time,
                           survHaz=survHaz, cenHaz=cenHaz, treatProb=treatProb$TreatProb, tau=tau,
                           timeIntMidPoint=coarsenedTime$timeIntMidPoint, timeIntLength=coarsenedTime$timeIntLength,
                           estimand=estimand, printIter=FALSE, printTau=TRUE, tempCompare=FALSE)

  }else if(method == "AIPW"){

    result <- estimateAIPW(treatment=treatment, eventObserved=outcome, time=time,
                           survHaz=survHaz, cenHaz=cenHaz, treatProb=treatProb$TreatProb, tau=tau,
                           timeIntMidPoint=coarsenedTime$timeIntMidPoint, timeIntLength=coarsenedTime$timeIntLength,
                           estimand=estimand, printTau=TRUE)

  }else if(method == "IPW"){

    result <- estimateIPW(treatment=treatment, eventObserved=outcome, time=time,
                          cenHaz=cenHaz, treatProb=treatProb$TreatProb, tau=tau,
                          timeIntMidPoint=coarsenedTime$timeIntMidPoint, timeIntLength=coarsenedTime$timeIntLength,
                          estimand=estimand, printTau=TRUE)

  }else if(method == "cox"){

    result <- strataCox(treatment=treatment, eventObserved=outcome,
                        time=time, treatProb=treatProb$TreatProb, cenHaz=NULL,
                        timeIntMidPoint=coarsenedTime$timeIntMidPoint,
                        timeIntLength=coarsenedTime$timeIntLength, breaks=coarsenedTime$breaks,
                        nsim=5000, printSim=TRUE, hazMethod=hazMethod)

  }else if(method == "weightedCox"){

    result <- strataCox(treatment=treatment, eventObserved=outcome, time=time,
                        treatProb=treatProb$TreatProb, cenHaz=cenHaz,
                        timeIntMidPoint=coarsenedTime$timeIntMidPoint,
                        timeIntLength=coarsenedTime$timeIntLength, breaks=coarsenedTime$breaks,
                        nsim=5000, printSim=TRUE, hazMethod=hazMethod)

  }

    if(method %in% c("TMLE", "AIPW", "IPW")){
      result_final <- data.frame(estimand = result$estimand1 - result$estimand0, SE = result$SE)
    }else{
      result_final <- data.frame(estimand_risk = result$S1 - result$S0,
                                 SE_risk = result$SE_S,
                                 estimand_rmst = result$rmst1 - result$rmst0,
                                 SE_rmst = result$SE_rmst)
    }

  return(result_final)

}
