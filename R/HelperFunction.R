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

transformData <- function(dwide, freqTime, type){

  ## Coarsen data
  if(freqTime > 1) dwide$time <- dwide$time %/% freqTime + 1

  ## number of subjects
  n <- dim(dwide)[1]

  if(type == "survival"){

  ## maximum follow-up time
  maxtime <- max(dwide$time[dwide$eventObserved == 1])
  ## time points for each subjects
  t <- rep(1:maxtime, n)

  Lt <- rep(NA, n*maxtime)
  It <- 1*(t == 1)

  for(i in 1:maxtime){
    Lt[t == i] <- dwide$eventObserved * (dwide$time == i)
    It[t == i] <- (dwide$time >= i)
  }

  ## Long-format dataset
  dlong <- data.frame(dwide[as.numeric(gl(n, maxtime)), ], t = t, It, Lt)

  }else if(type == "censoring"){

    ## maximum follow-up time
    maxtime <- max(dwide$time[dwide$eventObserved == 0])
    ## time points for each subjects
    t <- rep(1:maxtime, n)

    Rt <- rep(NA, n*maxtime)
    Jt <- 1*(t == 1)

    for(i in 1:maxtime){
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





