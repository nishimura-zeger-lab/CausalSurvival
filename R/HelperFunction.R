## bound number to (0, 1)
bound01 <- function(x, r = 1e-10){
  xx <- x
  xx[x < r] <- r
  xx[x > 1-r] <- 1-r
  return(as.numeric(xx))
}


## bound number to (0, Inf)
bound <- function(x, r = 0.001){
  xx <- x
  xx[x < r] <- r
  return(as.numeric(xx))
}


#' Transform survival data from wide-format to long-format
#'
#' @param dwide Wide-format survival data with columns: time (observed time), eventObserved (Observed event), id
#' @param freq.time Map time interval to coarser intervals
#'
#' @export dlong A long-format survival data (with coarsening if freq.time > 1)
#'               with columns: t (time points), It, Jt, Rt, Lt (four indicator functions) and other covariates

transformData <- function(dwide, freq.time){
  ## Coarsen data
  if(freq.time > 1) dwide$time <- dwide$time %/% freq.time + 1

  ## number of subjects
  n <- dim(dwide)[1]
  ## maximum follow-up time
  maxtime <- max(dwide$time)

  ## time points for each subjects
  t <- rep(1:maxtime, n)

  ## Indicator variables for each time point (see `Cross-fitted estimator for survival analysis` for definition)
  Lt <- Rt <- rep(NA, n*maxtime)
  It <- Jt <- 1*(t == 1)
  for(i in 1:maxtime){
    Rt[t == i] <- (1 - dwide$eventObserved) * (dwide$time == i)
    Lt[t == i] <- dwide$eventObserved * (dwide$time == i)
    It[t == i] <- (dwide$time >= i)
    Jt[t == i] <- (dwide$time > i) * dwide$eventObserved + (dwide$time >= i) * (1 - dwide$eventObserved)
  }

  ## Long-format dataset
  dlong <- data.frame(dwide[as.numeric(gl(n, maxtime)), ], t = t, It, Jt, Rt, Lt)
  return(dlong)
}



