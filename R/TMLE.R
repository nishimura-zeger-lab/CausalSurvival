#' Estimate cross-fitted TMLE of survival probability
#'
#' @param covariants Design matrix in triplet format (row index, col index, and value)
#' @param J For cross-fitting: random partition of subjects into J prediction sets of approximately the same size.
#' @param h.estimate Model for estimating nuisance parameter: survival hazards
#' @param gR.estimate Model for estimating nuisance parameter: censoring hazards
#' @param gA.estimate Model for estimating nuisance parameter: treatment probability
#' @export

estimateTMLEprob <- function(eventTime, censorTime, treatment, covariants,
                             J=5, h.estimate="glm",  gR.estimate="glm", gA.estimate="LASSO"){

  ## Indicator for event
  eventObserved <- ifelse(is.na(eventTime), 0, 1)
  ## Observed time
  censored <- is.na(eventTime)
  time <- eventTime
  time[censored] <- censorTime[censored]


  ## Get index for cross-fitting
  index_ls <- crossFit(eventObserved=eventObserved, J=J)


  ## transform data into long format




  ## Estimate cross-fitted nuisance parameter



  ## Update h



  ## S_tmle


}






#' Estimate cross-fitted nuisance parameter
#'
#' @param J For cross-fitting
#' @param freq.time Map time interval to coarser intervals
#' @export

estimateNuisance <- function(dlong, J=5, freq.time=90,
                             h.estimate="glm",  gR.estimate="glm", gA.estimate="LASSO"){

  ## container
  id <- m <- h1 <- h0 <- gR1 <- gR0 <- gA1 <- rep(0, length=)


  for (i in 1:J){
    ## training and testing sets
    idx_test <- index_ls[[i]]
    idx_train <- setdiff(unique(dlong$id), idx_test)
    d_test <- subset(dlong, id %in% idx_test)
    d_train <- subset(dlong, id %in% idx_train)

    ## Survival hazard: glm
    fitL <- glm(Lm ~ A * (m + age65 + cardiovascular + female + CHADS2),
                data = d_train, subset = Im == 1, family = binomial())
    ## predict
    h11 <- bound01(predict(fitL, newdata = mutate(d_test, A = 1), type = 'response'))
    h01 <- bound01(predict(fitL, newdata = mutate(d_test, A = 0), type = 'response'))

    ## store

    ## Censoring hazard


    ## propensity scores


  }
}


#' Transform survival data from wide-format to long-format
#'
#' @param dwide Wide-format survival data with columns: time (observed time), eventObserved (Observed event)
#' @param freq.time Map time interval to coarser intervals
#' @export dlong A long-format survival data (with coarsening if freq.time > 1)
#'
transformData <- function(dwide, freq.time){
  ## Transform a dataset from the short to the long form.
  if(freq.time > 1) dwide$time <- dwide$time %/% freq.time + 1

  n <- dim(dwide)[1]
  maxtime <- max(dwide$time)

  m <- rep(1:maxtime, n)
  Lm <- Rm <- rep(NA, n*maxtime)
  Im <- Jm <- 1*(m == 1)

  for(t in 1:maxtime){
    Rm[m == t] <- (1 - dwide$eventObserved) * (dwide$time == t)
    Lm[m == t] <- dwide$eventObserved * (dwide$time == t)
    Im[m == t] <- (dwide$time >= t)
    Jm[m == t] <- (dwide$time > t) * dwide$eventObserved + (dwide$time >= t) * (1 - dwide$eventObserved)
  }

  dlong <- data.frame(dwide[as.numeric(gl(n, maxtime)), ], m = m, Im, Jm, Rm, Lm)

  return(dlong)

}



#' Index for cross-fitting
#'
#' @param J For cross-fitting: random partition of subjects into J prediction sets of approximately the same size.
#' @export

crossFit <- function(eventObserved, J=5){

  ## divide data into J groups with equal percentage of events
  set.seed(08082021)
  n_folds <- J
  rowId <- 1:length(eventObserved)
  ## ID for subjects with or without observed events
  index_event <- rowId[which(eventObserved==1)]
  index_noevent <- rowId[which(eventObserved==0)]
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








