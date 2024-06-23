#' Estimate cross-fitted TMLE of survival probability
#'
#' @param covariants Design matrix in triplet format (row index, col index, and value)
#' @param J For cross-fitting: random partition of subjects into J prediction sets of approximately the same size.
#' @param h.estimate Model for estimating nuisance parameter: survival hazards
#' @param gR.estimate Model for estimating nuisance parameter: censoring hazards
#' @param gA.estimate Model for estimating nuisance parameter: treatment probability
#' @export

estimateTMLEprob <- function(eventTime, censorTime, treatment, covariates, covariates.names,
                             J=5, h.estimate="glm", gR.estimate="glm", gA.estimate="LASSO",
                            includeCovariatesInHestimate="",includeCovariatesIngRestimate="", includeCovariatesIngAestimate="all",
                             freq.time=90){

  ## Indicator for event
  eventObserved <- ifelse(is.na(eventTime), 0, 1)
  ## Observed time
  censored <- is.na(eventTime)
  time <- eventTime
  time[censored] <- censorTime[censored]
  ## subject id
  id <- 1:length(time)


  ## Get index for cross-fitting
  index_ls <- crossFit(eventObserved=eventObserved, id=id, J=J)


  ## Wide-form dataset
  dwide <- data.frame(id=id, time=time, eventObserved=eventObserved, treatment=treatment)


  ## transform data into long format
  dlong <- transformData(dwide=dwide, freq.time=freq.time)


  ## Estimate cross-fitted nuisance parameter



  ## Update h



  ## S_tmle


}






#' Estimate cross-fitted nuisance parameter survival hazards
#'
#' @param dlong Long-format survival data from function transformData(dwide, freq.time)
#' @param J For cross-fitting
#' @param h.estimate Model for estimating nuisance parameter: survival hazards
#' @export

estimateNuisanceH <- function(dlong, J, h.estimate="glm"){

  ## container
  id <- t <- h1 <- h0 <- rep(0, length=)


  for (i in 1:J){
    ## training and testing sets
    idx_test <- index_ls[[i]]
    idx_train <- setdiff(unique(dlong$id), idx_test)
    d_test <- subset(dlong, id %in% idx_test)
    d_train <- subset(dlong, id %in% idx_train)

    ## Survival hazard: glm
    fitL <- glm(Lm ~ A * (t + age65 + cardiovascular + female + CHADS2),
                data = d_train, subset = Im == 1, family = binomial())
    ## predict
    h11 <- bound01(predict(fitL, newdata = mutate(d_test, A = 1), type = 'response'))
    h01 <- bound01(predict(fitL, newdata = mutate(d_test, A = 0), type = 'response'))

    ## store


  }
}


#' Estimate cross-fitted nuisance parameter: censoring hazards
#'
#' @param dlong Long-format survival data from function transformData(dwide, freq.time)
#' @param J For cross-fitting
#' @param gR.estimate Model for estimating nuisance parameter: censoring hazards
#' @export

estimateNuisanceGR <- function(dlong, J, gR.estimate="glm"){

  ## container
  id <- t <- gR1 <- gR0 <- rep(0, length=)


  for (i in 1:J){
    ## training and testing sets
    idx_test <- index_ls[[i]]
    idx_train <- setdiff(unique(dlong$id), idx_test)
    d_test <- subset(dlong, id %in% idx_test)
    d_train <- subset(dlong, id %in% idx_train)

    ## Censoring hazard



  }
}



#' Estimate cross-fitted nuisance parameter: treatment probability
#'
#' @param dlong Long-format survival data from function transformData(dwide, freq.time)
#' @param J For cross-fitting
#' @param gA.estimate Model for estimating nuisance parameter: treatment probability
#' @export

estimateNuisanceGA <- function(dlong, J, id, treatment, covariates, gA.estimate="LASSO"){

  ## container
  ID <- gA1 <- c()
  ## outcomes
  outcomes <- data.frame(rowId = id, y = treatment)
  ## covariates to sparse matrix
  cov <- Matrix::sparseMatrix(i = covariates$i, j = covariates$j, x = covariates$val, repr = "T")


  for (i in 1:J){
    ## training and testing sets
    idx_test <- index_ls[[i]]
    idx_train <- setdiff(unique(dlong$id), idx_test)


    ## prepare training set into appropriate format
    outcomes_train <- outcomes[which(outcomes$rowId %in% idx_train), ]
    covariates_train <- covariates[which(covariates$i %in% idx_train), ]
    colnames(covariates_train) <- c("rowId", "covariateId", "covariateValue")
    # parameters
    floatingPoint <- getOption("floatingPoint")
    if (is.null(floatingPoint)) {
      floatingPoint <- 64
    }
    # convert to data structure
    cyclopsData <- Cyclops::convertToCyclopsData(outcomes_train, covariates_train, modelType = "lr", quiet = TRUE, floatingPoint = floatingPoint)


    # prior
    prior = createPrior("laplace", exclude = c(0), useCrossValidation = TRUE)
    control = createControl(noiseLevel = "silent", cvType = "auto", seed = 1, tolerance = 2e-07, cvRepetitions = 10, startingVariance = 0.01)


    ## model: logistic regression
    cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData, prior = prior, control = control)


    ## predict on testing set
    cov_test <- cov[idx_test,]
    c <- coef(cyclopsFit)
    gAtemp <- as.matrix(c %*% t(as.matrix(cbind(rep(1, length=length(idx_test)), cov_test))))
    gAtemp <- as.vector(gAtemp)
    gAtemp <- exp(gAtemp)/(1+exp(gAtemp))


    ## store
    ID <- c(ID, idx_test)
    gA1 <- c(gA1, gAtemp)
  }

    ## result
    out <- data.frame(id=ID, gA1=gA1)
    return(out)
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



#' Index for cross-fitting
#'
#' @param J For cross-fitting: random partition of subjects into J prediction sets of approximately the same size.
#' @export

crossFit <- function(eventObserved, id, J=5){

  ## divide data into J groups with equal percentage of events
  set.seed(08082021)
  n_folds <- J
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








