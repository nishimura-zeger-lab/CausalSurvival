#' Estimate treatment probability
#'
#' @param covIdTreatProb Covariates id to include in modeling treatment probability. If NULL, then include all covariates in the model
#' @param treatProbEstimate Model for estimating treatment probability. Options currently include logistic LASSO
#' @param maxCohortSizeForFitting If the target or comparator cohort are larger than this number, they
#'                                 will be downsampled before fitting the propensity model. The model
#'                                 will be used to compute propensity scores for all subjects. The
#'                                 purpose of the sampling is to gain speed.
#' @param index_ls Index for cross-fitting;
#'                 Default is no cross-fitting, index_ls = NULL
#' @param crossFitNum For cross-fitting: random partition of subjects into XXX prediction sets of approximately the same size.
#'                    If crossFitNum = 1, then no cross-fitting (default)
#' @return A data frame with columns: id, TreatProb

estimateTreatProb <- function(id, treatment, covariates, covIdTreatProb=NULL,
                              treatProbEstimate="LASSO", maxCohortSizeForFitting=25000,
                              index_ls=NULL, crossFitNum=1){

  ## container
  ID <- TreatProb <- c()
  ## outcomes
  outcomes <- data.frame(rowId = id, y = treatment)
  ## covariates
  if (!is.null(covIdTreatProb)){
    covariates <- covariates[which(covariates$covariateId %in% covIdTreatProb), ]
  }

  for (i in crossFitNum){

    if (is.null(index_ls)){
      ## training and testing sets are both the entire dataset
      idx_train <- idx_test <- id
    }else{
      ## training and testing sets
      idx_test <- index_ls[[i]]
      idx_train <- setdiff(id, idx_test)
    }

    ## prepare training set into appropriate format
    outcomes_train <- outcomes[which(outcomes$rowId %in% idx_train), ]
    covariates_train <- covariates[which(covariates$rowId %in% idx_train), ]
    ## downsize
    set.seed(0)
    targetRowIds <- outcomes_train$rowId[which(outcomes_train$y == 1)]
    targetRowIds <- sample(targetRowIds, size=min(maxCohortSizeForFitting, sum(outcomes_train$y == 1)), replace = FALSE)
    comparatorRowIds <- outcomes_train$rowId[which(outcomes_train$y == 0)]
    comparatorRowIds <- sample(comparatorRowIds, size=min(maxCohortSizeForFitting, sum(outcomes_train$y == 0)), replace = FALSE)

    outcomes_train_sub <- outcomes_train[which(outcomes_train$rowId %in% c(targetRowIds, comparatorRowIds)), ]
    covariates_train_sub <- covariates_train[which(covariates_train$rowId %in% c(targetRowIds, comparatorRowIds)), ]

    # parameters
    floatingPoint <- getOption("floatingPoint")
    if (is.null(floatingPoint)) {
      floatingPoint <- 64
    }
    # convert to data structure
    cyclopsData <- Cyclops::convertToCyclopsData(outcomes_train_sub, covariates_train_sub, modelType = "lr", quiet = TRUE, floatingPoint = floatingPoint)


    # prior
    prior = Cyclops::createPrior("laplace", exclude = c(0), useCrossValidation = TRUE)
    control = Cyclops::createControl(noiseLevel = "silent", cvType = "auto", seed = 1, tolerance = 2e-07, cvRepetitions = 10, startingVariance = 0.01)


    ## model: LASSO
    cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData, prior = prior, control = control)


    ## adjust intercept to full-training dataset
    y.bar <- mean(outcomes_train_sub$y)
    y.odds <- y.bar/(1 - y.bar)
    y.bar.new <- mean(outcomes_train$y)
    y.odds.new <- y.bar.new/(1 - y.bar.new)
    delta <- log(y.odds) - log(y.odds.new)
    cfs <- coef(cyclopsFit)
    cfs[1] <- cfs[1] - delta
    cyclopsFit$estimation$estimate[1] <- cfs[1]

    ## clear workspace
    rm(list=c("outcomes_train", "outcomes_train_sub", "covariates_train", "covariates_train_sub", "cyclopsData", "floatingPoint"))


    ## testing set
    outcomes_test <- outcomes[which(outcomes$rowId %in% idx_test), ]
    covariates_test <- covariates[which(covariates$rowId %in% idx_test), ]

    outcomes_test <- data.table::setDT(outcomes_test)
    covariates_test <- data.table::setDT(covariates_test)


    ## predict on testing set
    TreatProbTemp <- predict(cyclopsFit, newOutcomes = outcomes_test, newCovariates = covariates_test)


    ## store
    ID <- c(ID, idx_test[order(idx_test)])
    TreatProb <- c(TreatProb, TreatProbTemp)

    ## clear workspace
    rm(list=c("idx_test", "TreatProbTemp", "outcomes_test", "covariates_test", "cyclopsFit"))
  }

  ## result
  out <- data.frame(id=ID, TreatProb=TreatProb)
  out <- out[order(out$id), ]
  return(out)
}




#' Estimate discrete survival/censoring hazards
#'
#' @param id
#' @param treatment
#' @param eventObserved
#' @param time
#' @param covariates
#' @param covIdHaz Covariates id to include in modeling discrete hazards. If NULL, then include all covariates in the model
#' @param crossFitNum For cross-fitting: random partition of subjects into XXX prediction sets of approximately the same size.
#'                    If crossFitNum = 1 and index_ls = NULL, then no cross-fitting (default)
#'                    If crossFitNum = 1/2/3/4/5 and index_ls is not NULL, then perform crossFitNum'th fold cross-fit
#' @param index_ls Index for cross-fitting
#'                 Default is no cross-fitting, index_ls = NULL
#' @param timeEffect Functions of time in the discrete hazards model.
#'                   Options currently include "linear", "ns2", "ns3", "ns4", "ns5"
#' @param interactWithTime Data frame that include variables that interact with time in the hazards model
#' @param hazEstimate Model for estimating censoring hazards. Options currently include "glm", "ridge".
#'                    If hazEstimate = NULL, then must provide coef_Haz to have prediction of hazards
#' @param estimate_hazard "survival" or "censoring"
#' @param coef_H Supply output from pooledLogistic to have prediction of hazards
#' @param maxTimePredict Maximum time for prediction
#' @return A data frame with columns: ID, Haz1, Haz0

estimateHaz <- function(id, treatment, eventObserved, time, offset_t, offset_X, weight,
                           covariates, covIdHaz, sigma=exp(seq(log(1), log(0.01), length.out = 20)),
                           crossFitNum=1, index_ls=NULL, breaks, nInt,
                           timeEffect, evenKnot, penalizeTimeTreatment,
                           interactWithTime, hazEstimate, intercept,
                           estimate_hazard, getHaz, coef_H, robust, threshold){

  ## container
  if(getHaz){
    ID <- Haz1 <- Haz0 <- c()
  }else{
    coef_fit <- c()
    coef_var <- c()
    coef_robustVar <- c()
  }

  ## covariates to sparse matrix form, and delete unwanted covariates
  cov <- Matrix::sparseMatrix(i = covariates$rowId, j = covariates$covariateId, x = covariates$covariateValue, repr = "T")
  if (!is.null(covIdHaz)){
    cov <- cov[, covIdHaz]
  }

  ## coarsen data
  if(is.null(breaks)){
    cData <- coarseData(time=time, outcome=eventObserved, nInt=50)
    breaks <- cData$breaks
    timeIntMidPoint <- cData$timeIntMidPoint
    timeInt <- cData$timeInt
    rm(list=c("cData"))
  }else{
    timeIntMidPoint <- breaks[-length(breaks)] + (diff(breaks)/2)
    timeInt <- time
  }
  if(offset_X){
    offset_used_X <- log(as.double(as.character(cut(time, breaks=breaks, labels=diff(breaks)))))
  }else{
    offset_used_X <- NULL
  }



  for (i in crossFitNum){

    if (is.null(index_ls)){
      ## training and testing sets are both the entire dataset
      idx_train <- idx_test <- id
      ## offset
      offset_used_t <- offset_t
    }else{
      ## training and testing sets
      idx_test <- index_ls[[i]]
      idx_train <- setdiff(id, idx_test)
      idx_train <- idx_train[order(idx_train)]
      ## offset
      offset_used_t <- offset_t[[i]]
    }

    ## data
    X_baseline <- cbind(treatment, cov)[idx_train, ]
    if(is.null(dim(interactWithTime))){
      temporal_effect <- interactWithTime[idx_train]
    }else{
      temporal_effect <- interactWithTime[idx_train, ]
    }
    d_eventObserved <- eventObserved[idx_train]
    d_time <- timeInt[idx_train]
    d_offset_used_X <- offset_used_X[idx_train]

    ## reorder data
    time_indx <- order(d_time, 1-d_eventObserved, decreasing = TRUE)
    X_baseline <- X_baseline[time_indx, ]
    if(is.null(dim(interactWithTime))){
      temporal_effect <- temporal_effect[time_indx]
    }else{
      temporal_effect <- temporal_effect[time_indx, ]
    }
    d_eventObserved <- d_eventObserved[time_indx]
    d_time <- d_time[time_indx]
    d_offset_used_X <- d_offset_used_X[time_indx]
    rm(time_indx)


    ## model and prediction
    if(hazEstimate == "glm"){

      ## model: glm
      coef_Haz <- coef_pooled(X_baseline=X_baseline, temporal_effect=temporal_effect, time=d_time, eventObserved=d_eventObserved,
                              timeIntMidPoint=timeIntMidPoint, offset_t=offset_used_t, offset_X=d_offset_used_X, weight=weight,
                              timeEffect=timeEffect, is.temporal=TRUE, evenKnot=evenKnot, penalizeTimeTreatment=penalizeTimeTreatment,
                              intercept=intercept, estimate_hazard=estimate_hazard, sigma=NULL, maxiter=40, threshold=threshold, printIter=TRUE,
                              initial_coef=NULL, robust=robust)

      rm(list=c("X_baseline", "temporal_effect"))

    }else if(hazEstimate == "ridge"){

      ## model: ridge
      coef_Haz <- coef_ridge(X_baseline=X_baseline, temporal_effect=temporal_effect, eventObserved=d_eventObserved, time=d_time,
                             timeIntMidPoint=timeIntMidPoint, offset_t=offset_used_t, offset_X=d_offset_used_X, weight=weight,
                             timeEffect=timeEffect, is.temporal=TRUE, evenKnot=evenKnot, penalizeTimeTreatment=penalizeTimeTreatment,
                             intercept=intercept, estimate_hazard=estimate_hazard, sigma=sigma,
                             maxiter=40, threshold=threshold, printIter=TRUE)

      rm(list=c("X_baseline", "temporal_effect"))

    }

    ## prediction
    if(getHaz){

      ## parameter: maxTimeSplines for prediction if use ns(time, df=4)
      if(timeEffect == "ns" & estimate_hazard == "censoring"){
        indx_subset <- sapply(timeIntMidPoint, function(x) sum((d_time > x)*d_eventObserved+(d_time >= x)*(1 - d_eventObserved)), USE.NAMES = FALSE)
      }else if(timeEffect == "ns" & estimate_hazard == "survival"){
        indx_subset <- sapply(timeIntMidPoint, function(x) sum(d_time >= x), USE.NAMES = FALSE)
      }

    ## prepare dataset for prediction
    is.temporal <- TRUE
    if(intercept){
      X_baseline <- cbind(rep(1, dim(cov)[1]), rep(1, dim(cov)[1]), cov)[idx_test, , drop=FALSE]
    }else{
      X_baseline <- cbind(rep(1, dim(cov)[1]), cov)[idx_test, , drop=FALSE]
    }
    X_baseline <- Matrix::sparseMatrix(i = Matrix::summary(X_baseline)$i, j = Matrix::summary(X_baseline)$j,
                                       x = Matrix::summary(X_baseline)$x, repr = "R")
    if(is.null(interactWithTime)){
      temporal_effect <- interactWithTime[idx_test]
    }else{
      temporal_effect <- rep(1, length(idx_test))
    }
    d_offset_used_X <- offset_used_X[idx_test]

    if(is.null(temporal_effect) & !is.temporal){
      temporal_effect <- cbind(rep(0, dim(X_baseline)[1]), temporal_effect)
      nsBase <- NULL
    }else if(is.temporal & timeEffect == "linear"){
      temporal_effect <- cbind(rep(1, dim(X_baseline)[1]), temporal_effect)
      nsBase <- NULL
    }else if(timeEffect == "ns"){
      temporal_effect <- cbind(rep(1, dim(X_baseline)[1]), rep(1, dim(X_baseline)[1]), rep(1, dim(X_baseline)[1]),
                               rep(1, dim(X_baseline)[1]), rep(1, dim(X_baseline)[1]), temporal_effect,
                               temporal_effect, temporal_effect, temporal_effect, temporal_effect)
      if(evenKnot){
        nsBase <- splines::ns(timeIntMidPoint, df=5)
        if(hazEstimate == "ridge"){if(penalizeTimeTreatment){nsBase <- apply(nsBase, 2, function(x) x/sd(x))}}
      }else{
        nsBase <- splines::ns(timeIntMidPoint, knots=quantile(rep(timeIntMidPoint, times=indx_subset), probs=c(0.2, 0.4, 0.6, 0.8)))
        if(hazEstimate == "ridge"){if(penalizeTimeTreatment){nsBase <- apply(nsBase, 2, function(x) x/sd(rep(x, times=indx_subset)))}}
      }
    }


    if(is.null(coef_H)){
      ## prediction with derived coefficients
      Haz1temp <- predict_pooled(coef=coef_Haz$estimates, X_baseline=X_baseline,
                                 temporal_effect=temporal_effect, offset_t=offset_used_t, offset_X=d_offset_used_X,
                                 timeIntMidPoint=timeIntMidPoint, timeEffect=timeEffect, nsBase=nsBase)
      if(intercept){X_baseline[, 2] <- 0}else{X_baseline[, 1] <- 0}
      X_baseline <- Matrix::sparseMatrix(i = Matrix::summary(X_baseline)$i, j = Matrix::summary(X_baseline)$j,
                                         x = Matrix::summary(X_baseline)$x, repr = "R")
      if(is.null(interactWithTime)){
        temporal_effect <- interactWithTime[idx_test]
      }else{
        temporal_effect <- rep(0, length(idx_test))
      }

      if(is.null(temporal_effect) & !is.temporal){
        temporal_effect <- cbind(rep(0, dim(X_baseline)[1]), temporal_effect)
      }else if(is.temporal & timeEffect == "linear"){
        temporal_effect <- cbind(rep(1, dim(X_baseline)[1]), temporal_effect)
      }else if(timeEffect == "ns"){
        temporal_effect <- cbind(rep(1, dim(X_baseline)[1]), rep(1, dim(X_baseline)[1]), rep(1, dim(X_baseline)[1]),
                                 rep(1, dim(X_baseline)[1]), rep(1, dim(X_baseline)[1]),
                                 temporal_effect, temporal_effect, temporal_effect, temporal_effect, temporal_effect)
      }

      Haz0temp <- predict_pooled(coef=coef_Haz$estimates, X_baseline=X_baseline,
                                 temporal_effect=temporal_effect, offset_t=offset_used_t, offset_X=d_offset_used_X,
                                 timeIntMidPoint=timeIntMidPoint, timeEffect=timeEffect, nsBase=nsBase)

      rm(coef_Haz)

    }else{
      ## prediction with input coefficients
      Haz1temp <- predict_pooled(coef=coef_H[, i], X_baseline=X_baseline,
                                 temporal_effect=temporal_effect, offset_t=offset_used_t, offset_X=d_offset_used_X,
                                 timeIntMidPoint=timeIntMidPoint, timeEffect=timeEffect, nsBase=nsBase)

      if(intercept){X_baseline[, 2] <- 0}else{X_baseline[, 1] <- 0}
      X_baseline <- Matrix::sparseMatrix(i = Matrix::summary(X_baseline)$i, j = Matrix::summary(X_baseline)$j,
                                         x = Matrix::summary(X_baseline)$x, repr = "R")
      if(is.null(interactWithTime)){
        temporal_effect <- interactWithTime[idx_test]
      }else{
        temporal_effect <- rep(0, length(idx_test))
      }

      if(is.null(temporal_effect) & !is.temporal){
        temporal_effect <- cbind(rep(0, dim(X_baseline)[1]), temporal_effect)
      }else if(is.temporal & timeEffect == "linear"){
        temporal_effect <- cbind(rep(1, dim(X_baseline)[1]), temporal_effect)
      }else if(timeEffect == "ns"){
        temporal_effect <- cbind(rep(1, dim(X_baseline)[1]), rep(1, dim(X_baseline)[1]), rep(1, dim(X_baseline)[1]),
                                 rep(1, dim(X_baseline)[1]), rep(1, dim(X_baseline)[1]), temporal_effect, temporal_effect,
                                 temporal_effect, temporal_effect, temporal_effect)
      }
      Haz0temp <- predict_pooled(coef=coef_H[, i], X_baseline=X_baseline,
                                 temporal_effect=temporal_effect, offset_t=offset_used_t, offset_X=d_offset_used_X,
                                 timeIntMidPoint=timeIntMidPoint, timeEffect=timeEffect, nsBase=nsBase)
    }
    ## clear workspace
    rm(list=c("X_baseline", "temporal_effect", "idx_train"))

    ## store
    ID <- c(ID, rep(idx_test, each=length(timeIntMidPoint)))
    Haz1 <- c(Haz1, Haz1temp)
    Haz0 <- c(Haz0, Haz0temp)

    ## clear workspace
    rm(list=c("Haz1temp", "Haz0temp", "idx_test"))

    }else{
    ## no prediction, output coefficients
      coef_fit <- cbind(coef_fit, coef_Haz$estimates)
      if(hazEstimate == "glm"){
        coef_var <- cbind(coef_var, coef_Haz$var)
        if(robust){
        coef_robustVar <- cbind(coef_robustVar, coef_Haz$robustVar)
        }
        }
    }
  }

  ## result
  if(getHaz){
    out <- data.frame(ID=ID, Haz1=Haz1, Haz0=Haz0)
    out <- out[order(out$ID), ]
    rownames(out) <- NULL
    return(out)
  }else{
    if(hazEstimate == "glm"){
      if(robust){
      return(list(coef_fit=coef_fit, coef_var=coef_var, coef_robustVar=coef_robustVar))
      }else{
      return(list(coef_fit=coef_fit, coef_var=coef_var))
      }
    }else if(hazEstimate == "ridge"){
      return(coef_fit)
    }
  }
}




