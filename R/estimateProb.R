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

estimateTreatProb <- function(id, treatment, covariates, covIdTreatProb,
                              treatProbEstimate, maxCohortSizeForFitting,
                              index_ls=NULL, crossFitNum=1){

  ## container
  ID <- TreatProb <- c()
  ## outcomes
  outcomes <- data.frame(rowId = id, y = treatment)
  ## covariates
  if (!is.null(covIdTreatProb)){
    covariates <- covariates[which(covariates$j %in% covIdTreatProb), ]
  }

  for (i in 1:crossFitNum){
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
    covariates_train <- covariates[which(covariates$i %in% idx_train), ]
    ## downsize
    set.seed(0)
    targetRowIds <- outcomes_train$rowId[which(outcomes_train$y == 1)]
    targetRowIds <- sample(targetRowIds, size=min(maxCohortSizeForFitting, sum(outcomes_train$y == 1)), replace = FALSE)
    comparatorRowIds <- outcomes_train$rowId[which(outcomes_train$y == 0)]
    comparatorRowIds <- sample(comparatorRowIds, size=min(maxCohortSizeForFitting, sum(outcomes_train$y == 0)), replace = FALSE)

    outcomes_train_sub <- outcomes_train[which(outcomes_train$rowId %in% c(targetRowIds, comparatorRowIds)), ]
    covariates_train_sub <- covariates_train[which(covariates_train$i %in% c(targetRowIds, comparatorRowIds)), ]

    colnames(covariates_train_sub) <- c("rowId", "covariateId", "covariateValue")
    colnames(covariates_train) <- c("rowId", "covariateId", "covariateValue")


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
    covariates_test <- covariates[which(covariates$i %in% idx_test), ]

    outcomes_test <- data.table::setDT(outcomes_test)
    colnames(covariates_test) <- c("rowId", "covariateId", "covariateValue")



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
#'                    If crossFitNum = 1, then no cross-fitting (default)
#' @param index_ls Index for cross-fitting
#'                 Default is no cross-fitting, index_ls = NULL
#' @param timeEffect Functions of time in the discrete hazards model.
#'                   Options currently include "linear", "ns"
#' @param interactWithTime Data frame that include variables that interact with time in the hazards model
#' @param hazEstimate Model for estimating censoring hazards. Options currently include "glm", "ridge"
#' @param estimate_hazard "survival" or "censoring"
#' @return A data frame with columns: ID, Haz1, Haz0

estimateHaz <- function(id, treatment, eventObserved, time,
                           covariates, covIdHaz,
                           crossFitNum=1, index_ls=NULL,
                           timeEffect, interactWithTime, hazEstimate, estimate_hazard){

  ## container
  ID <- Haz1 <- Haz0 <- c()
  ## covariates to sparse matrix form, and delete unwanted covariates
  cov <- Matrix::sparseMatrix(i = covariates$i, j = covariates$j, x = covariates$val, repr = "T")
  if (!is.null(covIdHaz)){
    cov <- cov[, covIdHaz]
  }

  for (i in 1:crossFitNum){

    if (is.null(index_ls)){
      ## training and testing sets are both the entire dataset
      idx_train <- idx_test <- id
    }else{
      ## training and testing sets
      idx_test <- index_ls[[i]]
      idx_train <- setdiff(id, idx_test)
    }

    ## data
    X_baseline <- cbind(treatment, cov)[idx_train, ]
    if(is.null(dim(interactWithTime))){
      temporal_effect <- interactWithTime[idx_train]
    }else{
      temporal_effect <- interactWithTime[idx_train, ]
    }
    eventObserved <- eventObserved[idx_train]
    time <- time[idx_train]

    ## reorder data
    time_indx <- order(time, 1-eventObserved, decreasing = TRUE)
    X_baseline <- X_baseline[time_indx, ]
    if(is.null(dim(interactWithTime))){
      temporal_effect <- temporal_effect[time_indx]
    }else{
      temporal_effect <- temporal_effect[time_indx, ]
    }
    eventObserved <- eventObserved[time_indx]
    time <- time[time_indx]
    rm(time_indx)

    ## model and prediction
    if (hazEstimate == "glm"){

      ## model: glm
      coef_Haz <- coef_pooled(X_baseline=X_baseline, temporal_effect=temporal_effect, is.temporal=TRUE,
                                 timeEffect=timeEffect, eventObserved=eventObserved, time=time,
                                 estimate_hazard=estimate_hazard, sigma=NULL,
                                 maxiter=40, threshold=1e-14, printIter=TRUE)

      rm(list=c("X_baseline", "temporal_effect"))

    }else if(hazEstimate == "ridge"){

      ## model: ridge
      coef_Haz <- coef_ridge(X_baseline=X_baseline, temporal_effect=temporal_effect, is.temporal=TRUE,
                                timeEffect=timeEffect, eventObserved=eventObserved, time=time,
                                estimate_hazard=estimate_hazard, sigma=exp(seq(log(0.01), log(2), length.out = 100)),
                                maxiter=40, threshold=1e-14, printIter=TRUE)

      rm(list=c("X_baseline", "temporal_effect"))

    }

    ## prediction
    is.temporal <- TRUE
    X_baseline <- cbind(rep(1, dim(cov)[1]), rep(1, dim(cov)[1]), cov)[idx_test, , drop=FALSE]
    temporal_effect <- interactWithTime[idx_test]
    if(is.null(temporal_effect) & !is.temporal){
      temporal_effect <- cbind(rep(0, dim(X_baseline)[1]), temporal_effect)
    }else if(is.temporal & timeEffect == "linear"){
      temporal_effect <- cbind(rep(1, dim(X_baseline)[1]), temporal_effect)
    }else if(timeEffect == "ns"){
      temporal_effect <- cbind(rep(1, dim(X_baseline)[1]), rep(1, dim(X_baseline)[1]), rep(1, dim(X_baseline)[1]),
                               rep(1, dim(X_baseline)[1]), rep(1, dim(X_baseline)[1]), temporal_effect,
                               temporal_effect, temporal_effect, temporal_effect, temporal_effect)
    }

    ## parameter
    maxTime <- min(max(time[eventObserved == 1]), max(time[eventObserved == 0]))
    if(timeEffect == "ns" & estimate_hazard == "censoring"){
      maxTimeSplines <- max(time[eventObserved == 0])
    }else if(timeEffect == "ns" & estimate_hazard == "survival"){
      maxTimeSplines <- max(time[eventObserved == 1])
    }else{
      maxTimeSplines <- NULL
    }

    Haz1temp <- predict_pooled(coef=coef_Haz$estimates, X_baseline=X_baseline,
                                  temporal_effect=temporal_effect, timeEffect=timeEffect,
                                  maxTime=maxTime, maxTimeSplines=maxTimeSplines)

    X_baseline[, 2] <- 0
    Haz0temp <- predict_pooled(coef=coef_Haz$estimates, X_baseline=X_baseline,
                                  temporal_effect=temporal_effect, timeEffect=timeEffect,
                                  maxTime=maxTime, maxTimeSplines=maxTimeSplines)

    ## clear workspace
    rm(list=c("X_baseline", "temporal_effect", "idx_train", "coef_Haz"))

    ## store
    ID <- c(ID, rep(idx_test, each=maxTime))
    Haz1 <- c(Haz1, Haz1temp)
    Haz0 <- c(Haz0, Haz0temp)

    ## clear workspace
    rm(list=c("Haz1temp", "Haz0temp", "idx_test"))

  }

  ## result
  out <- data.frame(ID=ID, Haz1=Haz1, Haz0=Haz0)
  out <- out[order(out$ID), ]
  return(out)
}




