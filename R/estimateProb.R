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

estimateTreatProb <- function(id, treatment, covariates, covIdTreatProb, treatProbEstimate, maxCohortSizeForFitting, index_ls=NULL, crossFitNum=1){

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




#' Estimate discrete censoring hazards
#'
#' @param dlong Long-format survival data from function transformData(dwide, freqTime)
#' @param covIdCenHaz Covariates id to include in modeling discrete censoring hazards. If NULL, then include all covariates in the model
#' @param crossFitNum For cross-fitting: random partition of subjects into XXX prediction sets of approximately the same size.
#'                    If crossFitNum = 1, then no cross-fitting (default)
#' @param index_ls Index for cross-fitting
#'                 Default is no cross-fitting, index_ls = NULL
#' @param timeEffect Functions of time in the discrete censoring hazards model.
#'                   Options currently include "linear", "ns"
#' @param interactWithTime Data frame that include variables that interact with time in the censoring hazards model
#' @param cenHazEstimate Model for estimating censoring hazards. Options currently include "LASSO", "glm"
#' @return A data frame with columns: ID, CenHaz1, CenHaz0

estimateCenHaz <- function(dlong, covariates, covIdCenHaz, crossFitNum=1, index_ls=NULL, timeEffect, interactWithTime, cenHazEstimate){

  ## container
  ID <- CenHaz1 <- CenHaz0 <- c()
  ## covariates to sparse matrix form, and delete unwanted covariates
  cov <- Matrix::sparseMatrix(i = covariates$i, j = covariates$j, x = covariates$val, repr = "T")
  if (!is.null(covIdCenHaz)){
    cov <- cov[, covIdCenHaz]
  }
  ## parameter
  maxTime <- min(max(dlong$time[dlong$eventObserved == 1]), max(dlong$time[dlong$eventObserved == 0]))


  for (i in 1:crossFitNum){

    if (is.null(index_ls)){
      ## training and testing sets are both the entire dataset
      idx_train <- idx_test <- unique(dlong$id)
    }else{
      ## training and testing sets
      idx_test <- index_ls[[i]]
      idx_train <- setdiff(unique(dlong$id), idx_test)
    }

    ## model and prediction
    if (cenHazEstimate == "glm"){
      ## data
      X_baseline <- cbind(dlong$treatment[dlong$t==1], cov)[idx_train, ]
      if(is.null(dim(interactWithTime))){
        temporal_effect <- interactWithTime[idx_train]
      }else{
        temporal_effect <- interactWithTime[idx_train, ]
      }
      eventObserved <- (dlong$eventObserved[dlong$t==1])[idx_train]
      time <- (dlong$time[dlong$t==1])[idx_train]
      estimate_hazard <- "censoring"

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

      ## model: glm
      coef_CenHaz <- coef_pooled(X_baseline=X_baseline, is.temporal=TRUE, temporal_effect=temporal_effect,
                                 timeEffect=timeEffect, eventObserved=eventObserved, time=time,
                                 estimate_hazard=estimate_hazard, maxiter=40, threshold=1e-14)

      ## prediction
      is.temporal <- TRUE

      X_baseline <- cbind(rep(1, dim(cov)[1]), rep(1, dim(cov)[1]), cov)
      if(is.null(temporal_effect) & !is.temporal){
        temporal_effect <- cbind(rep(0, dim(X_baseline)[1]), temporal_effect)
      }else if(is.temporal & timeEffect == "linear"){
        temporal_effect <- cbind(rep(1, dim(X_baseline)[1]), temporal_effect)
      }else if(timeEffect == "ns"){
        temporal_effect <- cbind(rep(1, dim(X_baseline)[1]), rep(1, dim(X_baseline)[1]), rep(1, dim(X_baseline)[1]),
                                 rep(1, dim(X_baseline)[1]), rep(1, dim(X_baseline)[1]), temporal_effect,
                                 temporal_effect, temporal_effect, temporal_effect, temporal_effect)
      }

      CenHaz1temp <- predict_pooled(coef=coef_CenHaz$estimates, X_baseline=X_baseline[idx_test, , drop=FALSE],
                                    temporal_effect=temporal_effect[idx_test, , drop=FALSE], timeEffect=timeEffect,
                                    maxTime=maxTime)

      X_baseline <- cbind(rep(1, dim(cov)[1]), rep(0, dim(cov)[1]), cov)
      CenHaz0temp <- predict_pooled(coef=coef_CenHaz$estimates, X_baseline=X_baseline[idx_test, , drop=FALSE],
                                    temporal_effect=temporal_effect[idx_test, , drop=FALSE], timeEffect=timeEffect,
                                    maxTime=maxTime)

      ## clear workspace
      rm(list=c("X_baseline", "temporal_effect", "eventObserved", "time", "idx_train", "coef_CenHaz"))

    }else if(cenHazEstimate == "LASSO"){
      ## data
      inx_train_LASSO <- which(dlong$id %in% idx_train & dlong$Jt == 1)
      outcomes <- data.frame(rowId = 1:dim(dlong[inx_train_LASSO, ])[1], y = dlong$Rt[inx_train_LASSO])
      if(timeEffect == "linear"){
        cov_long <- cbind(dlong$treatment[inx_train_LASSO], dlong$t[inx_train_LASSO], interactWithTime[dlong$id[inx_train_LASSO]], cov[dlong$id[inx_train_LASSO], ])
      }else if(timeEffect == "cubic"){
        cov_long <- cbind(dlong$treatment[inx_train_LASSO],
                          dlong$t[inx_train_LASSO], (dlong$t[inx_train_LASSO])^2, (dlong$t[inx_train_LASSO])^3,
                          dlong$t[inx_train_LASSO] * interactWithTime[dlong$id[inx_train_LASSO]],
                          (dlong$t[inx_train_LASSO])^2 * interactWithTime[dlong$id[inx_train_LASSO]],
                          (dlong$t[inx_train_LASSO])^3 * interactWithTime[dlong$id[inx_train_LASSO]],
                          cov[dlong$id[inx_train_LASSO], ])
      }else if(timeEffect == "ns"){
        cov_long <- cbind(dlong$treatment[inx_train_LASSO], ns(dlong$t[inx_train_LASSO], df=5),
                          ns(dlong$t[inx_train_LASSO], df=5) * interactWithTime[dlong$id[inx_train_LASSO]],
                          cov[dlong$id[inx_train_LASSO], ])
      }
      covariates <- as.data.frame(summary(cov_long))
      colnames(covariates) <- c("rowId", "covariateId", "covariateValue")
      ## convert to CyclopsData
      cyclopsData <- Cyclops::convertToCyclopsData(outcomes, covariates, modelType = "lr", quiet = TRUE)

      ## clear workspace
      rm(list=c("inx_train_LASSO", "outcomes", "cov_long", "covariates"))

      # prior
      if(timeEffect == "linear"){
        prior = Cyclops::createPrior("laplace", exclude = 1:(2+!is.null(interactWithTime)), useCrossValidation = TRUE)
      }else if(timeEffect == "cubic"){
        prior = Cyclops::createPrior("laplace", exclude = 1:(4+3*!is.null(interactWithTime)), useCrossValidation = TRUE)
      }else if(timeEffect == "ns"){
        prior = Cyclops::createPrior("laplace", exclude = 1:(6+5*!is.null(interactWithTime)), useCrossValidation = TRUE)
      }
      control = Cyclops::createControl(noiseLevel = "silent", cvType = "auto", fold = 5, seed = 1, tolerance = 1e-06, cvRepetitions = 1, startingVariance = 0.01)

      ## model: LASSO
      cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData, prior = prior, control = control)

      ## clear workspace
      rm(list=c("cyclopsData", "prior", "control"))

      ## predict
      Intercept <- rep(1, length(idx_test))
      TreatmentIndi1 <- rep(1, length(idx_test))
      TreatmentIndi0 <- rep(0, length(idx_test))

      coef_CenHaz <- coef(cyclopsFit)

      if(timeEffect == "linear"){
        timeIndepLP1 <- rep(cbind(Intercept, TreatmentIndi1, cov[idx_test, ]) %*% coef_CenHaz[-3:-(3+!is.null(interactWithTime))], each=maxTime)
        timeDepenLP1 <- cbind(rep(1:maxTime, length(idx_test)), rep(1:maxTime, length(idx_test)) * TreatmentIndi1 * interactWithTime) %*% coef_CenHaz[3:(3+!is.null(interactWithTime))]
      }else if(timeEffect == "cubic"){
        timeIndepLP1 <- rep(cbind(Intercept, TreatmentIndi1, cov[idx_test, ]) %*% coef_CenHaz[-3:-(3+3*!is.null(interactWithTime))], each=maxTime)
        timeDepenLP1 <- cbind(rep(1:maxTime, length(idx_test)), (rep(1:maxTime, length(idx_test)))^2, (rep(1:maxTime, length(idx_test)))^3,
                              rep(1:maxTime, length(idx_test)) * TreatmentIndi1 * interactWithTime, (rep(1:maxTime, length(idx_test)))^2 * TreatmentIndi1 * interactWithTime,
                              (rep(1:maxTime, length(idx_test)))^3 * TreatmentIndi1 * interactWithTime) %*% coef_CenHaz[3:(3+3*!is.null(interactWithTime))]
      }else if(timeEffect == "ns"){
        timeIndepLP1 <- rep(cbind(Intercept, TreatmentIndi1, cov[idx_test, ]) %*% coef_CenHaz[-3:-(3+5*!is.null(interactWithTime))], each=maxTime)
        timeDepenLP1 <- cbind(ns(rep(1:maxTime, length(idx_test)), df=5), ns(rep(1:maxTime, length(idx_test)), df=5) * TreatmentIndi1 * interactWithTime) %*% coef_CenHaz[3:(3+5*!is.null(interactWithTime))]
      }

      LP1 <- timeIndepLP1 + timeDepenLP1
      CenHaz1temp <- exp(LP1)/(1+exp(LP1))

      if(timeEffect == "linear"){
        timeIndepLP0 <- rep(cbind(Intercept, TreatmentIndi0, cov[idx_test, ]) %*% coef_CenHaz[-3:-(3+!is.null(interactWithTime))], each=maxTime)
        timeDepenLP0 <- cbind(rep(1:maxTime, length(idx_test)), rep(1:maxTime, length(idx_test)) * TreatmentIndi0 * interactWithTime) %*% coef_CenHaz[3:(3+!is.null(interactWithTime))]
      }else if(timeEffect == "cubic"){
        timeIndepLP0 <- rep(cbind(Intercept, TreatmentIndi0, cov[idx_test, ]) %*% coef_CenHaz[-3:-(3+3*!is.null(interactWithTime))], each=maxTime)
        timeDepenLP0 <- cbind(rep(1:maxTime, length(idx_test)), (rep(1:maxTime, length(idx_test)))^2, (rep(1:maxTime, length(idx_test)))^3,
                              rep(1:maxTime, length(idx_test)) * TreatmentIndi0 * interactWithTime, (rep(1:maxTime, length(idx_test)))^2 * TreatmentIndi0 * interactWithTime,
                              (rep(1:maxTime, length(idx_test)))^3 * TreatmentIndi0 * interactWithTime) %*% coef_CenHaz[3:(3+3*!is.null(interactWithTime))]
      }else if(timeEffect == "ns"){
        timeIndepLP0 <- rep(cbind(Intercept, TreatmentIndi0, cov[idx_test, ]) %*% coef_CenHaz[-3:-(3+5*!is.null(interactWithTime))], each=maxTime)
        timeDepenLP0 <- cbind(ns(rep(1:maxTime, length(idx_test)), df=5), ns(rep(1:maxTime, length(idx_test)), df=5) * TreatmentIndi0 * interactWithTime) %*% coef_CenHaz[3:(3+5*!is.null(interactWithTime))]
      }

      LP0 <- timeIndepLP0 + timeDepenLP0
      CenHaz0temp <- exp(LP0)/(1+exp(LP0))

      rm(list=c("cyclopsFit", "LP1", "LP0", "idx_train", "coef_CenHaz", "Intercept",
                "TreatmentIndi1", "TreatmentIndi0", "timeIndepLP1", "timeDepenLP1",
                "timeIndepLP0", "timeDepenLP0"))
    }


    ## store
    ID <- c(ID, rep(idx_test, each=maxTime))
    CenHaz1 <- c(CenHaz1, CenHaz1temp)
    CenHaz0 <- c(CenHaz0, CenHaz0temp)

    ## clear workspace
    rm(list=c("CenHaz1temp", "CenHaz0temp", "idx_test"))
  }

  ## result
  out <- data.frame(ID=ID, CenHaz1=CenHaz1, CenHaz0=CenHaz0)
  out <- out[order(out$ID), ]
  return(out)
}



#' Estimate cross-fitted discrete survival hazards
#'
#' @param dlong Long-format survival data from function transformData(dwide, freqTime)
#' @param covIdSurvHaz Covariates id to include in modeling discrete survival hazards. If NULL, then include all covariates in the model
#' @param crossFitNum For cross-fitting: random partition of subjects into XXX prediction sets of approximately the same size.
#'                    If crossFitNum = 1, then no cross-fitting (default)
#' @param index_ls Index for cross-fitting
#'                 Default is no cross-fitting, index_ls = NULL
#' @param interactWithTime Data frame that include variables that interact with time in the censoring hazards model
#' @param survHazEstimate Model for estimating survival hazards. Options currently include logistic LASSO, glm
#' @return A data frame with columns: ID, SurvHaz1, SurvHaz0

estimateSurvHaz <- function(dlong, covariates, covIdSurvHaz, crossFitNum=1, index_ls=NULL, timeEffect, interactWithTime, survHazEstimate){

  ## container
  ID <- SurvHaz1 <- SurvHaz0 <- c()
  ## covariates to sparse matrix form, and delete unwanted covariates
  cov <- Matrix::sparseMatrix(i = covariates$i, j = covariates$j, x = covariates$val, repr = "T")
  if (!is.null(covIdSurvHaz)){
    cov <- cov[, covIdSurvHaz]
  }
  ## parameter
  maxTime <- min(max(dlong$time[dlong$eventObserved == 1]), max(dlong$time[dlong$eventObserved == 0]))

  for (i in 1:crossFitNum){

    if (is.null(index_ls)){
      ## training and testing sets are both the entire dataset
      idx_train <- idx_test <- unique(dlong$id)
    }else{
      ## training and testing sets
      idx_test <- index_ls[[i]]
      idx_train <- setdiff(unique(dlong$id), idx_test)
    }

    ## model and prediction
    if (survHazEstimate == "glm"){
      ## data
      X_baseline <- cbind(dlong$treatment[dlong$t==1], cov)[idx_train, ]
      temporal_effect <- interactWithTime[idx_train, ]
      eventObserved <- (dlong$eventObserved[dlong$t==1])[idx_train]
      time <- (dlong$time[dlong$t==1])[idx_train]
      estimate_hazard <- "survival"

      ## reorder data
      time_indx <- order(time, decreasing = TRUE)
      X_baseline <- X_baseline[time_indx, ]
      temporal_effect <- temporal_effect[time_indx, ]
      eventObserved <- eventObserved[time_indx]
      time <- time[time_indx]
      rm(time_indx)

      ## model: glm
      coef_SurvHaz <- coef_pooled(X_baseline=X_baseline, temporal_effect=temporal_effect,
                                  is.temporal=TRUE, timeEffect=timeEffect,
                                  eventObserved=eventObserved, time=time,
                                  estimate_hazard=estimate_hazard, maxiter=40, threshold=1e-14)
      ## prediction
      is.temporal <- TRUE

      rm(list=c("X_baseline", "temporal_effect", "eventObserved", "time"))

      X_baseline <- cbind(rep(1, dim(cov)[1]), rep(1, dim(cov)[1]), cov)
      if(is.null(temporal_effect) & !is.temporal){
        temporal_effect <- cbind(rep(0, dim(X_baseline)[1]), temporal_effect)
      }else if(is.temporal & timeEffect == "linear"){
        temporal_effect <- cbind(rep(1, dim(X_baseline)[1]), temporal_effect)
      }else if(timeEffect == "ns"){
        temporal_effect <- cbind(rep(1, dim(X_baseline)[1]), rep(1, dim(X_baseline)[1]), rep(1, dim(X_baseline)[1]),
                                 rep(1, dim(X_baseline)[1]), rep(1, dim(X_baseline)[1]), temporal_effect,
                                 temporal_effect, temporal_effect, temporal_effect, temporal_effect)
      }

      SurvHaz1temp <- predict_pooled(coef=coef_SurvHaz$estimates, X_baseline=X_baseline[idx_test, , drop=FALSE],
                                     temporal_effect=temporal_effect[idx_test, , drop=FALSE], timeEffect=timeEffect,
                                     maxTime=maxTime)

      X_baseline <- cbind(rep(1, dim(cov)[1]), rep(0, dim(cov)[1]), cov)
      SurvHaz0temp <- predict_pooled(coef=coef_SurvHaz$estimates, X_baseline=X_baseline[idx_test, , drop=FALSE],
                                     temporal_effect=temporal_effect[idx_test, , drop=FALSE], timeEffect=timeEffect,
                                     maxTime=maxTime)

      rm(list=c("X_baseline", "temporal_effect", "eventObserved", "time", "idx_train", "coef_SurvHaz"))

    }else if(survHazEstimate == "LASSO"){
      ## data
      inx_train_LASSO <- which(dlong$id %in% idx_train & dlong$It == 1)
      outcomes <- data.frame(rowId = 1:dim(dlong[inx_train_LASSO, ])[1], y = dlong$Lt[inx_train_LASSO])
      if(timeEffect == "linear"){
        cov_long <- cbind(dlong$treatment[inx_train_LASSO], dlong$t[inx_train_LASSO], interactWithTime[dlong$id[inx_train_LASSO]], cov[dlong$id[inx_train_LASSO], ])
      }else if(timeEffect == "cubic"){
        cov_long <- cbind(dlong$treatment[inx_train_LASSO],
                          dlong$t[inx_train_LASSO], (dlong$t[inx_train_LASSO])^2, (dlong$t[inx_train_LASSO])^3,
                          dlong$t[inx_train_LASSO] * interactWithTime[dlong$id[inx_train_LASSO]],
                          (dlong$t[inx_train_LASSO])^2 * interactWithTime[dlong$id[inx_train_LASSO]],
                          (dlong$t[inx_train_LASSO])^3 * interactWithTime[dlong$id[inx_train_LASSO]],
                          cov[dlong$id[inx_train_LASSO], ])
      }else if(timeEffect == "ns"){
        cov_long <- cbind(dlong$treatment[inx_train_LASSO], ns(dlong$t[inx_train_LASSO], df=5),
                          ns(dlong$t[inx_train_LASSO], df=5) * interactWithTime[dlong$id[inx_train_LASSO]],
                          cov[dlong$id[inx_train_LASSO], ])
      }
      covariates <- as.data.frame(summary(cov_long))
      colnames(covariates) <- c("rowId", "covariateId", "covariateValue")
      ## convert to CyclopsData
      cyclopsData <- Cyclops::convertToCyclopsData(outcomes, covariates, modelType = "lr", quiet = TRUE)

      rm(list=c("inx_train_LASSO", "outcomes", "cov_long", "covariates"))

      # prior
      if(timeEffect == "linear"){
        prior = Cyclops::createPrior("laplace", exclude = 1:(2+!is.null(interactWithTime)), useCrossValidation = TRUE)
      }else if(timeEffect == "cubic"){
        prior = Cyclops::createPrior("laplace", exclude = 1:(4+3*!is.null(interactWithTime)), useCrossValidation = TRUE)
      }else if(timeEffect == "ns"){
        prior = Cyclops::createPrior("laplace", exclude = 1:(6+5*!is.null(interactWithTime)), useCrossValidation = TRUE)
      }
      control = Cyclops::createControl(noiseLevel = "silent", cvType = "auto", fold = 5, seed = 1, tolerance = 1e-06, cvRepetitions = 1, startingVariance = 0.01)

      ## model
      cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData, prior = prior, control = control)

      rm(list=c("cyclopsData", "prior", "control"))

      ## predict
      Intercept <- rep(1, length(idx_test))
      TreatmentIndi1 <- rep(1, length(idx_test))
      TreatmentIndi0 <- rep(0, length(idx_test))

      coef_SurvHaz <- coef(cyclopsFit)

      if(timeEffect == "linear"){
        timeIndepLP1 <- rep(cbind(Intercept, TreatmentIndi1, cov[idx_test, ]) %*% coef_SurvHaz[-3:-(3+!is.null(interactWithTime))], each=maxTime)
        timeDepenLP1 <- cbind(rep(1:maxTime, length(idx_test)), rep(1:maxTime, length(idx_test)) * TreatmentIndi1 * interactWithTime) %*% coef_SurvHaz[3:(3+!is.null(interactWithTime))]
      }else if(timeEffect == "cubic"){
        timeIndepLP1 <- rep(cbind(Intercept, TreatmentIndi1, cov[idx_test, ]) %*% coef_SurvHaz[-3:-(5+3*!is.null(interactWithTime))], each=maxTime)
        timeDepenLP1 <- cbind(rep(1:maxTime, length(idx_test)), (rep(1:maxTime, length(idx_test)))^2, (rep(1:maxTime, length(idx_test)))^3,
                              rep(1:maxTime, length(idx_test)) * TreatmentIndi1 * interactWithTime, (rep(1:maxTime, length(idx_test)))^2 * TreatmentIndi1 * interactWithTime,
                              (rep(1:maxTime, length(idx_test)))^3 * TreatmentIndi1 * interactWithTime) %*% coef_SurvHaz[3:(5+3*!is.null(interactWithTime))]
      }else if(timeEffect == "ns"){
        timeIndepLP1 <- rep(cbind(Intercept, TreatmentIndi1, cov[idx_test, ]) %*% coef_SurvHaz[-3:-(7+5*!is.null(interactWithTime))], each=maxTime)
        timeDepenLP1 <- cbind(ns(rep(1:maxTime, length(idx_test)), df=5), ns(rep(1:maxTime, length(idx_test)), df=5) * TreatmentIndi1 * interactWithTime) %*% coef_SurvHaz[3:(7+5*!is.null(interactWithTime))]
      }

      LP1 <- timeIndepLP1 + timeDepenLP1
      SurvHaz1temp <- exp(LP1)/(1+exp(LP1))


      if(timeEffect == "linear"){
        timeIndepLP0 <- rep(cbind(Intercept, TreatmentIndi0, cov[idx_test, ]) %*% coef_SurvHaz[-3:-(3+!is.null(interactWithTime))], each=maxTime)
        timeDepenLP0 <- cbind(rep(1:maxTime, length(idx_test)), rep(1:maxTime, length(idx_test)) * TreatmentIndi0 * interactWithTime) %*% coef_SurvHaz[3:(3+!is.null(interactWithTime))]
      }else if(timeEffect == "cubic"){
        timeIndepLP0 <- rep(cbind(Intercept, TreatmentIndi0, cov[idx_test, ]) %*% coef_SurvHaz[-3:-(5+3*!is.null(interactWithTime))], each=maxTime)
        timeDepenLP0 <- cbind(rep(1:maxTime, length(idx_test)), (rep(1:maxTime, length(idx_test)))^2, (rep(1:maxTime, length(idx_test)))^3,
                              rep(1:maxTime, length(idx_test)) * TreatmentIndi0 * interactWithTime, (rep(1:maxTime, length(idx_test)))^2 * TreatmentIndi0 * interactWithTime,
                              (rep(1:maxTime, length(idx_test)))^3 * TreatmentIndi0 * interactWithTime) %*% coef_SurvHaz[3:(5+3*!is.null(interactWithTime))]
      }else if(timeEffect == "ns"){
        timeIndepLP0 <- rep(cbind(Intercept, TreatmentIndi0, cov[idx_test, ]) %*% coef_SurvHaz[-3:-(7+5*!is.null(interactWithTime))], each=maxTime)
        timeDepenLP0 <- cbind(ns(rep(1:maxTime, length(idx_test)), df=5), ns(rep(1:maxTime, length(idx_test)), df=5) * TreatmentIndi0 * interactWithTime) %*% coef_SurvHaz[3:(7+5*!is.null(interactWithTime))]
      }

      LP0 <- timeIndepLP0 + timeDepenLP0
      SurvHaz0temp <- exp(LP0)/(1+exp(LP0))

      ## clear workspace
      rm(list=c("cyclopsFit", "LP1", "LP0", "idx_train", "coef_SurvHaz", "Intercept",
                "TreatmentIndi1", "TreatmentIndi0", "timeIndepLP1", "timeDepenLP1",
                "timeIndepLP0", "timeDepenLP0"))
    }

    ## store
    ID <- c(ID, rep(idx_test, each=maxTime))
    SurvHaz1 <- c(SurvHaz1, SurvHaz1temp)
    SurvHaz0 <- c(SurvHaz0, SurvHaz0temp)

    ## clear workspace
    rm(list=c("SurvHaz1temp", "SurvHaz0temp", "idx_test"))
  }

  ## result
  out <- data.frame(ID=ID, SurvHaz1=SurvHaz1, SurvHaz0=SurvHaz0)
  out <- out[order(out$ID), ]
  return(out)
}




