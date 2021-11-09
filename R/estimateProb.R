#' Estimate cross-fitted treatment probability
#'
#' @param covID.TreatProb Covariates id to include in modeling treatment probability. If NULL, then include all covariates in the model
#' @param TreatProb.estimate Model for estimating nuisance parameter: treatment probability. Options currently include logistic LASSO
#' @param maxCohortSizeForFitting If the target or comparator cohort are larger than this number, they
#'                                 will be downsampled before fitting the propensity model. The model
#'                                 will be used to compute propensity scores for all subjects. The
#'                                 purpose of the sampling is to gain speed.
#' @param index_ls Index for cross-fitting
#' @param crossFitnum For cross-fitting: random partition of subjects into XXX prediction sets of approximately the same size.
#' @return A data frame with columns: id, TreatProb

estimateTreatProb <- function(id, treatment, covariates, covID.TreatProb, TreatProb.estimate, maxCohortSizeForFitting, index_ls, crossFitnum){

  ## container
  ID <- TreatProb <- c()
  ## outcomes
  outcomes <- data.frame(rowId = id, y = treatment)
  ## covariates
  if (!is.null(covID.TreatProb)){
    covariates <- covariates[which(covariates$j %in% covID.TreatProb), ]
  }

  for (i in 1:crossFitnum){
    ## training and testing sets
    idx_test <- index_ls[[i]]
    idx_train <- setdiff(id, idx_test)


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
    cfs <- Cyclops::coef(cyclopsFit)
    cfs[1] <- cfs[1] - delta
    cyclopsFit$estimation$estimate[1] <- cfs[1]

    ## clear workspace
    rm(list=c("outcomes_train", "outcomes_train_sub", "covariates_train", "covariates_train_sub", "cyclopsData"))



    ## testing set
    outcomes_test <- outcomes[which(outcomes$rowId %in% idx_test), ]
    covariates_test <- covariates[which(covariates$i %in% idx_test), ]

    outcomes_test <- data.table::setDT(outcomes_test)
    colnames(covariates_test) <- c("rowId", "covariateId", "covariateValue")



    ## predict on testing set
    TreatProbTemp <- Cyclops::predict(cyclopsFit, newOutcomes = outcomes_test, newCovariates = covariates_test)


    ## store
    ID <- c(ID, idx_test[order(idx_test)])
    TreatProb <- c(TreatProb, TreatProbTemp)

    ## clear workspace
    rm(list=c("idx_test", "TreatProbTemp", "outcomes_test", "covariates_test", "cyclopsFit"))
  }

  ## result
  out <- data.frame(id=ID, TreatProb=TreatProb)
  return(out)
}



#' Estimate treatment probability (no cross-fitting)
#'
#' @param covID.TreatProb Covariates id to include in modeling treatment probability. If NULL, then include all covariates in the model
#' @param TreatProb.estimate Model for estimating nuisance parameter: treatment probability. Options currently include logistic LASSO
#' @param maxCohortSizeForFitting If the target or comparator cohort are larger than this number, they
#'                                 will be downsampled before fitting the propensity model. The model
#'                                 will be used to compute propensity scores for all subjects. The
#'                                 purpose of the sampling is to gain speed.
#' @return A data frame with columns: id, TreatProb

estimateTreatProb2 <- function(id, treatment, covariates, covID.TreatProb, TreatProb.estimate, maxCohortSizeForFitting){

  ## outcomes
  outcomes <- data.frame(rowId = id, y = treatment)
  ## covariates
  if (!is.null(covID.TreatProb)){
    covariates <- covariates[which(covariates$j %in% covID.TreatProb), ]
  }

  ## downsize
  set.seed(0)
  targetRowIds <- outcomes$rowId[which(outcomes$y == 1)]
  targetRowIds <- sample(targetRowIds, size=min(maxCohortSizeForFitting, sum(outcomes$y == 1)), replace = FALSE)
  comparatorRowIds <- outcomes$rowId[which(outcomes$y == 0)]
  comparatorRowIds <- sample(comparatorRowIds, size=min(maxCohortSizeForFitting, sum(outcomes$y == 0)), replace = FALSE)

  outcomes_sub <- outcomes[which(outcomes$rowId %in% c(targetRowIds, comparatorRowIds)), ]
  covariates_sub <- covariates[which(covariates$i %in% c(targetRowIds, comparatorRowIds)), ]

  colnames(covariates_sub) <- c("rowId", "covariateId", "covariateValue")
  colnames(covariates) <- c("rowId", "covariateId", "covariateValue")


  # parameters
  floatingPoint <- getOption("floatingPoint")
  if (is.null(floatingPoint)) {
    floatingPoint <- 64
  }
  # convert to data structure
  cyclopsData <- Cyclops::convertToCyclopsData(outcomes_sub, covariates_sub, modelType = "lr", quiet = TRUE, floatingPoint = floatingPoint)


  # prior
  prior = Cyclops::createPrior("laplace", exclude = c(0), useCrossValidation = TRUE)
  control = Cyclops::createControl(noiseLevel = "silent", cvType = "auto", seed = 1, tolerance = 2e-07, cvRepetitions = 10, startingVariance = 0.01)


  ## model: LASSO
  cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData, prior = prior, control = control)



  ## adjust intercept to full-training dataset
  y.bar <- mean(outcomes_sub$y)
  y.odds <- y.bar/(1 - y.bar)
  y.bar.new <- mean(outcomes$y)
  y.odds.new <- y.bar.new/(1 - y.bar.new)
  delta <- log(y.odds) - log(y.odds.new)
  cfs <- Cyclops::coef(cyclopsFit)
  cfs[1] <- cfs[1] - delta
  cyclopsFit$estimation$estimate[1] <- cfs[1]

  ## clear workspace
  rm(list=c("outcomes_sub", "covariates_sub", "cyclopsData"))



  ## predict on testing set
  outcomes <- data.table::setDT(outcomes)
  TreatProb <- Cyclops::predict(cyclopsFit, newOutcomes = outcomes, newCovariates = covariates)


  ## result
  out <- data.frame(id=id, TreatProb=TreatProb)
  return(out)
}



#' Estimate cross-fitted discrete censoring hazards
#'
#' @param dlong Long-format survival data from function transformData(dwide, freq.time)
#' @param covID.CenHaz Covariates id to include in modeling discrete censoring hazards. If NULL, then include all covariates in the model
#' @param crossFitnum For cross-fitting: random partition of subjects into XXX prediction sets of approximately the same size.
#' @param index_ls Index for cross-fitting
#' @param CenHaz.estimate Model for estimating censoring hazards. Options currently include logistic LASSO, glm
#' @return A data frame with columns: ID, CenHaz1, CenHaz0

estimateCenHaz <- function(dlong, covariates, covID.CenHaz, crossFitnum, index_ls, CenHaz.estimate){

  ## container
  ID <- CenHaz1 <- CenHaz0 <- c()
  ## covariates to sparse matrix form, and delete unwanted covariates
  cov <- Matrix::sparseMatrix(i = covariates$i, j = covariates$j, x = covariates$val, repr = "T")
  if (!is.null(covID.CenHaz)){
    cov <- cov[, covID.CenHaz]
  }


  for (i in 1:crossFitnum){
    ## training and testing sets
    idx_test <- index_ls[[i]]
    idx_train <- setdiff(unique(dlong$id), idx_test)

    ## model and prediction
    if (CenHaz.estimate == "glm"){
      ## data
      X_baseline <- cbind(dlong$treatment[dlong$t==1], cov)
      temporal_effect <- NULL
      eventObserved <- dlong$eventObserved[dlong$t==1]
      time <- dlong$time[dlong$t==1]
      id <- dlong$id[dlong$t==1]
      estimate_hazard <- "censoring"

      ## model: glm
      coef_CenHaz <- coef_pooled(X_baseline=X_baseline[idx_train, ], temporal_effect=temporal_effect[idx_train], is.temporal=TRUE,
                       eventObserved=eventObserved[idx_train], time=time[idx_train], id=id[idx_train],
                       estimate_hazard=estimate_hazard, maxiter=40, threshold=1e-8)
      ## prediction
      X_baseline <- cbind(rep(1, dim(cov)[1]), cov)
      temporal_effect <- NULL
      CenHaz1temp <- predict_pooled(coef=coef_CenHaz, X_baseline=X_baseline[idx_test, ],
                             temporal_effect=temporal_effect, is.temporal = TRUE, maxTime=max(time))
      CenHaz1temp <- bound01(CenHaz1temp)

      X_baseline <- cbind(rep(0, dim(cov)[1]), cov)
      temporal_effect <- NULL
      CenHaz0temp <- predict_pooled(coef=coef_CenHaz, X_baseline=X_baseline[idx_test, ],
                                    temporal_effect=temporal_effect, is.temporal = TRUE, maxTime=max(time))
      CenHaz0temp <- bound01(CenHaz0temp)

      ## clear workspace
      rm(list=c("X_baseline", "temporal_effect", "eventObserved", "time", "id", "idx_train", "coef_CenHaz"))

    }else if(CenHaz.estimate == "LASSO"){
      ## data
      inx_train_LASSO <- which(dlong$id %in% idx_train & dlong$Jt == 1)
      outcomes <- data.frame(rowId = 1:dim(dlong[inx_train_LASSO, ])[1], y = dlong$Rt[inx_train_LASSO])
      cov_long <- cbind(dlong$treatment[inx_train_LASSO], dlong$t[inx_train_LASSO], cov[dlong$id[inx_train_LASSO], ])
      covariates <- as.data.frame(summary(cov_long))
      colnames(covariates) <- c("rowId", "covariateId", "covariateValue")
      ## convert to CyclopsData
      cyclopsData <- Cyclops::convertToCyclopsData(outcomes, covariates, modelType = "lr", quiet = TRUE)

      ## clear workspace
      rm(list=c("inx_train_LASSO", "outcomes", "cov_long", "covariates"))

      # prior
      prior = Cyclops::createPrior("laplace", exclude = c(1, 2), useCrossValidation = TRUE)
      control = Cyclops::createControl(noiseLevel = "silent", cvType = "auto", fold = 10, seed = 1, tolerance = 1e-06, cvRepetitions = 2, startingVariance = 0.01)

      ## model
      cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData, prior = prior, control = control)

      ## clear workspace
      rm(list=c("cyclopsData", "prior", "control"))

      ## predict
      coef_CenHaz <- Cyclops::coef(cyclopsFit)
      LP1 <- rep(cbind(rep(1, length(idx_test)), rep(1, length(idx_test)), cov[idx_test, ]) %*% coef_CenHaz[-3], each=70)+rep(1:70, length(idx_test))*coef_CenHaz[3]
      CenHaz1temp <- exp(LP1)/(1+exp(LP1))
      CenHaz1temp <- bound01(CenHaz1temp)
      LP0 <- rep(cbind(rep(1, length(idx_test)), rep(0, length(idx_test)), cov[idx_test, ]) %*% coef_CenHaz[-3], each=70)+rep(1:70, length(idx_test))*coef_CenHaz[3]
      CenHaz0temp <- exp(LP0)/(1+exp(LP0))
      CenHaz0temp <- bound01(CenHaz0temp)

      ## clear workspace
      rm(list=c("cyclopsFit", "LP1", "LP0", "idx_train", "coef_CenHaz"))
    }


    ## store
    ID <- c(ID, rep(idx_test, each=70))
    CenHaz1 <- c(CenHaz1, CenHaz1temp)
    CenHaz0 <- c(CenHaz0, CenHaz0temp)

    ## clear workspace
    rm(list=c("CenHaz1temp", "CenHaz0temp", "idx_test"))
  }

  ## result
  out <- data.frame(ID=ID, CenHaz1=CenHaz1, CenHaz0=CenHaz0)
  return(out)
}



#' Estimate discrete censoring hazards (no cross-fitting)
#'
#' @param dlong Long-format survival data from function transformData(dwide, freq.time)
#' @param covID.CenHaz Covariates id to include in modeling discrete censoring hazards. If NULL, then include all covariates in the model
#' @param CenHaz.estimate Model for estimating censoring hazards. Options currently include logistic LASSO, glm
#' @return A data frame with columns: CenHaz1, CenHaz0

estimateCenHaz2 <- function(dlong, covariates, covID.CenHaz, CenHaz.estimate){

  ## container
  ID <- CenHaz1 <- CenHaz0 <- c()
  ## covariates to sparse matrix form, and delete unwanted covariates
  cov <- Matrix::sparseMatrix(i = covariates$i, j = covariates$j, x = covariates$val, repr = "T")
  if (!is.null(covID.CenHaz)){
    cov <- cov[, covID.CenHaz]
  }

  ## model and prediction
  if (CenHaz.estimate == "glm"){
    ## data
    X_baseline <- cbind(dlong$treatment[dlong$t==1], cov)
    temporal_effect <- NULL
    eventObserved <- dlong$eventObserved[dlong$t==1]
    time <- dlong$time[dlong$t==1]
    id <- dlong$id[dlong$t==1]
    estimate_hazard <- "survival"

    ## model: glm
    coef_CenHaz <- coef_pooled(X_baseline=X_baseline, temporal_effect=temporal_effect, is.temporal=TRUE,
                     eventObserved=eventObserved, time=time, id=id,
                     estimate_hazard=estimate_hazard, maxiter=40, threshold=1e-8)
    ## prediction
    X_baseline <- cbind(rep(1, dim(cov)[1]), cov)
    temporal_effect <- NULL
    CenHaz1 <- predict_pooled(coef=coef_CenHaz, X_baseline=X_baseline,
                                  temporal_effect=temporal_effect, is.temporal = TRUE, maxTime=max(time))
    CenHaz1 <- bound01(CenHaz1)

    X_baseline <- cbind(rep(0, dim(cov)[1]), cov)
    temporal_effect <- NULL
    CenHaz0 <- predict_pooled(coef=coef_CenHaz, X_baseline=X_baseline,
                                  temporal_effect=temporal_effect, is.temporal = TRUE, maxTime=max(time))
    CenHaz0 <- bound01(CenHaz0)

    ## clear workspace
    rm(list=c("X_baseline", "temporal_effect", "eventObserved", "time", "id", "coef_CenHaz"))

  }else if(CenHaz.estimate == "LASSO"){
    ## data
    inx_train_LASSO <- which(dlong$Jt == 1)
    outcomes <- data.frame(rowId = 1:dim(dlong[inx_train_LASSO, ])[1], y = dlong$Rt[inx_train_LASSO])
    cov_long <- cbind(dlong$treatment[inx_train_LASSO], dlong$t[inx_train_LASSO], cov[dlong$id[inx_train_LASSO], ])
    covariates <- as.data.frame(summary(cov_long))
    colnames(covariates) <- c("rowId", "covariateId", "covariateValue")
    ## convert to CyclopsData
    cyclopsData <- Cyclops::convertToCyclopsData(outcomes, covariates, modelType = "lr", quiet = TRUE)

    ## clear workspace
    rm(list=c("inx_train_LASSO", "outcomes", "cov_long", "covariates"))

    # prior
    prior = Cyclops::createPrior("laplace", exclude = c(1, 2), useCrossValidation = TRUE)
    control = Cyclops::createControl(noiseLevel = "silent", cvType = "auto", fold = 10, seed = 1, tolerance = 1e-06, cvRepetitions = 2, startingVariance = 0.01)

    ## model
    cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData, prior = prior, control = control)

    ## clear workspace
    rm(list=c("cyclopsData", "prior", "control"))

    ## predict
    coef_CenHaz <- Cyclops::coef(cyclopsFit)
    n <- length(dlong$id[dlong$t==1])
    LP1 <- rep(cbind(rep(1, n), rep(1, n), cov) %*% coef_CenHaz[-3], each=70)+rep(1:70, n)*coef_CenHaz[3]
    CenHaz1 <- exp(LP1)/(1+exp(LP1))
    CenHaz1 <- bound01(CenHaz1)
    LP0 <- rep(cbind(rep(1, n), rep(0, n), cov) %*% coef_CenHaz[-3], each=70)+rep(1:70, n)*coef_CenHaz[3]
    CenHaz0 <- exp(LP0)/(1+exp(LP0))
    CenHaz0 <- bound01(CenHaz0)

    ## clear workspace
    rm(list=c("cyclopsFit", "LP1", "LP0", "coef_CenHaz"))
  }

  ## result
  out <- data.frame(CenHaz1=CenHaz1, CenHaz0=CenHaz0)
  return(out)
}




#' Estimate cross-fitted discrete survival hazards
#'
#' @param dlong Long-format survival data from function transformData(dwide, freq.time)
#' @param covID.SurvHaz Covariates id to include in modeling discrete survival hazards. If NULL, then include all covariates in the model
#' @param crossFitnum For cross-fitting: random partition of subjects into XXX prediction sets of approximately the same size.
#' @param index_ls Index for cross-fitting
#' @param SurvHaz.estimate Model for estimating survival hazards. Options currently include logistic LASSO, glm
#' @return A data frame with columns: ID, SurvHaz1, SurvHaz0

estimateSurvHaz <- function(dlong, covariates, covID.SurvHaz, crossFitnum, SurvHaz.estimate){

  ## container
  ID <- SurvHaz1 <- SurvHaz0 <- c()
  ## covariates to sparse matrix form, and delete unwanted covariates
  cov <- Matrix::sparseMatrix(i = covariates$i, j = covariates$j, x = covariates$val, repr = "T")
  if (!is.null(covID.SurvHaz)){
    cov <- cov[, covID.SurvHaz]
  }


  for (i in 1:crossFitnum){
    ## training and testing sets
    idx_test <- index_ls[[i]]
    idx_train <- setdiff(unique(dlong$id), idx_test)

    ## model and prediction
    if (SurvHaz.estimate == "glm"){
      ## data
      X_baseline <- cbind(dlong$treatment[dlong$t==1], cov)
      temporal_effect <- NULL
      eventObserved <- dlong$eventObserved[dlong$t==1]
      time <- dlong$time[dlong$t==1]
      id <- dlong$id[dlong$t==1]
      estimate_hazard <- "survival"

      ## model: glm
      coef_SurvHaz <- coef_pooled(X_baseline=X_baseline[idx_train, ], temporal_effect=temporal_effect[idx_train], is.temporal=TRUE,
                                 eventObserved=eventObserved[idx_train], time=time[idx_train], id=id[idx_train],
                                 estimate_hazard=estimate_hazard, maxiter=40, threshold=1e-8)
      ## prediction
      X_baseline <- cbind(rep(1, dim(cov)[1]), cov)
      temporal_effect <- NULL
      SurvHaz1temp <- predict_pooled(coef=coef_CenHaz, X_baseline=X_baseline[idx_test, ],
                                    temporal_effect=temporal_effect, is.temporal = TRUE, maxTime=max(time))
      SurvHaz1temp <- bound01(SurvHaz1temp)

      X_baseline <- cbind(rep(0, dim(cov)[1]), cov)
      temporal_effect <- NULL
      SurvHaz0temp <- predict_pooled(coef=coef_SurvHaz, X_baseline=X_baseline[idx_test, ],
                                    temporal_effect=temporal_effect, is.temporal = TRUE, maxTime=max(time))
      SurvHaz0temp <- bound01(SurvHaz0temp)

      ## clear workspace
      rm(list=c("X_baseline", "temporal_effect", "eventObserved", "time", "id", "idx_train", "coef_SurvHaz"))

    }else if(SurvHaz.estimate == "LASSO"){
      ## data
      inx_train_LASSO <- which(dlong$id %in% idx_train & dlong$It == 1)
      outcomes <- data.frame(rowId = 1:dim(dlong[inx_train_LASSO, ])[1], y = dlong$Lt[inx_train_LASSO])
      cov_long <- cbind(dlong$treatment[inx_train_LASSO], dlong$t[inx_train_LASSO], cov[dlong$id[inx_train_LASSO], ])
      covariates <- as.data.frame(summary(cov_long))
      colnames(covariates) <- c("rowId", "covariateId", "covariateValue")
      ## convert to CyclopsData
      cyclopsData <- Cyclops::convertToCyclopsData(outcomes, covariates, modelType = "lr", quiet = TRUE)

      ## clear workspace
      rm(list=c("inx_train_LASSO", "outcomes", "cov_long", "covariates"))

      # prior
      prior = Cyclops::createPrior("laplace", exclude = c(1, 2), useCrossValidation = TRUE)
      control = Cyclops::createControl(noiseLevel = "silent", cvType = "auto", fold = 10, seed = 1, tolerance = 1e-06, cvRepetitions = 2, startingVariance = 0.01)

      ## model
      cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData, prior = prior, control = control)

      ## clear workspace
      rm(list=c("cyclopsData", "prior", "control"))

      ## predict
      coef_SurvHaz <- Cyclops::coef(cyclopsFit)
      LP1 <- rep(cbind(rep(1, length(idx_test)), rep(1, length(idx_test)), cov[idx_test, ]) %*% coef_SurvHaz[-3], each=70)+rep(1:70, length(idx_test))*coef_SurvHaz[3]
      SurvHaz1temp <- exp(LP1)/(1+exp(LP1))
      SurvHaz1temp <- bound01(SurvHaz1temp)
      LP0 <- rep(cbind(rep(1, length(idx_test)), rep(0, length(idx_test)), cov[idx_test, ]) %*% coef_SurvHaz[-3], each=70)+rep(1:70, length(idx_test))*coef_SurvHaz[3]
      SurvHaz0temp <- exp(LP0)/(1+exp(LP0))
      SurvHaz0temp <- bound01(SurvHaz0temp)

      ## clear workspace
      rm(list=c("cyclopsFit", "LP1", "LP0", "idx_train", "coef_SurvHaz"))
    }

    ## store
    ID <- c(ID, rep(idx_test, each=70))
    SurvHaz1 <- c(SurvHaz1, SurvHaz1temp)
    SurvHaz0 <- c(SurvHaz0, SurvHaz0temp)

    ## clear workspace
    rm(list=c("SurvHaz1temp", "SurvHaz0temp", "idx_test"))
  }

  ## result
  out <- data.frame(ID=ID, SurvHaz1=SurvHaz1, SurvHaz0=SurvHaz0)
  return(out)
}




