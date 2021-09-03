#' Estimate cross-fitted nuisance parameter: treatment probability
#'
#' @param J For cross-fitting
#' @param gA.estimate Model for estimating nuisance parameter: treatment probability
#' @param maxCohortSizeForFitting
#' @export

estimateNuisanceGA <- function(J, id, treatment, covariates, gA.estimate="LASSO", maxCohortSizeForFitting){

  ## container
  ID <- gA1 <- c()
  ## outcomes
  outcomes <- data.frame(rowId = id, y = treatment)


  for (i in 1:J){
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
    prior = createPrior("laplace", exclude = c(0), useCrossValidation = TRUE)
    control = createControl(noiseLevel = "silent", cvType = "auto", seed = 1, tolerance = 2e-07, cvRepetitions = 10, startingVariance = 0.01)


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



    ## testing set
    outcomes_test <- outcomes[which(outcomes$rowId %in% idx_test), ]
    covariates_test <- covariates[which(covariates$i %in% idx_test), ]

    outcomes_test <- data.table::setDT(outcomes_test)
    colnames(covariates_test) <- c("rowId", "covariateId", "covariateValue")



    ## predict on testing set
    gAtemp <- Cyclops::predict(cyclopsFit, newOutcomes = outcomes_test, newCovariates = covariates_test)


    ## store
    ID <- c(ID, idx_test)
    gA1 <- c(gA1, gAtemp)
  }

  ## result
  out <- data.frame(id=ID, gA1=gA1)
  return(out)
}
