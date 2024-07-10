#' Estimate nuisance parameter: treatment probability
#'
#' @param gA.estimate Model for estimating nuisance parameter: treatment probability
#' @param maxCohortSizeForFitting
#' @export

estimateNuisanceGA2 <- function(id, treatment, covariates, gA.estimate="LASSO", maxCohortSizeForFitting){

  ## outcomes
  outcomes <- data.frame(rowId = id, y = treatment)


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
  prior = createPrior("laplace", exclude = c(0), useCrossValidation = TRUE)
  control = createControl(noiseLevel = "silent", cvType = "auto", seed = 1, tolerance = 2e-07, cvRepetitions = 10, startingVariance = 0.01)


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


  ## predict on testing set
  gA1 <- Cyclops::predict(cyclopsFit, newOutcomes = outcomes, newCovariates = covariates)


  ## result
  out <- data.frame(id=id, gA1=gA1)
  return(out)
}


