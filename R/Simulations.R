#' Simulate time-to-event data with survival and censoring hazards
#' @param survHaz Output from estimateHaz
#' @param cenHaz Output from estimateHaz
#' @param nInt number of time intervals for coarsening the data
#' @export

simulateData <- function(treatment, survHaz, cenHaz, coarsenedTime, seed){

  n <- length(treatment)
  nInt <- length(coarsenedTime$timeIntMidPoint)

  survHaz_all <- survHaz$Haz1 * treatment[survHaz$ID] + survHaz$Haz0 * (1 - treatment[survHaz$ID])
  cenHaz_all <- cenHaz$Haz1 * treatment[cenHaz$ID] + cenHaz$Haz0 * (1 - treatment[cenHaz$ID])

  set.seed(seed)
  rS <- rbinom(length(survHaz_all), 1, survHaz_all)
  rG <- rbinom(length(cenHaz_all), 1, cenHaz_all)

  tS <- tapply(rS, rep(1:n, each=nInt), function(x){which(x == 1)[1]})
  tG <- tapply(rG, rep(1:n, each=nInt), function(x){which(x == 1)[1]})
  tG <- ifelse(is.na(tG), nInt+1, tG)

  ObservedTime <- apply(cbind(tS, tG), 1, function(x){min(x, na.rm = TRUE)})
  ObservedEvent <- 1*(tS<=tG)
  ObservedEvent <- ifelse(is.na(ObservedEvent), 0, ObservedEvent)

  ObservedEvent <- as.double(ObservedEvent)
  ObservedTime <- as.double(as.character(factor(ObservedTime, labels = c(coarsenedTime$timeIntMidPoint, max(coarsenedTime$timeInt)))))

  return(list(ObservedEvent=ObservedEvent,
              ObservedTime=ObservedTime))

}



#' Calculate counterfactuals from simulated data
#' @param coarsenedTime Output from coarsenData
#' @param survHaz Output from estimateHaz
#' @export

calculateSurvCurve <- function(coarsenedTime, survHaz){

  n <- length(coarsenedTime$timeInt)
  nInt <- length(coarsenedTime$timeIntMidPoint)

  Sm0 <- unlist(tapply(1-survHaz$Haz0, rep(1:n, each=nInt), cumprod, simplify = FALSE), use.names = FALSE)
  S0 <- tapply(Sm0, rep(1:nInt, n), mean)
  rm(list=c("Sm0"))
  Sm1 <- unlist(tapply(1-survHaz$Haz1, rep(1:n, each=nInt), cumprod, simplify = FALSE), use.names = FALSE)
  S1 <- tapply(Sm1, rep(1:nInt, n), mean)
  rm(list=c("Sm1"))

  return(data.frame(S0=S0, S1=S1))

}

#' Calculate counterfactuals from simulated data
#' @param coarsenedTime Output from coarsenData
#' @param survCurve Output from calculateSurvCurve
#' @export

calculateRMST <- function(coarsenedTime, survCurve){

  rmst0 <- cumsum(coarsenedTime$timeIntLength * survCurve$S0)
  rmst1 <- cumsum(coarsenedTime$timeIntLength * survCurve$S1)

  return(data.frame(RMST0=rmst0, RMST1=rmst1))

}



#' Simulate censoring time, nonproportional
#'

simCenTime <- function(treatment, covariates, lambda1, nu1, lambda0, nu0, seed){

  n <- length(treatment)
  nvar <- length(unique(covariates$covariateId))
  cov <- Matrix::sparseMatrix(i = covariates$rowId, j = covariates$covariateId, x = covariates$covariateValue, repr = "T")
  time <- rep(0, n)

  set.seed(seed)
  u <- runif(n)
  beta <- log(runif(nvar, min=1, max=3))

  time[treatment == 0] <- (- log(u)[which(treatment == 0)]/(lambda0 * exp(cov[which(treatment == 0), ] %*% beta)[, 1]))^(1/nu0)
  time[treatment == 1] <- (- log(u)[which(treatment == 1)]/(lambda1 * exp(log(2) + cov[which(treatment == 1), ] %*% beta)[, 1]))^(1/nu1)

  return(time)

}

























