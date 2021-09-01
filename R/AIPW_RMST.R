#' Estimate cross-fitted Augmented IPW of RMST
#'
#' @param covariants Design matrix in triplet format (row index, col index, and value)
#' @param J For cross-fitting: random partition of subjects into J prediction sets of approximately the same size.
#' @param h.estimate Model for estimating nuisance parameter: survival hazards. Options currently include glm
#' @param gR.estimate Model for estimating nuisance parameter: censoring hazards. Options currently include glm
#' @param gA.estimate Model for estimating nuisance parameter: treatment probability. Options currently include LASSO
#' @param maxCohortSizeForFitting If the target or comparator cohort are larger than this number, they
#'                                 will be downsampled before fitting the propensity model. The model
#'                                 will be used to compute propensity scores for all subjects. The
#'                                 purpose of the sampling is to gain speed.
#' @param tau Time of interest
#' @export

estimateTMLErmst <- function(eventTime, censorTime, treatment, covariates, covariates.names,
                             J=5, h.estimate="glm", gR.estimate="glm",
                             gA.estimate="LASSO", maxCohortSizeForFitting=250000,
                             #                           includeCovariatesInHestimate="",includeCovariatesIngRestimate="", includeCovariatesIngAestimate="all",
                             freq.time=90, tau){


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



  ## Temporary: add covariates to dlong for estimating h and gR
  cov_name <- c("gender = FEMALE", "CHADS2", "Subgroup: Elderly (age >=65)", "condition_era group during day -30 through 0 days relative to index: Acute disease of cardiovascular system")
  index_cov <- which(covariates.names %in% cov_name)
  cov <- Matrix::sparseMatrix(i = covariates$i, j = covariates$j, x = covariates$val, repr = "T")
  cov_new <- cov[, index_cov]
  dwide <- cbind(dwide, as.matrix(cov_new))
  colnames(dwide)[5:8] <- c("age65", "cardiovascular", "female", "CHADS2")


  ## transform data into long format
  dlong <- transformData(dwide=dwide, freq.time=freq.time)



  ## Estimate cross-fitted nuisance parameter
  dH <- estimateNuisanceH(dlong, J=J, h.estimate="glm")
  dGR <- estimateNuisanceGR(dlong, J=J, gR.estimate="glm")
  dGA <- estimateNuisanceGA(J=J, id=id, treatment=treatment, covariates=covariates, gA.estimate="LASSO", maxCohortSizeForFitting=maxCohortSizeForFitting)
  d <- cbind(dH, dGR)



  ## Update h
  h1 <- d$h1
  h0 <- d$h0
  gR1 <- d$gR1
  gR0 <- d$gR0
  gA1 <- dGA$gA1
  ID <- dlong$id
  A <- dlong$treatment

  gA0 <- 1 - gA1
  h  <- A * h1 + (1-A) * h0
  gR <- A * gR1 + (1 - A) * gR0

  ## number of subjects
  n <- length(unique(ID))
  ## time points
  m <- as.numeric(dlong$t)
  ## max follow-up time
  K <- max(m)

  ind <- outer(m, 1:K, "<=")


  S1 <- tapply(1 - h1, ID, cumprod, simplify = FALSE)
  S0 <- tapply(1 - h0, ID, cumprod, simplify = FALSE)

  G1 <- tapply(1 - gR1, ID, cumprod, simplify = FALSE)
  G0 <- tapply(1 - gR0, ID, cumprod, simplify = FALSE)

  St1 <- do.call("rbind", S1[ID])
  St0 <- do.call("rbind", S0[ID])

  Sm1 <- unlist(S1)
  Sm0 <- unlist(S0)

  Gm1 <- unlist(G1)
  Gm0 <- unlist(G0)


  H1 <- - rowSums((ind * St1)[, 1:(tau-1)]) / bound(Sm1 * gA1[ID] * Gm1)
  H0 <- - rowSums((ind * St0)[, 1:(tau-1)]) / bound(Sm0 * gA0[ID] * Gm0)

  DT1 <- with(dlong, tapply(Im * A * H1 * (Lm - h), ID, sum))
  DT0 <- with(dlong, tapply(Im * (1 - A) * H0 * (Lm - h), ID, sum))

  DW1 <- with(dlong, rowSums(St1[t == 1, 1:(tau-1)]))
  DW0 <- with(dlong, rowSums(St0[t == 1, 1:(tau-1)]))

  ## AIPW
  aipw <- 1+c(mean(DT0 + DW0), mean(DT1 + DW1))
  ## standard error
  D <- DT1 - DT0 + DW1 - DW0
  sdn <- sqrt(var(D) / n)

  ## result
  out <- list(rmst1=aipw[2], rmst0=aipw[1], std.error.diff=sdn)
  return(out)
}














