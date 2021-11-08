#' Estimate cross-fitted TMLE of survival probability at time tau
#'
#' @param covariates Design matrix in triplet format (row index, col index, and value)
#' @param crossFit For cross-fitting: random partition of subjects into XXX prediction sets of approximately the same size.
#' @param SurvHaz.estimate Model for estimating nuisance parameter: survival hazards. Options currently include logistic LASSO
#' @param CenHaz.estimate Model for estimating nuisance parameter: censoring hazards. Options currently include logistic LASSO
#' @param TreatProb.estimate Model for estimating nuisance parameter: treatment probability. Options currently include logistic LASSO
#' @param maxCohortSizeForFitting If the target or comparator cohort are larger than this number, they
#'                                 will be downsampled before fitting the propensity model. The model
#'                                 will be used to compute propensity scores for all subjects. The
#'                                 purpose of the sampling is to gain speed.
#' @param freq.time Coarsen observed time to XXX days intervals
#' @param tau Time of interest
#' @return A list S1: S(1, tau); S0: S(0, tau); std.error.diff: standard error of S1-S0

estimateTMLEprob <- function(eventTime, censorTime, treatment, covariates, covariates.names,
                             crossFit=5, SurvHaz.estimate="LASSO", CenHaz.estimate="LASSO",
                             TreatProb.estimate="LASSO", maxCohortSizeForFitting=250000,
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
  d <- dplyr::inner_join(dH, dGR, by=c("id", "t"))


  ## into the same order
  dlong <- dlong[order(dlong$id, dlong$t), ]
  d <- d[order(d$id, d$t), ]
  dGA <- dGA[order(dGA$id), ]



  ## Update h
  h1 <- d$h1
  h0 <- d$h0
  h  <- A*h1 + (1-A)*h0
  gR1 <- d$gR1
  gR0 <- d$gR0
  gA1 <- dGA$gA1
  gA0 <- 1 - gA1
  ID <- dlong$id
  A <- dlong$treatment

  ## clear workspace
  rm(list=c("dGA"))

  ## number of subjects
  n <- length(unique(ID))
  ## time points
  m <- as.numeric(dlong$t)
  ## max follow-up time
  K <- max(m)

  ind <- (m <= tau)

  crit <- TRUE
  iter <- 1

  G1 <- tapply(1 - gR1, ID, cumprod, simplify = FALSE)
  G0 <- tapply(1 - gR0, ID, cumprod, simplify = FALSE)

  Gm1 <- unlist(G1, use.names = FALSE)
  Gm0 <- unlist(G0, use.names = FALSE)

  ## clear workspace
  rm(list=c("G1", "G0"))

  while(crit && iter <= 20){

    S1 <- tapply(1 - h1, ID, cumprod, simplify = FALSE)
    S0 <- tapply(1 - h0, ID, cumprod, simplify = FALSE)

    Sm1 <- unlist(S1, use.names = FALSE)
    Sm0 <- unlist(S0, use.names = FALSE)

    ## clear workspace
    rm(list=c("S1", "S0"))

    ## clever covariate for survival hazard
    H1 <- - (ind * Sm1) / bound(Sm1 * gA1[ID] * Gm1)
    H0 <- - (ind * Sm0) / bound(Sm0 * gA0[ID] * Gm0)
    H <- A * H1 + (1-A) * H0

    ## update for survival hazard
    eps   <- coef(glm2::glm2(Lt ~ 0 + offset(qlogis(h)) + H,
                       family = binomial(), subset = It == 1, data = dlong))

    ## NA as 0 for the new values
    eps[is.na(eps)] <- 0

    ## update values
    h1 <- bound01(plogis(qlogis(h1) + eps * H1))
    h0 <- bound01(plogis(qlogis(h0) + eps * H0))
    h  <- A * h1  + (1 - A) * h0

    iter <-  iter + 1
    crit <- abs(eps) > 1e-3/n^(0.6)
  }

  ## S_tmle
  S1 <- tapply(1 - h1, ID, cumprod, simplify = FALSE)
  S0 <- tapply(1 - h0, ID, cumprod, simplify = FALSE)

  G1 <- tapply(1 - gR1, ID, cumprod, simplify = FALSE)
  G0 <- tapply(1 - gR0, ID, cumprod, simplify = FALSE)

  Sm1 <- unlist(S1, use.names = FALSE)
  Sm0 <- unlist(S0, use.names = FALSE)

  Gm1 <- unlist(G1, use.names = FALSE)
  Gm0 <- unlist(G0, use.names = FALSE)

  H1 <- - (ind * Sm1) / bound(Sm1 * gA1[ID] * Gm1)
  H0 <- - (ind * Sm0) / bound(Sm0 * gA0[ID] * Gm0)
  DT <- with(dlong, tapply(It * (A * H1 - (1 - A) * H0) * (Lt - h), ID, sum))

  DW1 <- Sm1[which(m == tau)]
  DW0 <- Sm0[which(m == tau)]
  ## S1
  theta1 <- mean(DW1)
  ## S0
  theta0 <- mean(DW0)
  ## d_S
  theta <- theta1 - theta0
  ## standard error of d_S
  D <- DT + DW1 - DW0 - theta
  sdn <- sqrt(var(D) / n)

  ## result
  out <- list(S1=theta1, S0=theta0, std.error.diff=sdn)
  return(out)
}





