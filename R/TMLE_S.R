#' Estimate cross-fitted TMLE of survival probability
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

estimateTMLEprob <- function(eventTime, censorTime, treatment, covariates, covariates.names,
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

  ## number of subjects
  n <- length(unique(ID))
  ## time points
  m <- as.numeric(dlong$t)
  ## max follow-up time
  K <- max(m)

  ind <- outer(m, 1:K, "<=")

  crit <- TRUE
  iter <- 1

  while(crit && iter <= 20){

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

    ## clever covariate for survival hazard
    H1 <- - (ind * St1)[, tau] / bound(Sm1 * gA1[ID] * Gm1)
    H0 <- - (ind * St0)[, tau] / bound(Sm0 * gA0[ID] * Gm0)
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

  St1 <- do.call('rbind', S1[ID])
  St0 <- do.call('rbind', S0[ID])

  Sm1 <- unlist(S1)
  Sm0 <- unlist(S0)
  Gm1 <- unlist(G1)
  Gm0 <- unlist(G0)

  H1 <- - (ind * St1)[, tau] / bound(Sm1 * gA1[ID] * Gm1)
  H0 <- - (ind * St0)[, tau] / bound(Sm0 * gA0[ID] * Gm0)
  DT <- with(dlong, tapply(Im * (A * H1 - (1 - A) * H0) * (Lm - h), ID, sum))

  DW1 <- with(dlong, St1[t == 1, tau])
  DW0 <- with(dlong, St0[t == 1, tau])
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



#' Estimate cross-fitted nuisance parameter survival hazards
#'
#' @param dlong Long-format survival data from function transformData(dwide, freq.time)
#' @param J For cross-fitting
#' @param h.estimate Model for estimating nuisance parameter: survival hazards
#' @export

estimateNuisanceH <- function(dlong, J, h.estimate="glm"){

  ## container
  ID <- Time <- h1 <- h0 <- c()


  for (i in 1:J){
    ## training and testing sets
    idx_test <- index_ls[[i]]
    idx_train <- setdiff(unique(dlong$id), idx_test)
    d_test <- subset(dlong, id %in% idx_test)
    d_train <- subset(dlong, id %in% idx_train)


    ## Survival hazard: glm
    fitL <- glm(Lt ~ treatment * (t + age65 + cardiovascular + female + CHADS2),
                data = d_train, subset = It == 1, family = binomial())


    ## predict
    h1temp <- bound01(predict(fitL, newdata = mutate(d_test, treatment = 1), type = 'response'))
    h0temp <- bound01(predict(fitL, newdata = mutate(d_test, treatment = 0), type = 'response'))


    ## store
    ID <- c(ID, idx_test)
    Time <- c(Time, d_test$t)
    h1 <- c(h1, h1temp)
    h0 <- c(h0, h0temp)
  }
    ## result
    out <- data.frame(id=ID, t=Time, h1=h1, h0=h0)
    return(out)
}


#' Estimate cross-fitted nuisance parameter: censoring hazards
#'
#' @param dlong Long-format survival data from function transformData(dwide, freq.time)
#' @param J For cross-fitting
#' @param gR.estimate Model for estimating nuisance parameter: censoring hazards
#' @export

estimateNuisanceGR <- function(dlong, J, gR.estimate="glm"){

  ## container
  ID <- Time <- gR1 <- gR0 <- c()


  for (i in 1:J){
    ## training and testing sets
    idx_test <- index_ls[[i]]
    idx_train <- setdiff(unique(dlong$id), idx_test)
    d_test <- subset(dlong, id %in% idx_test)
    d_train <- subset(dlong, id %in% idx_train)


    ## model: glm
    fitR <- glm(Rt ~ treatment * (t + age65 + cardiovascular + female + CHADS2),
                data = d_train, subset = Jt == 1, family = binomial())


    ## predict
    gR1temp <- bound01(predict(fitR, newdata = mutate(d_test, treatment = 1), type = 'response'))
    gR0temp <- bound01(predict(fitR, newdata = mutate(d_test, treatment = 0), type = 'response'))


    ## store
    ID <- c(ID, idx_test)
    Time <- c(Time, d_test$t)
    gR1 <- c(gR1, gR1temp)
    gR0 <- c(gR0, gR0temp)
  }
    ## result
    out <- data.frame(id=ID, t=Time, gR1=gR1, gR0=gR0)
    return(out)
}



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

crossFit <- function(eventObserved, id, J){

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








