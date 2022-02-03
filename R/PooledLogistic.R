#' Coefficients for pooled logistic regression (via iterative reweighted least square)
#'
#' @param X_baseline Baseline variables that won't interact with time in regression,
#'                  sparse matrix of class "dgTMatrix".
#'                  Rows are ordered into decreasing survival time.
#' @param is.temporal Whether there is temporal effect, i.e. whether time t is a variable in the regression
#' @param temporal_effect Baseline variables that will interact with time in regression,
#'                        sparse matrix of class "dgTMatrix" or matrix.
#'                        Rows are ordered into decreasing survival time
#' @param timeEffect Functions of time in the discrete censoring hazards model.
#'                   Options currently include "linear", "ns2", "ns3", "ns4", "ns5", NULL (if is.temporal = FALSE)
#' @param time Observed survival time. Ordered into decreasing observed survival time
#' @param eventObserved Event indicator. Ordered into decreasing observed survival time
#' @param estimate_hazard "survival" or "censoring"
#' @param sigma Penalized parameter for ridge regression, a scalar. If sigma = NULL, then no penalization.
#' @param maxiter Maximum iterations
#' @param threshold Threshold for convergence
#' @param printIter TRUE/FALSE. Whether to print iterations or not
#' @return A vector of coefficients in the order of: intercept, baseline covariates, time,
#'                                                   interaction term between baseline covariates and time.
#'         And the standard error of the estimated coefficients


coef_pooled <- function(X_baseline, is.temporal, temporal_effect, timeEffect,
                        time, eventObserved, estimate_hazard, sigma,
                        maxiter, threshold, printIter, initial_coef){

  ## check and warning for reorder


  ## lambda
  if(is.null(sigma)){
    lambda <- 0
  }else{
    lambda <- 1/(2*sigma^2)
  }

  ## subset index for each time point
  if (estimate_hazard == "survival"){
    maxTime <- max(time[eventObserved==1])
    indx_subset <- sapply(1:maxTime, function(x) sum(time >= x), USE.NAMES = FALSE)
  }else if(estimate_hazard == "censoring"){
    maxTime <- max(time[eventObserved==0])
    indx_subset <- sapply(1:maxTime, function(x) sum((time > x)*eventObserved+(time >= x)*(1 - eventObserved)), USE.NAMES = FALSE)
  }


  ## Add intercept term to X_baseline and temporal_effect
  X_baseline <- cbind(rep(1, dim(X_baseline)[1]), X_baseline)
  X_baseline <- Matrix::sparseMatrix(i = Matrix::summary(X_baseline)$i, j = Matrix::summary(X_baseline)$j,
                             x = Matrix::summary(X_baseline)$x, repr = "R")
  if(is.null(temporal_effect) & !is.temporal){
    temporal_effect <- cbind(rep(0, dim(X_baseline)[1]), temporal_effect)
  }else if(is.temporal & timeEffect == "linear"){
    temporal_effect <- cbind(rep(1, dim(X_baseline)[1]), temporal_effect)
  }else if(timeEffect == "ns2"){
    temporal_effect <- cbind(rep(1, dim(X_baseline)[1]), rep(1, dim(X_baseline)[1]),
                             temporal_effect, temporal_effect)
  }else if(timeEffect == "ns3"){
    temporal_effect <- cbind(rep(1, dim(X_baseline)[1]), rep(1, dim(X_baseline)[1]), rep(1, dim(X_baseline)[1]),
                             temporal_effect, temporal_effect, temporal_effect)
  }else if(timeEffect == "ns4"){
    temporal_effect <- cbind(rep(1, dim(X_baseline)[1]), rep(1, dim(X_baseline)[1]), rep(1, dim(X_baseline)[1]),
                             rep(1, dim(X_baseline)[1]), temporal_effect, temporal_effect, temporal_effect, temporal_effect)
  }else if(timeEffect == "ns5"){
    temporal_effect <- cbind(rep(1, dim(X_baseline)[1]), rep(1, dim(X_baseline)[1]), rep(1, dim(X_baseline)[1]),
                             rep(1, dim(X_baseline)[1]), rep(1, dim(X_baseline)[1]), temporal_effect, temporal_effect,
                             temporal_effect, temporal_effect, temporal_effect)
  }

  ## initial value
  converged <- FALSE
  iter <- 1
  if(is.null(initial_coef)){
    beta <- rep(0, length=dim(temporal_effect)[2] + dim(X_baseline)[2])
  }else{
    beta <- initial_coef
  }

  ## outcome
  Y <- outcomeY(time=time, eventObserved=eventObserved, estimate_hazard=estimate_hazard, maxTime=maxTime)


  ## calculate X^T y
  design_matvec_Xy <- pooled_design_matvec(X_baseline=X_baseline,
                                           temporal_effect=temporal_effect,
                                           timeEffect=timeEffect,
                                           Y=Y,
                                           indx_subset=indx_subset, maxTime=maxTime)

  ## parameter
  comp <- pooled_design_iter(X_baseline=X_baseline,
                             temporal_effect=temporal_effect, Y=Y,
                             timeEffect=timeEffect,
                             beta=beta, indx_subset=indx_subset, maxTime=maxTime)
  ## parameter
  Imop <- diag(dim(comp$fisher_info)[1])
  Imop[1, 1] <- 0
  Imop[2, 2] <- 0
  Imop[(dim(X_baseline)[2]+1):length(beta), ] <- 0

  ## initial (penalized) log-likelihood
  logLikelihood <- comp$logLik - lambda*sum(beta[2:dim(X_baseline)[2]]^2)

  ## iterate until converge
  while((!converged) && iter <= maxiter){


    ## beta_new
    beta_new <- solve((comp$fisher_info + 2*lambda*Imop), (comp$fisher_info %*% beta + design_matvec_Xy - comp$Xmu)[, 1])


    ## calculate iterative components
    comp <- pooled_design_iter(X_baseline=X_baseline,
                               temporal_effect=temporal_effect, Y=Y,
                               timeEffect=timeEffect,
                               beta=beta_new, indx_subset=indx_subset, maxTime=maxTime)

    ## initial (penalized) log-likelihood
    logLikelihood_new <- comp$logLik - lambda*sum(beta[2:dim(X_baseline)[2]]^2)

    ## stopping rule
    iter <-  iter + 1
    converged <- (abs(logLikelihood_new-logLikelihood)/abs(logLikelihood_new) <= threshold)

    ## update value
    beta <- beta_new
    logLikelihood <- logLikelihood_new

    rm(list=c("beta_new", "logLikelihood_new"))

    if(printIter){print(iter - 1)}
  }

  ## result
  return(list(estimates=beta,
              fisherInfo = comp$fisher_info + 2*lambda*Imop,
              logLik=comp$logLik))
}


#' Output outcome Y in the long-format, include outcome with all individuals at each time point
#' @param time Observed survival time. Ordered into decreasing observed survival time
#' @param eventObserved Event indicator. Ordered into decreasing observed survival time
#' @param estimate_hazard "survival" or "censoring"
#' @param maxTime Maximum time for estimation
#' @return A vector

outcomeY <- function(time, eventObserved, estimate_hazard, maxTime){

  ## container
  Y <- c()

  ## loop over each time point
  for (i in 1:maxTime){
    ## create outcome variable
    if (estimate_hazard == "survival"){
      y <- eventObserved * (time == i)
    }else if(estimate_hazard == "censoring"){
      y <- (1 - eventObserved) * (time == i)
    }
    Y <- c(Y, y)
  }

  ## result
  return(Y)
}


#' Calculate X^T y for pooled logistic regression (this part doesn't need to iterate)
#'
#' @param X_baseline Baseline variables that won't interact with time in regression,
#'                  sparse matrix of class "dgTMatrix",
#'                  reorder observations into decreasing survival time,
#'                  Intercept included.
#' @param temporal_effect Baseline variables that will interact with time in regression,
#'                        sparse matrix of class "dgTMatrix" or matrix,
#'                        reorder observations into decreasing survival time,
#'                        intercept included if there is temporal effect
#' @param timeEffect Functions of time in the discrete censoring hazards model.
#'                   Options currently include "linear", "cubic", NULL
#' @param Y Outcome variable in the pooled logistic regression.
#'          Long-format, include outcome with all individuals at each time point
#' @param indx_subset Subset index for each time point
#' @param maxTime Maximum time for estimation
#' @return A vector

pooled_design_matvec <- function(X_baseline, temporal_effect, timeEffect, Y, indx_subset, maxTime){

  ## container
  result_Xy <- rep(0, length=dim(X_baseline)[2])
  result_temporaly <- rep(0, length=dim(temporal_effect)[2])

  ## natural spline for time
  if(timeEffect %in% c("ns2", "ns3", "ns4", "ns5")){
    indx <- as.integer(substr(timeEffect, 3, 3))
    nsBase <- splines::ns(c(1:maxTime), df=indx)
  }

  ## loop over each time point
  for (i in 1:maxTime){
    ## parameter
    atRiskIndx <- indx_subset[i]
    n <- dim(X_baseline)[1]
    y <- Y[((i-1)*n+1):(i*n)]
    ## X^T y
    result_Xy <- result_Xy + computeSubsetSparseMatVec(X=X_baseline, v=y, subsetSize=atRiskIndx, transposed=TRUE)
    if(timeEffect == "linear" | is.null(timeEffect)){
      result_temporaly <- result_temporaly + computeSubsetMatVec(Y=temporal_effect, v=y, subsetSize=atRiskIndx, transposed=TRUE) * i
    }else if(timeEffect %in% c("ns2", "ns3", "ns4", "ns5")){
      timeNS <- c(nsBase[i, ], rep(nsBase[i, ], each=(dim(temporal_effect)[2]-indx)/indx))
      result_temporaly <- result_temporaly + computeSubsetMatVec(Y=temporal_effect, v=y, subsetSize=atRiskIndx, transposed=TRUE) * timeNS
    }
  }
  ## result
  return(c(result_Xy, result_temporaly))
}


#' Calculate X^T mu, D and X^T diag(D) X for pooled logistic regression
#'
#' @param X_baseline Baseline variables that won't interact with time in regression,
#'                  sparse matrix of class "dgTMatrix",
#'                  reorder observations into decreasing survival time,
#'                  intercept included
#' @param temporal_effect Baseline variables that will interact with time in regression,
#'                        sparse matrix of class "dgTMatrix" or matrix,
#'                        reorder observations into decreasing survival time,
#'                        intercept included
#' @param Y Outcome variable in the pooled logistic regression.
#'          Long-format, include outcome with all individuals at each time point
#' @param timeEffect Functions of time in the discrete censoring hazards model.
#'                   Options currently include "linear", "ns", NULL
#' @param beta Current iteration of coefficient value
#' @param indx_subset Subset index for each time point
#' @param maxTime Maximum time for estimation
#' @return A list

pooled_design_iter <- function(X_baseline, temporal_effect, Y, timeEffect, beta, indx_subset, maxTime){

  ## container
  fisher_info <- matrix(0, nrow=(dim(temporal_effect)[2] + dim(X_baseline)[2]),
                        ncol=(dim(temporal_effect)[2] + dim(X_baseline)[2]))
  baselineMu <- rep(0, length=dim(X_baseline)[2])
  temporalMu <- rep(0, length=dim(temporal_effect)[2])
  logLik <- 0

  ## natural spline for time
  if(timeEffect %in% c("ns2", "ns3", "ns4", "ns5")){
    indx <- as.integer(substr(timeEffect, 3, 3))
    nsBase <- splines::ns(c(1:maxTime), df=indx)
  }

  ## loop over each time point
  for (i in 1:maxTime){
    ## parameters
    atRiskIndx <- indx_subset[i]
    timeIndepCoef <- beta[1:dim(X_baseline)[2]]
    timeDepCoef <- beta[(dim(X_baseline)[2]+1):length(beta)]
    baselineEffect <- computeSubsetSparseMatVec(X=X_baseline, v=timeIndepCoef, subsetSize=atRiskIndx, transposed=FALSE)
    n <- dim(X_baseline)[1]
    y <- Y[((i-1)*n+1):(i*n)]
    ## mu
    if((timeEffect == "linear" | is.null(timeEffect))){
      temporalEffect <- computeSubsetMatVec(Y=temporal_effect, v=timeDepCoef, subsetSize=atRiskIndx, transposed=FALSE) * i
      temp_mu <- 1/(1 + exp(- baselineEffect - temporalEffect))
    }else if(timeEffect %in% c("ns2", "ns3", "ns4", "ns5")){
      timeNS <- c(nsBase[i, ], rep(nsBase[i, ], each=(dim(temporal_effect)[2]-indx)/indx))
      temporalEffect <- computeSubsetMatVec(Y=temporal_effect, v=(timeDepCoef*timeNS), subsetSize=atRiskIndx, transposed=FALSE)
      temp_mu <- 1/(1 + exp(- baselineEffect - temporalEffect))
    }

    ## mu(1-mu)
    diagmu <- temp_mu*(1-temp_mu)
    ## X^T mu
    baselineMu <- baselineMu + computeSubsetSparseMatVec(X=X_baseline, v=temp_mu, subsetSize=atRiskIndx, transposed=TRUE)
    if(timeEffect == "linear" | is.null(timeEffect)){
      temporalMu <- temporalMu + computeSubsetMatVec(Y=temporal_effect, v=temp_mu, subsetSize=atRiskIndx, transposed=TRUE) * i
    }else if(timeEffect %in% c("ns2", "ns3", "ns4", "ns5")){
      temporalMu <- temporalMu + computeSubsetMatVec(Y=temporal_effect, v=temp_mu, subsetSize=atRiskIndx, transposed=TRUE) * timeNS
    }
    ## X^T diag(D) X
    temp_X <- computeSubsetSparseInformationMatrix(X=X_baseline, weight=diagmu, subsetSize=atRiskIndx)
    if(timeEffect == "linear" | is.null(timeEffect)){
      temp_temporal <- computeSubsetInformationMatrix(Y=temporal_effect, weight=diagmu, subsetSize=atRiskIndx) * i^2
      temp_Xtemporal <- computeSubsetMatrix(X=X_baseline, weight=diagmu, Y=temporal_effect, subsetSize=atRiskIndx) * i
    }else if(timeEffect %in% c("ns2", "ns3", "ns4", "ns5")){
      temp_temporal <- computeSubsetInformationMatrix(Y=t(t(temporal_effect) * timeNS), weight=diagmu, subsetSize=atRiskIndx)
      temp_Xtemporal <- computeSubsetMatrix(X=X_baseline, weight=diagmu, Y=t(t(temporal_effect)*timeNS), subsetSize=atRiskIndx)
    }
    fisher_info <- fisher_info + cbind(rbind(temp_X, t(temp_Xtemporal)), rbind(temp_Xtemporal, temp_temporal))
    ## log likelihood
    logLik_temp <- y[1:atRiskIndx] * log(temp_mu) + (1-y[1:atRiskIndx]) * log(1-temp_mu)
    logLik_temp[is.infinite(logLik_temp)] <- NA
    logLik <- logLik + sum(logLik_temp, na.rm = TRUE)

    rm(list=c("temp_mu", "temp_X", "temp_temporal", "temp_Xtemporal",
              "atRiskIndx", "timeIndepCoef", "timeDepCoef",
              "baselineEffect", "temporalEffect", "diagmu", "logLik_temp"))
  }
  ## result
  return(list(Xmu=unname(c(baselineMu, temporalMu)),
              fisher_info=unname(fisher_info), logLik=logLik))
}




#' Prediction for pooled logistic regression
#'
#' @param X_baseline Baseline variables that won't interact with time in regression,
#'                  sparse matrix of class "dgTMatrix".
#'                  Need to include intercept
#' @param temporal_effect Baseline variables that will interact with time in regression,
#'                        sparse matrix of class "dgTMatrix" or matrix.
#'                        Need to include intercept if is.temporal = TRUE
#' @param timeEffect Functions of time in the discrete censoring hazards model.
#'                   Options currently include "linear", "cubic", NULL
#' @param indx_subset Subset index for each time point
#' @param maxTime Maximum time for prediction
#' @param maxTimeSplines Maximum time used in regression. Relavant when using ns(time, df=3)
#' @return Probability, in the order of: id1(t1, t2....), id2(t1, t2....)....

predict_pooled <- function(coef, X_baseline, temporal_effect, timeEffect, maxTime, maxTimeSplines){

  ## predict
  timeIndepLP <- computeSubsetSparseMatVec(X=X_baseline, v=coef[1:dim(X_baseline)[2]], subsetSize=dim(X_baseline)[1], transpose=FALSE)

  if(timeEffect == "linear" | is.null(timeEffect)){

    timeDepenLP <- computeSubsetMatVec(Y=temporal_effect, v=coef[1:dim(X_baseline)[2]], subsetSize=dim(X_baseline)[1], transpose=FALSE)
    logitProb <- rep(timeIndepLP, each=maxTime) + rep(1:maxTime, dim(X_baseline)[1]) * rep(timeDepenLP, each=maxTime)

  }else if(timeEffect %in% c("ns2", "ns3", "ns4", "ns5")){

    indx <- as.integer(substr(timeEffect, 3, 3))
    nsBase <- splines::ns(c(1:maxTimeSplines), df=indx)

    timeDepenLP <- c()
    for (i in 1:maxTime){
      timeDepenLP_temp <- computeSubsetMatVec(Y=temporal_effect,
                                              v=(coef[(dim(X_baseline)[2] + 1):length(coef)] * c(nsBase[i, ], rep(nsBase[i, ], each=(dim(temporal_effect)[2]-indx)/indx))),
                                              subsetSize=dim(X_baseline)[1], transpose=FALSE)
      timeDepenLP <- c(timeDepenLP, timeDepenLP_temp)
      rm(timeDepenLP_temp)
    }
    timeDepenLP <- timeDepenLP[order(rep(1:dim(X_baseline)[1], maxTime))]
    logitProb <- rep(timeIndepLP, each=maxTime) + timeDepenLP

  }

  predictedProb <- exp(logitProb) / (1 + exp(logitProb))

  ## result
  return(predictedProb)
}




#' “evidence maximization” approach for tunning the penalty parameter
#' @param X_baseline Baseline variables that won't interact with time in regression,
#'                  sparse matrix of class "dgTMatrix".
#'                  Rows are ordered into decreasing survival time.
#' @param is.temporal Whether there is temporal effect, i.e. whether time t is a variable in the regression
#' @param temporal_effect Baseline variables that will interact with time in regression,
#'                        sparse matrix of class "dgTMatrix" or matrix.
#'                        Rows are ordered into decreasing survival time
#' @param timeEffect Functions of time in the discrete censoring hazards model.
#'                   Options currently include "linear", "ns", NULL (if is.temporal = FALSE)
#' @param time Observed survival time. Ordered into decreasing observed survival time
#' @param eventObserved Event indicator. Ordered into decreasing observed survival time
#' @param estimate_hazard "survival" or "censoring"
#' @param sigma A range of penalized parameters for ridge regression, a vector.
#' @param maxiter Maximum iterations
#' @param threshold Threshold for convergence
#' @param printIter TRUE/FALSE. Whether to print iterations or not

coef_ridge <- function(X_baseline, is.temporal, temporal_effect,
                      timeEffect, eventObserved, time,
                      estimate_hazard, sigma, maxiter, threshold, printIter){

  r_coef <- c()
  r_lik <- c()

  for (i in 1:length(sigma)){

    if(i == 1){
      initial_coef <- NULL
    }else{initial_coef <- coef_temp$estimates}

    try({
    ## coef's for each sigma
    coef_temp <- coef_pooled(X_baseline=X_baseline, is.temporal=is.temporal, temporal_effect=temporal_effect,
                             timeEffect=timeEffect, eventObserved=eventObserved, time=time,
                             estimate_hazard=estimate_hazard, sigma=sigma[i],
                             maxiter=maxiter, threshold=threshold, printIter=printIter, initial_coef=initial_coef)

    ## marginal likelihood for each sigma
    marginalLogLik_temp <- marginalLogLik(betaMAP=coef_temp$estimates, sigma=sigma[i], p=dim(X_baseline)[2],
                                    logLik=coef_temp$logLik, fisherInfo=coef_temp$fisherInfo)

    ## result
    r_coef <- cbind(r_coef, coef_temp$estimates)
    r_lik <- c(r_lik, marginalLogLik_temp)

    })

    if(length(r_lik) > 3){
      crit <- (r_lik[length(r_lik)] < r_lik[length(r_lik)-1]) & (r_lik[length(r_lik)-1] < r_lik[length(r_lik)-2])
      if(crit){
        break
      }
    }
     }

  ## pick the coef's and sigma with the largest marginal likelihood
  pick <- which.max(r_lik)
  return(list(estimates=unname(r_coef[, pick]), all_like=r_lik))
}




#' Laplace’s method for estimating marginal likelihood

marginalLogLik <- function(betaMAP, sigma, p, logLik, fisherInfo){

  ## marginal log likelihood
  logLikelihood <- logLik - p * log(sqrt(2*pi)*sigma) - sum(betaMAP[2:(p+1)]^2)/(2*sigma^2) + p * log(sqrt(2*pi)) + determinant(solve(as.matrix(fisherInfo)))$modulus/2

  ## result
  return(logLikelihood)

}









