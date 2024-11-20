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
#'                   Options currently include "linear", "ns", NULL (if is.temporal = FALSE)
#' @param time Observed survival time. Ordered into decreasing observed survival time
#' @param eventObserved Event indicator. Ordered into decreasing observed survival time
#' @param estimate_hazard "survival" or "censoring"
#' @param lambda Penalized parameter for ridge regression, a scalar. If lambda = NULL, then no penalization.
#' @param maxiter Maximum iterations
#' @param threshold Threshold for convergence
#' @param printIter TRUE/FALSE. Whether to print iterations or not
#' @return A vector of coefficients in the order of: intercept, baseline covariates, time,
#'                                                   interaction term between baseline covariates and time.
#'         And the standard error of the estimated coefficients


coef_pooled <- function(X_baseline, is.temporal, temporal_effect, timeEffect,
                        time, eventObserved, estimate_hazard, lambda,
                        maxiter, threshold, printIter){

  ## check and warning for reorder



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
  if(is.null(temporal_effect) & !is.temporal){
    temporal_effect <- cbind(rep(0, dim(X_baseline)[1]), temporal_effect)
  }else if(is.temporal & timeEffect == "linear"){
    temporal_effect <- cbind(rep(1, dim(X_baseline)[1]), temporal_effect)
  }else if(timeEffect == "ns"){
    temporal_effect <- cbind(rep(1, dim(X_baseline)[1]), rep(1, dim(X_baseline)[1]), rep(1, dim(X_baseline)[1]),
                             rep(1, dim(X_baseline)[1]), rep(1, dim(X_baseline)[1]), temporal_effect,
                             temporal_effect, temporal_effect, temporal_effect, temporal_effect)
  }

  ## initial value
  converged <- FALSE
  iter <- 1
  beta <- rep(0, length=dim(temporal_effect)[2] + dim(X_baseline)[2])

  ## outcome
  Y <- outcomeY(time=time, eventObserved=eventObserved, estimate_hazard=estimate_hazard, maxTime=maxTime)


  ## calculate X^T y
  design_matvec_Xy <- pooled_design_matvec(X_baseline=X_baseline,
                                           temporal_effect=temporal_effect,
                                           timeEffect=timeEffect,
                                           Y=Y,
                                           indx_subset=indx_subset, maxTime=maxTime)

  comp <- pooled_design_iter(X_baseline=X_baseline,
                             temporal_effect=temporal_effect, Y=Y,
                             timeEffect=timeEffect,
                             beta=beta, indx_subset=indx_subset, maxTime=maxTime)

  ## initial residual deviance
  # dev_resid <- resid_pooled(coef=beta, X_baseline=X_baseline,
  #                           temporal_effect=temporal_effect,
  #                           timeEffect=timeEffect, Y=Y,
  #                           indx_subset=indx_subset, maxTime=maxTime,
  #                           lambda = lambda)

  ## initial (penalized) log-likelihood
  logLikelihood <- comp$logLik - lambda*sum(beta[2:dim(X_baseline)[2]]^2)

  ## iterate until converge
  while((!converged) && iter <= maxiter){


    ## beta_new
    Imop <- diag(dim(comp$fisher_info)[1])
    Imop[1, 1] <- 0
    Imop[(dim(X_baseline)[2]+1):length(beta), ] <- 0
    beta_new <- solve((comp$fisher_info + 2*lambda*Imop), (comp$fisher_info %*% beta + design_matvec_Xy - comp$Xmu)[, 1])

    ## new residual
    # dev_resid_new <- resid_pooled(coef=beta_new, X_baseline=X_baseline,
    #                               temporal_effect=temporal_effect, timeEffect=timeEffect,
    #                               Y=Y, indx_subset=indx_subset, maxTime=maxTime, lambda = lambda)

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
  r_var <- solve(comp$fisher_info + 2*lambda*Imop) %*% comp$fisher_info %*% solve(comp$fisher_info + 2*lambda*Imop)

  return(list(estimates=beta,
             sd=sqrt(diag(r_var)),
              fisherInfo = comp$fisher_info + 2*lambda*diag(dim(comp$fisher_info)[1]),
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
  if(timeEffect == "ns"){
    nsBase <- splines::ns(c(1:maxTime), df=5)
  }

  ## loop over each time point
  for (i in 1:maxTime){
    ## parameter
    atRiskIndx <- indx_subset[i]
    n <- dim(X_baseline)[1]
    y <- Y[((i-1)*n+1):(i*n)]
    ## X^T y
    result_Xy <- result_Xy + Matrix::t(X_baseline[1:atRiskIndx, ,drop=FALSE]) %*% y[1:atRiskIndx]
    if(timeEffect == "linear" | is.null(timeEffect)){
      result_temporaly <- result_temporaly + i * t(temporal_effect[1:atRiskIndx, ,drop=FALSE]) %*% y[1:atRiskIndx]
    }else if(timeEffect == "ns"){
      timeNS <- c(nsBase[i, ], rep(nsBase[i, ], each=(dim(temporal_effect)[2]-5)/5))
      result_temporaly <- result_temporaly + (t(temporal_effect[1:atRiskIndx, ,drop=FALSE]) %*% y[1:atRiskIndx]) * timeNS
    }
  }
  ## result
  return(c(result_Xy[, 1], result_temporaly[, 1]))
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
  if(timeEffect == "ns"){
    nsBase <- splines::ns(c(1:maxTime), df=5)
  }

  ## loop over each time point
  for (i in 1:maxTime){
    ## parameters
    atRiskIndx <- indx_subset[i]
    timeIndepCoef <- beta[1:dim(X_baseline)[2]]
    timeDepCoef <- beta[(dim(X_baseline)[2]+1):length(beta)]
    baselineEffect <- X_baseline[1:atRiskIndx, ,drop=FALSE] %*% timeIndepCoef
    n <- dim(X_baseline)[1]
    y <- Y[((i-1)*n+1):(i*n)]
    ## mu
    if((timeEffect == "linear" | is.null(timeEffect))){
      temporalEffect <- temporal_effect[1:atRiskIndx, ,drop=FALSE] %*% timeDepCoef * i
      temp_mu <- 1/(1 + exp(- baselineEffect - temporalEffect))
    }else if(timeEffect == "ns"){
      timeNS <- c(nsBase[i, ], rep(nsBase[i, ], each=(dim(temporal_effect)[2]-5)/5))
      temporalEffect <- temporal_effect[1:atRiskIndx, ,drop=FALSE] %*% (timeDepCoef * timeNS)
      temp_mu <- 1/(1 + exp(- baselineEffect - temporalEffect))
    }
    ## mu(1-mu)
    diagmu <- temp_mu[, 1]*(1-temp_mu[, 1])
    ## X^T mu
    baselineMu <- baselineMu + Matrix::t(X_baseline[1:atRiskIndx, ,drop=FALSE]) %*% temp_mu[, 1]
    if(timeEffect == "linear" | is.null(timeEffect)){
      temporalMu <- temporalMu + t(temporal_effect[1:atRiskIndx, ,drop=FALSE]) %*% temp_mu[, 1] * i
    }else if(timeEffect == "ns"){
      temporalMu <- temporalMu +  t(temporal_effect[1:atRiskIndx, ,drop=FALSE]) %*% temp_mu[, 1] * timeNS
    }
    ## X^T diag(D) X
    temp_X <- Matrix::t(X_baseline[1:atRiskIndx, ,drop=FALSE]) %*% (diagmu * X_baseline[1:atRiskIndx, ,drop=FALSE])
    if(timeEffect == "linear" | is.null(timeEffect)){
      temp_temporal <- (i^2 * t(temporal_effect[1:atRiskIndx, ,drop=FALSE])) %*% (diagmu * temporal_effect[1:atRiskIndx, ,drop=FALSE])
      temp_Xtemporal <- (i * Matrix::t(X_baseline[1:atRiskIndx, ,drop=FALSE])) %*% (diagmu * temporal_effect[1:atRiskIndx, ,drop=FALSE])
    }else if(timeEffect == "ns"){
      temp_temporal <- (t(temporal_effect[1:atRiskIndx, ,drop=FALSE]) * timeNS) %*% (diagmu * t(t(temporal_effect[1:atRiskIndx, ,drop=FALSE]) * timeNS))
      temp_Xtemporal <- Matrix::t(X_baseline[1:atRiskIndx, ,drop=FALSE]) %*% (diagmu * t(t(temporal_effect[1:atRiskIndx, ,drop=FALSE]) * timeNS))
    }
    fisher_info <- fisher_info + cbind(rbind(temp_X, Matrix::t(temp_Xtemporal)), rbind(temp_Xtemporal, temp_temporal))
    ## log likelihood
    logLik <- logLik + sum(y[1:atRiskIndx] * log(temp_mu[, 1])) + sum((1-y[1:atRiskIndx]) * log(1-temp_mu[, 1]))

    rm(list=c("temp_mu", "temp_X", "temp_temporal", "temp_Xtemporal",
              "atRiskIndx", "timeIndepCoef", "timeDepCoef",
              "baselineEffect", "temporalEffect", "diagmu"))
  }
  ## result
  return(list(Xmu=c(baselineMu[, 1], temporalMu[, 1]),
              fisher_info=fisher_info, logLik=logLik))
}



#' Residual deviance for pooled logistic regression
#'
#' @param X_baseline Baseline variables that won't interact with time in regression,
#'                  sparse matrix of class "dgTMatrix".
#'                  Need to include intercept
#' @param temporal_effect Baseline variables that will interact with time in regression,
#'                        sparse matrix of class "dgTMatrix" or matrix.
#'                        Need to include intercept if is.temporal = TRUE
#' @param timeEffect Functions of time in the discrete censoring hazards model.
#'                   Options currently include "linear", "cubic", NULL
#' @param Y Outcome variable in the pooled logistic regression.
#'          Long-format, include outcome with all individuals at each time point
#' @param indx_subset Subset index for each time point
#' @param maxTime Maximum time for estimation
#'

resid_pooled <- function(coef, X_baseline, temporal_effect, timeEffect, Y, indx_subset, maxTime, lambda){

  ## container
  resid <- 0

  ## natural spline for time
  if(timeEffect == "ns"){
    nsBase <- splines::ns(c(1:maxTime), df=5)
  }

  ## loop over each time point
  for (i in 1:maxTime){
    atRiskIndx <- indx_subset[i]
    n <- dim(X_baseline)[1]
    y <- Y[((i-1)*n+1):(i*n)]
    ## predict
    timeIndepLP_temp <- X_baseline[1:atRiskIndx, ,drop=FALSE] %*% coef[1:dim(X_baseline)[2]]
    if(timeEffect == "linear" | is.null(timeEffect)){
      timeDepenLP_temp <- temporal_effect[1:atRiskIndx, ,drop=FALSE] %*% coef[(dim(X_baseline)[2] + 1):length(coef)] * i
    }else if(timeEffect == "ns"){
      timeDepenLP_temp <- temporal_effect[1:atRiskIndx, ,drop=FALSE] %*% (coef[(dim(X_baseline)[2] + 1):length(coef)] * c(nsBase[i, ], rep(nsBase[i, ], each=(dim(temporal_effect)[2]-5)/5)))
    }
    logitProb <- timeIndepLP_temp + timeDepenLP_temp
    predictedProb <- exp(logitProb) / (1 + exp(logitProb))
    ## residuals
    resid_temp <- sum((y[1:atRiskIndx]-predictedProb)^2)
    ## store
    resid <- resid + resid_temp
  }

  resid <- resid + lambda * sum(coef[2:dim(X_baseline)[2]]^2)

  ## result
  return(resid)
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
#' @param maxTimeSplines Maximum time used in regression. Relavant when using ns(time, df=5)
#' @return Probability, in the order of: id1(t1, t2....), id2(t1, t2....)....

predict_pooled <- function(coef, X_baseline, temporal_effect, timeEffect, maxTime, maxTimeSplines){

  ## predict
  timeIndepLP <- X_baseline %*% coef[1:dim(X_baseline)[2]]

  if(timeEffect == "linear" | is.null(timeEffect)){
    timeDepenLP <- temporal_effect %*% coef[(dim(X_baseline)[2] + 1):length(coef)]
    logitProb <- rep(timeIndepLP, each=maxTime) + rep(1:maxTime, dim(X_baseline)[1]) * rep(timeDepenLP, each=maxTime)
  }else if(timeEffect == "ns"){
    nsBase <- splines::ns(c(1:maxTimeSplines), df=5)
    timeDepenLP <- unlist(lapply(1:maxTime, function(i) {temporal_effect %*% (coef[(dim(X_baseline)[2] + 1):length(coef)] * c(nsBase[i, ], rep(nsBase[i, ], each=(dim(temporal_effect)[2]-5)/5)))}), use.names = FALSE)
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
#' @param lambda A range of penalized parameters for ridge regression, a vector.
#' @param maxiter Maximum iterations
#' @param threshold Threshold for convergence
#' @param printIter TRUE/FALSE. Whether to print iterations or not

coef_ridge <- function(X_baseline, is.temporal, temporal_effect,
                      timeEffect, eventObserved, time,
                      estimate_hazard, lambda, maxiter, threshold, printIter){

  result_temp <- sapply(lambda, function(s){
    ## coef's for each lambda
    coef_temp <- coef_pooled(X_baseline=X_baseline, is.temporal=is.temporal, temporal_effect=temporal_effect,
                             timeEffect=timeEffect, eventObserved=eventObserved, time=time,
                             estimate_hazard=estimate_hazard, lambda=s, maxiter=maxiter, threshold=threshold, printIter=printIter)
    ## marginal likelihood for each lambda
    marginalLogLik_temp <- marginalLogLik(betaMAP=coef_temp$estimates, lambda=s, p=length(coef_temp$estimates),
                                    logLik=coef_temp$logLik, fisherInfo=coef_temp$fisherInfo)
    ## result
    return(c(coef_temp$estimates, marginalLogLik_temp))
  }, USE.NAMES = FALSE)

  ## pick the coef's and lambda with the largest marginal likelihood
  pick <- which.max(result_temp[dim(result_temp)[1],])
  return(result_temp[-dim(result_temp)[1], pick])
}




#' Laplace’s method for estimating marginal likelihood

marginalLogLik <- function(betaMAP, lambda, p, logLik, fisherInfo){

  sigma <- 1/sqrt(2*lambda)

  ## marginal log likelihood
  logLikelihood <- logLik - p * log(sqrt(2*pi)*sigma) - sum(betaMAP^2)/(2*sigma^2) + p * log(sqrt(2*pi)) + Matrix::det(fisherInfo, logarithm = TRUE)/2

  ## result
  return(logLikelihood)

}









