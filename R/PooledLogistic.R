#' Coefficients for pooled logistic regression (via iterative reweighted least square)
#'
#' @param X_baseline Baseline variables that won't interact with time in regression,
#'                  sparse matrix of class "dgTMatrix".
#'                  Rows are ordered into decreasing survival time.
#' @param is.temporal Whether there is temporal effect, i.e. whether time t is a variable in the regression
#' @param temporal_effect Baseline variables that will interact with time in regression,
#'                        sparse matrix of class "dgTMatrix" or matrix.
#'                        Rows are ordered into decreasing survival time
#' @param is.temporal Whether there is temporal effect, i.e. whether time t is a variable in the regression
#' @param time Observed survival time. Ordered into decreasing observed survival time
#' @param eventObserved Event indicator. Ordered into decreasing observed survival time
#' @param estimate_hazard "survival" or "censoring"
#' @return A vector of coefficients in the order of: baseline covariates, time, interaction term between baseline covariates and time


coef_pooled <- function(X_baseline, is.temporal, temporal_effect,
                        time, eventObserved, estimate_hazard,
                        maxiter, threshold){

  ## check and warning for reorder



  ## Add intercept term to X_baseline_reorder and temporal_effect_reorder
  X_baseline <- cbind(rep(1, dim(X_baseline)[1]), X_baseline)
  if(is.null(temporal_effect) & !is.temporal){
    temporal_effect <- cbind(rep(0, dim(X_baseline)[1]), temporal_effect)
  }else if(is.temporal){
    temporal_effect <- cbind(rep(1, dim(X_baseline)[1]), temporal_effect)
  }

  ## subset index for each time point
  maxTime <- max(time)
  if (estimate_hazard == "survival"){
    indx_subset <- sapply(1:maxTime, function(x) sum(time >= x), USE.NAMES = FALSE)
  }else if(estimate_hazard == "censoring"){
    indx_subset <- sapply(1:maxTime, function(x) sum((time > x)*eventObserved+(time >= x)*(1 - eventObserved)), USE.NAMES = FALSE)
  }

  ## initial value
  converged <- FALSE
  iter <- 1
  beta <- rep(0, length=dim(temporal_effect)[2] + dim(X_baseline)[2])

  ## outcome
  Y <- outcomeY(time=time, eventObserved=eventObserved, estimate_hazard=estimate_hazard)
  ## calculate X^T y
  design_matvec_Xy <- pooled_design_matvec(X_baseline=X_baseline,
                                           temporal_effect=temporal_effect,
                                           Y=Y,
                                           indx_subset=indx_subset, maxTime=maxTime)

  ## initial residual deviance
  dev_resid <- sum((Y - predict_pooled(coef=beta, X_baseline=X_baseline, temporal_effect=temporal_effect, maxTime=maxTime))^2)

  ## iterate until converge
  while((!converged) && iter <= maxiter){

    ## calculate iterative components
    comp <- pooled_design_iter(X_baseline=X_baseline,
                               temporal_effect=temporal_effect,
                               beta=beta, indx_subset=indx_subset, maxTime=maxTime)

    ## beta_new
    beta_new <- solve(comp$fisher_info, (comp$fisher_info %*% beta + design_matvec_Xy - comp$Xmu))

    ## new residual
    dev_resid_new <- sum((Y - predict_pooled(coef=beta_new, X_baseline=X_baseline, temporal_effect=temporal_effect, maxTime=maxTime))^2)

    ## stopping rule
    iter <-  iter + 1
    converged <- (abs(dev_resid_new-dev_resid)/abs(dev_resid) <= threshold)

    ## update value
    beta <- beta_new
    dev_resid <- dev_resid_new

    ## clear workspace
    rm(list=c("comp","beta_new"))

    print(iter)
  }

  ## result
  return(beta)
}

#' Output outcome Y in the long-format, include outcome with the at-risk group at each time point
#' @param time Observed survival time. Ordered into decreasing observed survival time
#' @param eventObserved Event indicator. Ordered into decreasing observed survival time
#' @param estimate_hazard "survival" or "censoring"
#' @return A vector

outcomeY <- function(time, eventObserved, estimate_hazard){

  ## container
  Y <- c()
  ## parameter
  maxTime <- max(time)

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
#' @param Y Outcome variable in the pooled logistic regression.
#'          Long-format, include outcome with the at-risk group at each time point
#' @param indx_subset Subset index for each time point
#' @return A vector

pooled_design_matvec <- function(X_baseline, temporal_effect, Y, indx_subset, maxTime){

  ## container
  result_Xy <- rep(0, length=dim(X_baseline)[2])
  result_temporaly <- rep(0, length=dim(temporal_effect)[2])


  ## loop over each time point
  for (i in 1:maxTime){
  ## parameter
  atRiskIndx <- indx_subset[i]
  n <- dim(X_baseline)[1]
  y <- Y[((i-1)*n+1):(i*n)]
  ## X^T y
  result_Xy <- result_Xy + t(X_baseline[1:atRiskIndx, , drop=FALSE]) %*% y[1:atRiskIndx]
  result_temporaly <- result_temporaly + i * t(temporal_effect[1:atRiskIndx, , drop=FALSE]) %*% y[1:atRiskIndx] ## could add functions to i
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
#' @param beta Current iteration of coefficient value
#' @param indx_subset Subset index for each time point
#' @return A list
pooled_design_iter <- function(X_baseline, temporal_effect, beta, indx_subset, maxTime){

  ## container
  fisher_info <- matrix(0, nrow=(dim(temporal_effect)[2] + dim(X_baseline)[2]),
                        ncol=(dim(temporal_effect)[2] + dim(X_baseline)[2]))
  baselineMu <- rep(0, length=dim(X_baseline)[2])
  temporalMu <- rep(0, length=dim(temporal_effect)[2])

  ## loop over each time point
  for (i in 1:maxTime){
    ## parameters
    atRiskIndx <- indx_subset[i]
    timeIndepCoef <- beta[1:dim(X_baseline)[2]]
    timeDepCoef <- beta[(dim(X_baseline)[2]+1):length(beta)]
    baselineEffect <- X_baseline[1:atRiskIndx, , drop=FALSE] %*% timeIndepCoef
    temporalEffect <- temporal_effect[1:atRiskIndx, , drop=FALSE] %*% timeDepCoef
    ## mu
    temp_mu <- 1/(1 + exp(- baselineEffect - i * temporalEffect))
    ## X^T mu
    baselineMu <- baselineMu + t(X_baseline[1:atRiskIndx, , drop=FALSE]) %*% temp_mu[, 1]
    temporalMu <- temporalMu + i * t(temporal_effect[1:atRiskIndx, , drop=FALSE]) %*% temp_mu[, 1]
    ## X^T diag(D) X
    temp_X <- t(X_baseline[1:atRiskIndx, , drop=FALSE]) %*% ((temp_mu[, 1]) * X_baseline[1:atRiskIndx, , drop=FALSE])
    temp_temporal <- (i^2) * t(temporal_effect[1:atRiskIndx, , drop=FALSE]) %*% ((temp_mu[, 1]) * temporal_effect[1:atRiskIndx, , drop=FALSE])
    temp_Xtemporal <- i * t(X_baseline[1:atRiskIndx, , drop=FALSE]) %*% ((temp_mu[, 1]) * temporal_effect[1:atRiskIndx, , drop=FALSE])
    fisher_info <- fisher_info + cbind(rbind(temp_X, t(temp_Xtemporal)), rbind(temp_Xtemporal, temp_temporal))

    rm(list=c("temp_mu", "temp_X", "temp_temporal", "temp_Xtemporal",
              "atRiskIndx", "timeIndepCoef", "timeDepCoef",
              "baselineEffect", "temporalEffect"))
  }
  ## result
  return(list(Xmu=c(baselineMu[, 1], temporalMu[, 1]),
              fisher_info=fisher_info))
}




#' Prediction for pooled logistic regression
#'
#' @param X_baseline Baseline variables that won't interact with time in regression,
#'                  sparse matrix of class "dgTMatrix".
#'                  Need to include intercept
#' @param is.temporal Whether there is temporal effect, i.e. whether time t is a variable in the regression
#' @param temporal_effect Baseline variables that will interact with time in regression,
#'                        sparse matrix of class "dgTMatrix" or matrix.
#'                        Need to include intercept if is.temporal = TRUE
#' @return Probability, in the order of: id1(t1, t2....), id2(t1, t2....)....

predict_pooled <- function(coef, X_baseline, is.temporal, temporal_effect, maxTime){

  ## predict
  timeIndepLP <- X_baseline %*% coef[1:dim(X_baseline)[2]]
  timeDepenLP <- temporal_effect %*% coef[(dim(X_baseline)[2] + 1):length(coef)]
  logitProb <- rep(timeIndepLP, each=maxTime) + rep(1:maxTime, dim(X_baseline)[1]) * rep(timeDepenLP, each=maxTime)
  predictedProb <- exp(logitProb) / (1 + exp(logitProb))

  ## result
  return(predictedProb)
}














