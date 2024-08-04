#' Coefficients for pooled logistic regression (via iterative reweighted least square)
#'
#' @param X_baseline Baseline variables that won't interact with time in regression,
#'                  sparse matrix of class "dgTMatrix"
#' @param temporal_effect Baseline variables that will interact with time in regression,
#'                        sparse matrix of class "dgTMatrix" or matrix
#' @param eventObserved Binary outcome variable
#' @param time Observed time
#' @param id Subject id
#' @param estimate_hazard "survival" or "censoring"
#' @export
coef_pooled <- function(X_baseline, temporal_effect, eventObserved, time, id, estimate_hazard){

  ## index for decreasing survival time
  indx <- order(time, decreasing = TRUE)
  ## order observed data into decreasing survival time
  id_reorder <- id[indx]
  eventObserved_reorder <- eventObserved[indx]
  time_reorder <- time[indx]
  X_baseline_reorder <- X_baseline[indx, ]
  if(is.null(temporal_effect)){
    temporal_effect_reorder <- NULL
  }else if(is.null(dim(temporal_effect))){
    temporal_effect_reorder <- temporal_effect[indx]
  }else{
    temporal_effect_reorder <- temporal_effect[indx, ]
  }
  ## clear workspace
  rm(list=c("id", "eventObserved", "time", "X_baseline", "temporal_effect"))

  ## Add intercept term to X_baseline_reorder and temporal_effect_reorder
  X_baseline_reorder <- cbind(rep(1, dim(X_baseline_reorder)[1]), X_baseline_reorder)
  if(is.null(temporal_effect_reorder)){
    temporal_effect_reorder <- cbind(rep(0, dim(X_baseline_reorder)[1]), temporal_effect_reorder)
  }else{
    temporal_effect_reorder <- cbind(rep(1, dim(X_baseline_reorder)[1]), temporal_effect_reorder)
  }

  ## subset index for each time point
  K <- max(time_reorder)
  if (estimate_hazard == "survival"){
    indx_subset <- sapply(1:K, function(x) sum(time_reorder>=x), USE.NAMES = FALSE)
  }else if(estimate_hazard == "censoring"){
    indx_subset <- sapply(1:K, function(x) sum((time_reorder>x)*eventObserved_reorder+(time_reorder>=x)*(1-eventObserved_reorder)), USE.NAMES = FALSE)
  }

  ## initial value
  crit <- TRUE
  iter <- 1
  beta <- rep(0, length=dim(temporal_effect_reorder)[2]+dim(X_baseline_reorder)[2])

  ## calculate X^T y
  design_matvec_Xy <- pooled_design_matvec(X_baseline_reorder=X_baseline_reorder,
                                           temporal_effect_reorder=temporal_effect_reorder,
                                           eventObserved_reorder=eventObserved_reorder,
                                           time_reorder=time_reorder, estimate_hazard=estimate_hazard,
                                           indx_subset=indx_subset, K=K)

  ## iterate until converge
  while(crit && iter <= 20){

    ## calculate iterative components
    comp <- pooled_design_iter(X_baseline_reorder=X_baseline_reorder,
                               temporal_effect_reorder=temporal_effect_reorder,
                               beta=beta, indx_subset=indx_subset, K=K)

    ## beta_new
    beta_new <- solve(comp$design_information) %*% (comp$design_information %*% beta + design_matvec_Xy - comp$design_matvec_Xmu)

    ## stopping rule
    iter <-  iter + 1
    crit <- max(abs(beta_new-beta)/abs(beta)) > 1e-5

    ## update value
    beta <- beta_new

    ## clear workspace
    rm(list=c("comp","beta_new"))

    print(iter)
  }

  ## result
  return(beta)
}



#' Calculate X^T y for pooled logistic regression (this part doesn't need to iterate)
#'
#' @param X_baseline_reorder Baseline variables that won't interact with time in regression,
#'                  sparse matrix of class "dgTMatrix",
#'                  reorder observations into decreasing survival time,
#'                  intercept included
#' @param temporal_effect_reorder Baseline variables that will interact with time in regression,
#'                        sparse matrix of class "dgTMatrix" or matrix,
#'                        reorder observations into decreasing survival time,
#'                        intercept included
#' @param y Outcome variable in the pooled logistic regression
#' @param indx_subset Subset index for each time point
#' @export result
pooled_design_matvec <- function(X_baseline_reorder, temporal_effect_reorder, eventObserved_reorder, time_reorder, estimate_hazard, indx_subset, K){

  ## container
  result_Xy <- rep(0, length=dim(X_baseline_reorder)[2])
  result_temporaly <- rep(0, length=dim(temporal_effect_reorder)[2])

  ## loop over each time point
  for (i in 1:K){
  ## create outcome variable
    if (estimate_hazard == "survival"){
      y <- eventObserved_reorder * (time_reorder == i)
    }else if(estimate_hazard == "censoring"){
      y <- (1 - eventObserved_reorder) * (time_reorder == i)
    }
  ## X^T y
  result_Xy <- result_Xy + t(X_baseline_reorder[1:indx_subset[i], , drop=FALSE])%*%y[1:indx_subset[i]]
  result_temporaly <- result_temporaly + i*t(temporal_effect_reorder[1:indx_subset[i], , drop=FALSE])%*%y[1:indx_subset[i]] ## could add functions to i
  }
  ## result
  return(c(result_Xy[, 1], result_temporaly[, 1]))
}



#' Calculate X^T y, X^T mu, D and X^T diag(D) X for pooled logistic regression
#'
#' @param X_baseline_reorder Baseline variables that won't interact with time in regression,
#'                  sparse matrix of class "dgTMatrix",
#'                  reorder observations into decreasing survival time,
#'                  intercept included
#' @param temporal_effect_reorder Baseline variables that will interact with time in regression,
#'                        sparse matrix of class "dgTMatrix" or matrix,
#'                        reorder observations into decreasing survival time,
#'                        intercept included
#' @param beta Current iteration of coefficient value
#' @param indx_subset Subset index for each time point
#' @export result A list
pooled_design_iter <- function(X_baseline_reorder, temporal_effect_reorder, beta, indx_subset, K){

  ## container
  result_info <- matrix(0, nrow=(dim(temporal_effect_reorder)[2]+dim(X_baseline_reorder)[2]),
                        ncol=(dim(temporal_effect_reorder)[2]+dim(X_baseline_reorder)[2]))
  result_Xmu <- rep(0, length=dim(X_baseline_reorder)[2])
  result_temporalmu <- rep(0, length=dim(temporal_effect_reorder)[2])

  ## loop over each time point
  for (i in 1:K){
    ## mu
    temp_mu <- 1/(1+exp(-X_baseline_reorder[1:indx_subset[i], , drop=FALSE]%*%beta[1:dim(X_baseline_reorder)[2]]-i*temporal_effect_reorder[1:indx_subset[i], , drop=FALSE]%*%beta[(dim(X_baseline_reorder)[2]+1):length(beta)]))
    ## X^T mu
    result_Xmu <- result_Xmu + t(X_baseline_reorder[1:indx_subset[i], , drop=FALSE])%*%temp_mu[, 1]
    result_temporalmu <- result_temporalmu + i*t(temporal_effect_reorder[1:indx_subset[i], , drop=FALSE])%*%temp_mu[, 1]
    ## X^T diag(D) X
    temp_X <- t(X_baseline_reorder[1:indx_subset[i], , drop=FALSE])%*%((temp_mu[, 1])*X_baseline_reorder[1:indx_subset[i], , drop=FALSE])
    temp_temporal <- (i^2)*t(temporal_effect_reorder[1:indx_subset[i], , drop=FALSE])%*%((temp_mu[, 1])*temporal_effect_reorder[1:indx_subset[i], , drop=FALSE])
    temp_Xtemporal <- i*t(X_baseline_reorder[1:indx_subset[i], , drop=FALSE])%*%((temp_mu[, 1])*temporal_effect_reorder[1:indx_subset[i], , drop=FALSE])
    result_info <- result_info + cbind(rbind(temp_X, t(temp_Xtemporal)), rbind(temp_Xtemporal, temp_temporal))

    ## clear workspace
    rm(list=c("temp_mu", "temp_X", "temp_temporal", "temp_Xtemporal"))
  }
  ## result
  return(list(design_matvec_Xmu=c(result_Xmu[, 1], result_temporalmu[, 1]),
              design_information=result_info))
}





#' Prediction for pooled logistic regression
#' @export result probability, in the order of: t1(id1, id2....), t2(id1, id2....),.....

predict_pooled <- function(coef, X_baseline, temporal_effect, K){

  ## Add intercept term to X_baseline and temporal_effect
  X_baseline <- cbind(rep(1, dim(X_baseline)[1]), X_baseline)
  temporal_effect <- cbind(rep(1, dim(X_baseline)[1]), temporal_effect)

  ## predict
  LP1 <- X_baseline%*%coef[1:dim(X_baseline)[2]]
  LP2 <- temporal_effect%*%coef[(dim(X_baseline)[2]+1):length(coef)]
  LP <- rep(LP1, K)+rep(1:K, each=dim(X_baseline)[1])*rep(LP2, K)
  p <- exp(LP)/(1+exp(LP))

  ## result
  return(p)
}














