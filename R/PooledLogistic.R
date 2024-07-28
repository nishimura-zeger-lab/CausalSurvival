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
  temporal_effect_reorder <- temporal_effect[indx, ]
  ## clear workspace
  rm(list=c("id", "eventObserved", "time", "X_baseline", "temporal_effect"))

  ## Add intercept term to X_baseline_reorder and temporal_effect_reorder
  X_baseline_reorder <- cbind(rep(1, dim(X_baseline_reorder)[1]), X_baseline_reorder)
  temporal_effect_reorder <- cbind(rep(1, dim(X_baseline_reorder)[1]), temporal_effect_reorder)

  ## subset index for each time point
  K <- max(time_reorder)
  if (estimate_hazard == "survival"){
    indx_subset <- sapply(1:K, function(x) sum(time_reorder>=x), USE.NAMES = FALSE)
    ## create outcome variable
    y <- eventObserved_reorder * (time_reorder == i)
  }else if(estimate_hazard == "censoring"){
    indx_subset <- sapply(1:K, function(x) sum((time_reorder>x)*eventObserved_reorder+(time_reorder>=x)*(1-eventObserved_reorder)), USE.NAMES = FALSE)
    ## create outcome variable
    y <- (1 - eventObserved_reorder) * (time_reorder == i)
  }

  ## initial value
  crit <- TRUE
  iter <- 1
  beta <- rep(0, length=dim(temporal_effect_reorder)[2]+dim(X_baseline_reorder)[2])

  ## calculate X^T y
  design_matvec_Xy <- pooled_design_matvec(X_baseline_reorder=X_baseline_reorder,
                                           temporal_effect_reorder=temporal_effect_reorder,
                                           y=y, indx_subset=indx_subset)

  ## iterate until converge
  while(crit && iter <= 20){

    ## calculate iterative components
    comp <- pooled_design_iter(X_baseline_reorder=X_baseline_reorder,
                               temporal_effect_reorder=temporal_effect_reorder,
                               y=y, beta=beta, indx_subset=indx_subset)

    ## beta_new
    beta_new <- solve(design_information) %*% (comp$design_information %*% beta + design_matvec_Xy - comp$design_matvec_Xmu)

    ## stopping rule
    iter <-  iter + 1
    crit <- max(abs(beta_new-beta)/abs(beta)) > 1e-5

    ## update value
    beta <- beta_new

    ## clear workspace
    rm(list=c("comp","beta_new"))
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
pooled_design_matvec <- function(X_baseline_reorder, temporal_effect_reorder, y, indx_subset){

  ## container
  result_Xy <- rep(0, length=dim(X_baseline_reorder)[2])
  result_temporaly <- rep(0, length=dim(temporal_effect_reorder)[2])

  ## loop over each time point
  for (i in 1:K){
  ## X^T y
  result_Xy <- result_Xy + t(X_baseline_reorder[indx_subset[i], ])%*%y[indx_subset[i]]
  result_temporaly <- result_temporaly + i*t(temporal_effect_reorder[indx_subset[i], ])%*%y[indx_subset[i]] ## could add functions to i
  }
  ## result
  return(c(result_Xy, result_temporaly))
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
#' @param y Outcome variable in the pooled logistic regression
#' @param beta Current iteration of coefficient value
#' @param indx_subset Subset index for each time point
#' @export result A list
pooled_design_iter <- function(X_baseline_reorder, temporal_effect_reorder, y, beta, indx_subset){

  ## container
  result_info <- matrix(0, nrow=(dim(temporal_effect_reorder)[2]+dim(X_baseline_reorder)[2]),
                        ncol=(dim(temporal_effect_reorder)[2]+dim(X_baseline_reorder)[2]))
  result_mu <- c()
  result_Xmu <- rep(0, length=dim(X_baseline_reorder)[2])
  result_temporalmu <- rep(0, length=dim(temporal_effect_reorder)[2])

  ## loop over each time point
  for (i in 1:K){
    ## mu
    temp_mu <- 1/(1+exp(-X_baseline_reorder[indx_subset[i], ]%*%beta[1:dim(X_baseline_reorder)[2]]-i*temporal_effect_reorder[indx_subset[i], ]%*%beta[(dim(X_baseline_reorder)[2]+1):length(beta)]))
    result_mu <- c(result_mu, temp_mu)
    ## X^T mu
    result_Xmu <- result_Xmu + t(X_baseline_reorder[indx_subset[i], ])%*%temp_mu
    result_temporalmu <- result_temporalmu + i*t(temporal_effect_reorder[indx_subset[i], ])%*%temp_mu
    ## X^T diag(D) X
    temp_X <- t(X_baseline_reorder[indx_subset[i], ])%*%diag(temp_mu)%*%X_baseline_reorder[indx_subset[i], ]
    temp_temporal <- (i^2)*t(temporal_effect_reorder[indx_subset[i], ])%*%diag(temp_mu)%*%temporal_effect_reorder[indx_subset[i], ]
    temp_Xtemporal <- i*t(X_baseline_reorder[indx_subset[i], ])%*%diag(temp_mu)%*%temporal_effect_reorder[indx_subset[i], ]
    result_info <- result_info + cbind(rbind(temp_X, t(temp_Xtemporal)), rbind(temp_Xtemporal, temp_temporal))

    ## clear workspace
    rm(list=c("temp_mu", "temp_X", "temp_temporal", "temp_Xtemporal"))
  }
  ## result
  return(list(design_matvec_Xmu=c(result_Xmu, result_temporalmu),
              design_information=result_info))
}


















