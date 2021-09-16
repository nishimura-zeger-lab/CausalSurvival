#' Calculate X^T y for pooled logistic regression
#'
#' @param X_baseline Baseline variables that won't interact with time in regression,
#'                  sparse matrix of class "dgTMatrix"
#' @param temporal_effect Baseline variables that will interact with time in regression,
#'                        sparse matrix of class "dgTMatrix" or matrix
#' @param eventObserved Binary outcome variable
#' @param time Observed time
#' @param id Subject id
#' @export result A vector
pooled_design_matvec <- function(X_baseline, temporal_effect, eventObserved, time, id){
  ## index for decreasing survival time
  indx <- order(time, decreasing = TRUE)
  ## order observed data into decreasing survival time
  id_reorder <- id[indx]
  eventObserved_reorder <- eventObserved[indx]
  time_reorder <- time[indx]
  X_baseline_reorder <- X_baseline[indx, ]
  temporal_effect_reorder <- temporal_effect[indx, ]
  ## Add intercept term to X_baseline_reorder and temporal_effect_reorder
  X_baseline_reorder <- cbind(rep(1, dim(X_baseline_reorder)[1]), X_baseline_reorder)
  temporal_effect_reorder <- cbind(rep(1, dim(X_baseline_reorder)[1]), temporal_effect_reorder)
  ## subset index for each time point
  K <- max(time_reorder)
  indx_subset <- sapply(1:K, function(x) sum(time_reorder>=x), USE.NAMES = FALSE)
  ## container
  result_X <- rep(0, length=dim(temporal_effect_reorder)[2])
  result_temporal <- rep(0, length=dim(X_baseline_reorder)[2])
  ## loop over each time point
  for (i in 1:K){
    result_X <- result_X + t(X_baseline_reorder[indx_subset[i], ])%*%time_reorder[indx_subset[i]]
    result_temporal <- result_temporal + i*t(result_temporal[indx_subset[i], ])%*%time_reorder[indx_subset[i]] ## could add functions to i
  }
  ## result
  return(c(result_X, result_temporal))
}



#' Calculate X^T diag(D) X for pooled logistic regression
#'
#' @param X_baseline Baseline variables that won't interact with time in regression,
#'                  sparse matrix of class "dgTMatrix"
#' @param temporal_effect Baseline variables that will interact with time in regression,
#'                        sparse matrix of class "dgTMatrix" or matrix
#' @param eventObserved Binary outcome variable
#' @param time Observed time
#' @param id Subject id
#' @export result A matrix
pooled_design_information <- function(X_baseline, temporal_effect, eventObserved, time, id){
  ## index for decreasing survival time
  indx <- order(time, decreasing = TRUE)
  ## order observed data into decreasing survival time
  id_reorder <- id[indx]
  eventObserved_reorder <- eventObserved[indx]
  time_reorder <- time[indx]
  X_baseline_reorder <- X_baseline[indx, ]
  temporal_effect_reorder <- temporal_effect[indx, ]
  ## Add intercept term to X_baseline_reorder and temporal_effect_reorder
  X_baseline_reorder <- cbind(rep(1, dim(X_baseline_reorder)[1]), X_baseline_reorder)
  temporal_effect_reorder <- cbind(rep(1, dim(X_baseline_reorder)[1]), temporal_effect_reorder)
  ## subset index for each time point
  K <- max(time_reorder)
  indx_subset <- sapply(1:K, function(x) sum(time_reorder>=x), USE.NAMES = FALSE)
  ## calculate weights in IRLS

  ## container

  ## loop over each time point

  ## result


}



#' Coefficients for pooled logistic regression (via iterative reweighted least square)

















