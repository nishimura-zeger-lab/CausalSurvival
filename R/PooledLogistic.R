#' Calculate X^T y, X^T mu, D and X^T diag(D) X for pooled logistic regression
#'
#' @param X_baseline Baseline variables that won't interact with time in regression,
#'                  sparse matrix of class "dgTMatrix"
#' @param temporal_effect Baseline variables that will interact with time in regression,
#'                        sparse matrix of class "dgTMatrix" or matrix
#' @param eventObserved Binary outcome variable
#' @param time Observed time
#' @param id Subject id
#' @export result A vector
pooled_design <- function(X_baseline, temporal_effect, beta, eventObserved, time, id){

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
  result_Xy <- rep(0, length=dim(X_baseline_reorder)[2])
  result_temporaly <- rep(0, length=dim(temporal_effect_reorder)[2])
  result_info <- matrix(0, nrow=(dim(temporal_effect_reorder)[2]+dim(X_baseline_reorder)[2]),
                        ncol=(dim(temporal_effect_reorder)[2]+dim(X_baseline_reorder)[2]))
  result_mu <- c()
  result_Xmu <- rep(0, length=dim(X_baseline_reorder)[2])
  result_temporalmu <- rep(0, length=dim(temporal_effect_reorder)[2])

  ## loop over each time point
  for (i in 1:K){
    ## X^T y
    result_Xy <- result_Xy + t(X_baseline_reorder[indx_subset[i], ])%*%time_reorder[indx_subset[i]]
    result_temporaly <- result_temporaly + i*t(temporal_effect_reorder[indx_subset[i], ])%*%time_reorder[indx_subset[i]] ## could add functions to i
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
  return(list(design_matvec=c(result_Xy, result_temporaly),
              design_matvec2=c(result_Xmu, result_temporalmu),
              design_information=result_info, mu=result_mu,
              indx_reorder=indx, indx_subset=indx_subset))
}



#' Coefficients for pooled logistic regression (via iterative reweighted least square)
#'
#'
coef_pooled <- function(design_matvec, design_matvec2, design_information, mu, eventObserved, indx_reorder, indx_subset){

  ## order observed data into decreasing survival time
  eventObserved_reorder <- eventObserved[indx_reorder]

  ## iterate until converge
  crit <- TRUE
  iter <- 1
  beta <- rep(0, length=dim(design_information)[1])
  while(crit && iter <= 20){
    beta_new <- solve(design_information) %*% (design_information %*% beta + design_matvec - design_matvec2)

    iter <-  iter + 1
    crit <- abs(sum(beta_new-beta)/sum(beta)) > 1e-5

    ## update value
    beta <- beta_new
  }

  ## result
  return(beta)
}
















