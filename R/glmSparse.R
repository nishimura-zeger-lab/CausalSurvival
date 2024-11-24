#' Coefficients for logistic regression with sparse matrix (via iterative reweighted least square)
#'

coefSparse <- function(outcome, offset, H, subset, maxiter, threshold, initial_coef, printIter){

  ## subset
  outcome <- outcome[subset]
  offset <- qlogis(offset)[subset]
  H <- H[subset, ]

  ## initial value
  converged <- FALSE
  iter <- 1
  if(is.null(initial_coef)){
    beta <- rep(0, length=dim(H)[2])
  }else{
    beta <- initial_coef
  }

  ## calculate X^T y
  design_matvec_Xy <- Matrix::t(H) %*% outcome

  ## initial residual deviance
  logitProb <- offset + H %*% beta
  predictedProb <- exp(logitProb) / (1 + exp(logitProb))
  logLik <- sum(outcome * log(predictedProb)) + sum((1 - outcome) * log(1-predictedProb))

  rm(list=c("logitProb", "predictedProb"))

  ## iterate until converge
  while((!converged) && iter <= maxiter){

    ## calculate iterative components
    logitProb_mu <- offset + H %*% beta
    mu <- exp(logitProb) / (1 + exp(logitProb))
    ## mu(1-mu)
    diagmu <- mu[, 1]*(1-mu[, 1])
    fisher_info <- Matrix::t(H) %*% (diagmu * H)

    ## beta_new
    beta_new <- solve(fisher_info, (fisher_info %*% beta + design_matvec_Xy - H %*% mu)[, 1])

    ## new residual
    logitProb <- offset + H %*% beta_new
    predictedProb <- exp(logitProb) / (1 + exp(logitProb))
    logLik_new <- sum(outcome * log(predictedProb)) + sum((1 - outcome) * log(1-predictedProb))

    ## stopping rule
    iter <-  iter + 1
    converged <- (abs(logLik_new-logLik)/abs(logLik_new) <= threshold)

    ## update value
    beta <- beta_new
    logLik <- logLik_new

    rm(list=c("beta_new", "logLik_new"))

    if(printIter){print(iter - 1)}
  }

  ## result
  return(beta)
}



