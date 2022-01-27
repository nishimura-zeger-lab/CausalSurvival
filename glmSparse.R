#' Coefficients for logistic regression with sparse matrix (via iterative reweighted least square)
#'

coefSparse <- function(outcome, offset, H, maxiter, threshold, initial_coef, printIter){

  # ## subset
  # outcome <- outcome[subset]
  # offset <- qlogis(offset)[subset]
  # H <- H[subset, ]

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

  ## initial log likelihood
  logitProb <- offset + H %*% beta
  predictedProb <- 1 / (1 + exp(-logitProb))
  logLik_comp <- outcome * log(predictedProb) + (1 - outcome) * log(1-predictedProb)
  logLik_comp[is.infinite(logLik_comp)] <- NA
  logLik <- sum(logLik_comp, na.rm=TRUE)

  rm(list=c("logitProb", "predictedProb", "logLik_comp"))

  ## iterate until converge
  while((!converged) && iter <= maxiter){

    ## calculate iterative components
    logitProb_mu <- offset + H %*% beta
    mu <- 1 / (1 + exp(-logitProb_mu))
    ## mu(1-mu)
    diagmu <- mu[, 1]*(1-mu[, 1])
    fisher_info <- Matrix::t(H) %*% (diagmu * H)

    ## beta_new
    beta_new <- solve(fisher_info, (fisher_info %*% beta + design_matvec_Xy - Matrix::t(H) %*% mu[, 1])[, 1])

    ## new log likelihood
    logitProb <- offset + H %*% beta_new
    predictedProb <- 1 / (1 + exp(-logitProb))
    logLik_comp <- outcome * log(predictedProb) + (1 - outcome) * log(1-predictedProb)
    logLik_comp[is.infinite(logLik_comp)] <- NA
    logLik_new <- sum(logLik_comp, na.rm=TRUE)


    ## stopping rule
    iter <-  iter + 1
    converged <- (abs(logLik_new-logLik)/abs(logLik_new) <= threshold)

    ## update value
    beta <- beta_new
    logLik <- logLik_new

    rm(list=c("beta_new", "logLik_new", "logLik_comp", "logitProb_mu", "mu", "fisher_info", "logitProb", "predictedProb"))

    if(printIter){print(iter - 1)}
  }

  ## result
  return(beta)
}



