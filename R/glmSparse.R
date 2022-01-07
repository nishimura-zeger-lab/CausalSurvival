#' Coefficients for logistic regression with sparse matrix (via iterative reweighted least square)
#'

coefSparse <- function(outcome, offset, H, maxiter, threshold){

  ## subset
  outcome <- outcome[which(dlong$It == 1)]
  offset <- qlogis(offset)[which(dlong$It == 1)]
  H <- H[which(dlong$It == 1), ]
  
  ## initial value
  converged <- FALSE
  iter <- 1
  beta <- rep(0, length=dim(H)[2])
  
  ## calculate X^T y
  design_matvec_Xy <- Matrix::t(H) %*% outcome
  
  ## initial residual deviance
  logitProb <- offset + H %*% beta
  predictedProb <- exp(logitProb) / (1 + exp(logitProb))
  dev_resid <- sum((outcome-predictedProb)^2)
  
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
    dev_resid_new <- sum((outcome-predictedProb)^2)

    ## stopping rule
    iter <-  iter + 1
    converged <- (abs(dev_resid_new-dev_resid)/abs(dev_resid_new) <= threshold)
    
    ## update value
    beta <- beta_new
    dev_resid <- dev_resid_new
    
    rm(list=c("beta_new", "dev_resid_new"))
    
    print(iter)
  }
  
  ## result
  return(list(estimates=beta, sd=sqrt(diag(solve(fisher_info)))))
}



