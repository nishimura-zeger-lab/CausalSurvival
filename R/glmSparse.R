#' Coefficients for logistic regression with sparse matrix (via iterative reweighted least square)
#'

coefSparse <- function(outcome, offset, H, maxiter, threshold, initial_coef, printIter){

  H <- Matrix::sparseMatrix(i = Matrix::summary(H)$i, j = Matrix::summary(H)$j,
                            x = Matrix::summary(H)$x, repr = "R")

  ## initial value
  converged <- FALSE
  iter <- 1
  if(is.null(initial_coef)){
    beta <- rep(0, length=dim(H)[2])
  }else{
    beta <- initial_coef
  }

  ## calculate X^T y
  design_matvec_Xy <- computeSubsetSparseMatVec(X=H, v=outcome, subsetSize=dim(H)[1], transposed=TRUE)

  ## initial log likelihood
  logitProb <- offset + computeSubsetSparseMatVec(X=H, v=beta, subsetSize=dim(H)[1], transposed=FALSE)
  predictedProb <- 1 / (1 + exp(-logitProb))
  logLik_comp <- outcome * log(predictedProb) + (1 - outcome) * log(1-predictedProb)
  logLik_comp[is.infinite(logLik_comp)] <- NA
  logLik <- sum(logLik_comp, na.rm=TRUE)

  rm(list=c("logitProb", "predictedProb", "logLik_comp"))

  ## iterate until converge
  while((!converged) && iter <= maxiter){

    ## calculate iterative components
    logitProb_mu <- offset + computeSubsetSparseMatVec(X=H, v=beta, subsetSize=dim(H)[1], transposed=FALSE)
    mu <- 1 / (1 + exp(-logitProb_mu))
    ## mu(1-mu)
    diagmu <- mu*(1-mu)
    fisher_info <- computeSubsetSparseInformationMatrix(X=H, weight=diagmu, subsetSize=dim(H)[1])

    ## beta_new
    beta_new <- solve(fisher_info, (fisher_info %*% beta + design_matvec_Xy - computeSubsetSparseMatVec(X=H, v=mu, subsetSize=dim(H)[1], transposed=TRUE))[, 1])

    ## new log likelihood
    logitProb <- offset + computeSubsetSparseMatVec(X=H, v=beta_new, subsetSize=dim(H)[1], transposed=FALSE)
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



