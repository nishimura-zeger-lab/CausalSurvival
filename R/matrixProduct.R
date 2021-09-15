#' Calculate Xv
#'
#' @param temporal_effect vector
#' @param X_baseline sparse Matrix of class "dgTMatrix"
#' @export result A vector
pooled_design_matvec <- function(temporal_effect, X_baseline){
    result <- rep(0, dim(X_baseline)[1])
    n_nonzero <- length(X_baseline@x)
    for (k in 1:n_nonzero) {
      result[X_baseline@i[k] + 1] <- result[X_baseline@i[k] + 1] + X_baseline@x[k] * temporal_effect[X_baseline@j[k] + 1]
    }
    return(result)
}



#' Calculate X^T diag(W) X
#'
#' @param X_baseline sparse Matrix of class "dgTMatrix"
#' @param weight vector
#' @export result A matrix
pooled_design_information <- function(X_baseline, weight){
  result <- as(sparseMatrix(i={},j={},dims=c(dim(X_baseline)[2],dim(X_baseline)[2])),"dgCMatrix")
  X_baseline2 <- as(sparseMatrix(i={},j={},dims=c(dim(X_baseline)[1],dim(X_baseline)[2])),"dgCMatrix")
  n_nonzero <- length(X_baseline@x)
  for (k in 1:n_nonzero){
    X_baseline2[(X_baseline@i[k] + 1), (X_baseline@j[k] + 1)] <- X_baseline2[(X_baseline@i[k] + 1), (X_baseline@j[k] + 1)] + X_baseline@x[k]*weight[X_baseline@j[k] + 1]
  }
  # for (m in ){
  #   for (n in )
  #     result[(X_baseline@j[k] + 1), (X_baseline@j[k] + 1)] <- result[(X_baseline@j[k] + 1), (X_baseline@j[k] + 1)] + X_baseline@x[k] *
  # }
  return(result)
  }













