#' Estimate cross-fitted TMLE of survival probability
#'
#' @param pscoreAdjustment Takes a value of either 'stratification' or 'spline'.
#'   'match' may be supported in the future.
#' @export

estimateTMLEprob <- function(eventTime, censorTime, treatment, covariants,
                             J, h.estimate="glm",  gR.estimate="glm", gA.estimate="LASSO"){

}






#' Estimate cross-fitted nuisance parameter
#'
#' @param J For cross-fitting
#' @export

estimateNuisance <- function(eventTime, censorTime, treatment, covariants,
                             J, h.estimate="glm",  gR.estimate="glm", gA.estimate="LASSO"){

}








#' Index for cross-fitting
#'
#' @param J For cross-fitting
#' @export

crossFit <- function(eventTime, censorTime, treatment, J){

  ## outcome
  outcome <- ifelse(is.na(eventTime), 0, 1)

  ## divide data into J groups with equal percentage of outcome
  set.seed(08082021)

  n_folds <- J
  rowId <- 1:length(outcome)

  index_event <- rowId[which(outcome==1)]
  index_noevent <- rowId[which(outcome==0)]
  nid_event <- length(index_event)
  nid_noevent <- length(index_noevent)

  index_ls_event <- split(sample(index_event, size=nid_event, replace=FALSE), rep(1:n_folds, each=ceiling(nid_event/n_folds))[1:nid_event])
  index_ls_noevent <- split(sample(index_noevent, size=nid_noevent, replace=FALSE), rep(1:n_folds, each=ceiling(nid_noevent/n_folds))[1:nid_noevent])

  index_ls <- lapply(1:n_folds, function(x) c(index_ls_event[[x]], index_ls_noevent[[x]]))

  return(index_ls)
}








