#' Estimate cross-fitted TMLE of survival probability
#'
#' @param covariants Design matrix in triplet format (row index, col index, and value)
#' @param J For cross-fitting: random partition of subjects into J prediction sets of approximately the same size.
#' @param h.estimate Model for estimating nuisance parameter: survival hazards. Options currently include glm, xgboost
#' @param gR.estimate Model for estimating nuisance parameter: censoring hazards. Options currently include glm, xgboost
#' @param gA.estimate Model for estimating nuisance parameter: treatment probability. Options currently include LASSO
#' @param tau Time of interest
#' @export
