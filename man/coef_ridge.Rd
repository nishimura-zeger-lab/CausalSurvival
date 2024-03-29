% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PooledLogistic.R
\name{coef_ridge}
\alias{coef_ridge}
\title{“evidence maximization” approach for tunning the penalty parameter}
\usage{
coef_ridge(
  X_baseline,
  temporal_effect,
  time,
  eventObserved,
  timeIntMidPoint,
  offset_t,
  offset_X,
  weight,
  is.temporal,
  timeEffect,
  evenKnot,
  penalizeTimeTreatment,
  intercept,
  estimate_hazard,
  sigma,
  maxiter,
  threshold,
  printIter
)
}
\arguments{
\item{X_baseline}{Baseline variables that won't interact with time in regression,
sparse matrix of class "dgTMatrix".
Rows are ordered into decreasing survival time.}

\item{temporal_effect}{Baseline variables that will interact with time in regression,
sparse matrix of class "dgTMatrix" or matrix.
Rows are ordered into decreasing survival time}

\item{time}{Observed survival time. Ordered into decreasing observed survival time}

\item{eventObserved}{Event indicator. Ordered into decreasing observed survival time}

\item{is.temporal}{Whether there is temporal effect, i.e. whether time t is a variable in the regression}

\item{timeEffect}{Functions of time in the discrete censoring hazards model.
Options currently include "linear", "ns", NULL (if is.temporal = FALSE)}

\item{estimate_hazard}{"survival" or "censoring"}

\item{sigma}{A range of penalized parameters for ridge regression, a vector.}

\item{maxiter}{Maximum iterations}

\item{threshold}{Threshold for convergence}

\item{printIter}{TRUE/FALSE. Whether to print iterations or not}
}
\description{
“evidence maximization” approach for tunning the penalty parameter
}
