% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HelperFunction.R
\name{transformData}
\alias{transformData}
\title{Transform survival data from wide-format to long-format}
\usage{
transformData(dwide, timeIntMidPoint, type)
}
\arguments{
\item{dwide}{Wide-format survival data with columns: time (observed time), eventObserved (Observed event), id}

\item{type}{"survival" or "censoring"}

\item{freqTime}{Coarsen observed time to XXX days intervals}
}
\value{
A long-format survival data (with coarsening if freqTime > 1)
with columns: t (time points), It, Jt, Rt, Lt (four indicator functions) and other covariates
}
\description{
Transform survival data from wide-format to long-format
}
