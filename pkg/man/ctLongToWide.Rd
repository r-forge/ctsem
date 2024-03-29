% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/ctLongtowide.r
\name{ctLongToWide}
\alias{ctLongToWide}
\title{ctLongToWide
Restructures time series / panel data from long format to wide format for ctsem analysis}
\usage{
ctLongToWide(datalong, id, time, manifestNames, TDpredNames = NULL,
  TIpredNames = NULL)
}
\arguments{
\item{datalong}{dataset in long format, including subject/id column, observation time
(or change in observation time, with 0 for first observation) column,
indicator (manifest / observed) variables,
any time dependent predictors, and any time independent predictors.}

\item{id}{character string giving column name of the subject/id column}

\item{time}{character string giving column name of the time columnn}

\item{manifestNames}{vector of character strings giving column names of manifest indicator variables}

\item{TDpredNames}{vector of character strings giving column names of time dependent predictor variables}

\item{TIpredNames}{vector of character strings giving column names of time independent predictor variables}
}
\description{
ctLongToWide
Restructures time series / panel data from long format to wide format for ctsem analysis
}
\details{
Time column must be numeric
}
\examples{
#First load the long format data with absolute times
 data('longexample')

 #Then convert to wide format
 wideexample <- ctLongToWide(datalong = longexample, id = "subject",
 time = "Time", manifestNames = c("Y1", "Y2", "Y3"),
 TDpredNames = "TD1", TIpredNames = c("TI1", "TI2"))

 #Then convert the absolute times to intervals, using the Tpoints reported from the prior step.
 wide <- ctIntervalise(datawide = wideexample, Tpoints = 3, n.manifest = 3,
 n.TDpred = 1, n.TIpred = 2, manifestNames = c("Y1", "Y2", "Y3"),
 TDpredNames = "TD1", TIpredNames = c("TI1", "TI2") )
}
\seealso{
\code{\link{ctIntervalise}}
}

