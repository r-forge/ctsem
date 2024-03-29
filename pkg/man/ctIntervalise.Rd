% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/ctIntervalise.R
\name{ctIntervalise}
\alias{ctIntervalise}
\title{Converts absolute times to intervals for wide format ctsem panel data}
\usage{
ctIntervalise(datawide, Tpoints, n.manifest, n.TDpred = 0, n.TIpred = 0,
  imputedefs = F, manifestNames = "auto", TDpredNames = "auto",
  TIpredNames = "auto", digits = 5, mininterval = 0.001,
  individualRelativeTime = TRUE, startoffset = 0)
}
\arguments{
\item{datawide}{Wide format data, containing absolute time measurements,
to convert to interval time scale. Otherwise as used in \code{\link{ctFit}}.
See \code{\link{ctLongToWide}} to easily convert long format data.}

\item{Tpoints}{Maximum number of discrete time points (waves of data, or measurement occasions)
for an individual in the input data structure.}

\item{n.manifest}{number of manifest variables per time point in the data.}

\item{n.TDpred}{number of time dependent predictors in the data structure.}

\item{n.TIpred}{number of time independent predictors in the data structure.}

\item{imputedefs}{if TRUE, impute time intervals based on the measurement occasion (i.e. column)
they are in, if FALSE (default), set related observations to NA.
FALSE is recommended unless you are certain that the imputed value
(mean of the relevant time column) is appropriate.
Noise and bias in estimates will result if wrongly set to TRUE.}

\item{manifestNames}{vector of character strings giving variable names of manifest
indicator variables (without _Tx suffix for measurement occasion).}

\item{TDpredNames}{vector of character strings giving variable names of time
dependent predictor variables (without _Tx suffix for measurement occasion).}

\item{TIpredNames}{vector of character strings giving variable names of time
independent predictor variables.}

\item{digits}{How many digits to round to for interval calculations.}

\item{mininterval}{set to lower than any possible observed measurement interval,
but above 0 - this is used for filling NA values where necessary and has no
impact on estimates when set in the correct range.
(If all observed intervals are greater than 1, mininterval=1 may be a good choice)}

\item{individualRelativeTime}{if TRUE (default), the first measurement for each individual is
assumed to be taken at time 0, and all other times are adjusted accordingly.
If FALSE, new columns for an initial wave are created, consisting only of observations
which occurred at the earliest observation time of the entire sample.}

\item{startoffset}{if 0 (default) uses earliest observation as start time.
If greater than 0, all first observations are NA, with distance of
startoffset to first recorded observation.}
}
\description{
Converts absolute times to intervals for wide format ctsem panel data
}
\details{
Time column must be numeric!
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

 print(wide)
}

