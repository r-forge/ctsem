% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/ctPSMfit.R
\name{ctPSMfit}
\alias{ctPSMfit}
\title{ctPSMfit}
\usage{
ctPSMfit(datawide, ctmodelobj, omxStartValues = NULL, ...)
}
\arguments{
\item{datawide}{a dataset in the wide format used by ctsem.}

\item{ctmodelobj}{A ctsem model specified using \code{\link{ctModel}}. As predictors (covariates) are
not implemented currently, ensure the model contains none of these or the fit may fail.}

\item{omxStartValues}{Named vector of start values, equivalent to this argument for \code{\link{ctFit}}.}

\item{...}{Additional parameters to pass to PSM.estimate}
}
\value{
PSM fit data
}
\description{
Fits a specified ctsem model (without predictors) using the PSM package.
}
\examples{
generatingModel <- ctModel(n.latent = 1, n.manifest = 1, Tpoints = 10,
LAMBDA = diag(1), DRIFT = matrix(-.3, nrow = 1),
MANIFESTVAR = diag(1),
CINT = matrix(3, 1, 1),
DIFFUSION = t(chol(diag(5, 1))))

dat <- ctGenerate(generatingModel, n.subjects=10, burnin=300)

### ctsem model and fit
ctsemModel <- ctModel(n.latent=1, n.manifest = 1, Tpoints = 10,
  LAMBDA = diag(1))
ctsemFit <- ctFit(dat, ctsemModel, stationary = c('T0VAR'))

### fit with PSM
psmFit <- ctPSMfit(dat, omxStartValues =
    omxGetParameters(ctsemFit$mxobj), ctsemModel)
}

