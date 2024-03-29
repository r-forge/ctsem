% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/summary.ctsemFit.R
\name{summary.ctsemFit}
\alias{summary.ctsemFit}
\title{Summary function for ctsemFit object}
\usage{
\method{summary}{ctsemFit}(object, ridging = FALSE, timeInterval = 1,
  verbose = FALSE, ...)
}
\arguments{
\item{object}{ctsemFit object as generated by ctFit.}

\item{ridging}{if TRUE, adds a small amount of variance to diagonals when calculating standardised (correlation) matrices,
should only be used if standardised matrices return NAN.}

\item{timeInterval}{positive numeric value specifying time interval to use for discrete parameter matrices, defaults to 1.}

\item{verbose}{Logical. If TRUE, displays the raw, internally transformed (when fitting with default arguments) OpenMx parameters and corresponding standard errors, as well as
additional summary matrices. Parameter transforms are described in the vignette, \code{vignette('ctsem')}. Additional summary matrices
include: 'discrete' matrices -- matrices representing
the effect for the given time interval (default of 1); 'asymptotic' matrices -- represents the effect as time interval
approaches infinity (therefore asymCINT describes mean level of processes at the asymptote, asymDIFFUSION describes total within-
subject variance at the asymptote, etc); 'standardised' matrices -- transforms covariance matrices to correlation matrices, and transforms
discreteDRIFT based on DIFFUSION, to give effect sizes.}

\item{...}{additional parameters to pass.}
}
\value{
Summary of ctsemFit object
}
\description{
Provides summary details for ctsemFit objects.
}
\details{
Important: Although \code{\link{ctModel}} takes cholesky decomposed variance-covariance matrices as input,
the summary function displays the full variance-covariance matrices. These can be cholesky decomposed for comparison purposes using
\code{t(chol(summary(ctfitobject)$covariancematrix))}.
}
\examples{
### example from Driver, Oud, Voelkle (2015),
### simulated happiness and leisure time with unobserved heterogeneity.
data(ctExample1)
traitmodel <- ctModel(n.manifest=2, n.latent=2, Tpoints=6, LAMBDA=diag(2),
  manifestNames=c('LeisureTime', 'Happiness'),
  latentNames=c('LeisureTime', 'Happiness'), TRAITVAR="auto")
traitfit <- ctFit(datawide=ctExample1, ctmodelobj=traitmodel)
summary(traitfit,timeInterval=1)
}

