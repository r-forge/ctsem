#' ctCI
#' Computes confidence intervals on specified parameters / matrices for already fitted ctsem fit object.
#' @param ctfitobj Already fit ctsem fit object (class: ctsemFit) to estimate confidence intervals for.
#' @param confidenceintervals character vector of matrices and or parameters for which to estimate 95\% confidence intervals for.
#' @param optimizer character vector. Defaults to NPSOL (recommended), but other optimizers available within OpenMx (e.g. 'SLSQP') may be specified.
#' @param verbose Integer between 0 and 3 reflecting amount of output while calculating.
#' @return ctfitobj, with confidence intervals included.
#' @details Confidence intervals typically estimate more reliably using the proprietary NPSOL optimizer available within OpenMx only when
#' installing directly from OpenMx website. Use command " source('http://openmx.psyc.virginia.edu/getOpenMx.R') " to install OpenMx with NPSOL.
#' If estimating for a multigroup model, specify confidence intervals as normal, e.g. \code{confidenceintervals = c('DRIFT', 'diffusion_Y1_Y1')} .
#' The necessary group prefixes are added internally.
#' 
#' @examples
#' data('ctExample1')
#' 
#' example1model <- ctModel(n.latent = 2, n.manifest = 2, Tpoints = 6, 
#'   manifestNames = c('LeisureTime', 'Happiness'), 
#'   latentNames = c('LeisureTime', 'Happiness'), LAMBDA = diag(2))
#'   
#' example1fit <- ctFit(datawide = ctExample1, ctmodelobj = example1model)
#' 
#' example1cifit <- ctCI(example1fit, confidenceintervals = 'DRIFT')
#' 
#` summary(example1cifit)$confidenceIntervals
#' @export
ctCI<-function(ctfitobj, confidenceintervals, optimizer='NPSOL', verbose=0){
  
  if(class(ctfitobj)!= 'ctsemFit' & class(ctfitobj)!= 'ctsemMultigroupFit') stop('Not a ctsem fit object!')
  
  if(class(ctfitobj)=='ctsemMultigroupFit') confidenceintervals<-unlist(lapply(confidenceintervals,function(x){
    paste0(names(ctfitobj$mxobj$submodels),'.', confidenceintervals)
  }))
  
  originalOptimizer<- mxOption(NULL, "Default optimizer")
  
  if(optimizer=='NPSOL'){
  if(imxHasNPSOL()==TRUE) mxOption(NULL,'Default optimizer', 'NPSOL')
  if(imxHasNPSOL()==FALSE) warning("NPSOL optimizer not available - recommend installing OpenMx using command:  source('http://openmx.psyc.virginia.edu/getOpenMx.R') ")
  }
  
  if(optimizer!='NPSOL') {
    mxOption(NULL,'Default optimizer', optimizer)
  }
    
  # if(class(ctfitobj)=='ctsemFit'){
  ctfitobj$mxobj <- OpenMx::mxModel(ctfitobj$mxobj, 
    mxCI(confidenceintervals, interval = 0.95, type = "both"),
  mxComputeSequence(list(
    mxComputeConfidenceInterval(plan=mxComputeGradientDescent(nudgeZeroStarts=FALSE, 
      maxMajorIter=3000),
      constraintType=ifelse(mxOption(NULL, "Default optimizer") == 'NPSOL','none','ineq')),
    mxComputeNumericDeriv(), mxComputeStandardError(), 
    mxComputeReportDeriv()))
  )
  # }
  
#   if(class(ctfitobj)=='ctsemMultigroupFit'){
#     for(i in names(ctfitobj$mxobj$submodels)){
#     submodel <- OpenMx::mxModel(ctfitobj$mxobj$submodels[[i]], 
#       mxCI(confidenceintervals, interval = 0.95, type = "both"),
#       mxComputeSequence(list(
#         mxComputeConfidenceInterval(plan=mxComputeGradientDescent(nudgeZeroStarts=FALSE, 
#           verbose=verbose,
#           maxMajorIter=3000),
#           constraintType=ifelse(mxOption(NULL, "Default optimizer") == 'NPSOL','none','ineq')),
#         mxComputeNumericDeriv(), mxComputeStandardError(), 
#         mxComputeReportDeriv()))
#     )
#     ctfitobj$mxobj<-OpenMx::mxModel(ctfitobj$mxobj, remove=T, i)
#     ctfitobj$mxobj<-OpenMx::mxModel(ctfitobj$mxobj, submodel)
#     }
#   }

  ctfitobj$mxobj<-try(mxRun(ctfitobj$mxobj,intervals=TRUE,suppressWarnings=T,silent=T))
  ctfitobj$mxobj@compute<-NULL
  mxOption(NULL,'Default optimizer', originalOptimizer)
  
  return(ctfitobj)
}