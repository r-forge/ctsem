utils::globalVariables(c('DRIFTHATCH','invDRIFT','II','bigI','Ilatent','Alatent','Atrait','Strait',
  'Amanifesttrait','Smanifesttrait','Mpred','Apred','Spred','Amanifest',
  'invIminusAlatent','Mmanifest','expMean','datarow','existenceVector',
  'rowResults','ctsem.fitfunction','ctsem.penalties','FIMLpenaltyweight',
  'DRIFT','n.latent','DIFFUSION','TRAITVAR','n.TDpred','TDPREDEFFECT',
  'TDPREDMEANS','TDPREDVAR','TRAITTDPREDCOV','n.TIpred','TIPREDEFFECT',
  'TIPREDMEANS','TIPREDVAR','CINT','n.manifest','MANIFESTTRAITVAR','LAMBDA',
  'MANIFESTVAR','MANIFESTMEANS','sumMatrix',  'L3VarianceParameters',
  'groupcount','bigMeans',  'L3VarParDeviance','S','DRIFT','CINT',  
  'DIFFUSION','TRAITVAR', 'TIPREDEFFECT','TIPREDVAR','asymDIFFUSION',
  'expCov',
  'cvectorize','mxData','mxMatrix','mxAlgebra','mxOption','mxExpectationRAM',
  'mxFitFunctionML','Amanifestcov','Smanifest','mxExpectationNormal',
  'omxSelectCols','mxFitFunctionRow','mxExpectationStateSpace',
  'mxFitFunctionAlgebra','mxCI','mxMatrix','mxAlgebra','tr',
  'mxFitFunctionAlgebra','mxFitFunctionMultigroup','mxCI','omxGetParameters',
  'mxRun','omxSetParameters'))


#' Fit a ctsem object
#' 
#' This function fits continuous time SEM models specified via \code{\link{ctModel}}
#' to a dataset containing one or more subjects.
#' @param datawide the data you wish to fit a ctsem model to.
#' @param ctmodelobj the ctsem model object you wish to use, specified via the \code{\link{ctModel}} function.
#' @param nofit if TRUE, output only openmx model without fitting
#' @param confidenceintervals vector of character strings of matrices 
#' to calculate 95\% confidence intervals for.  e.g. c("DRIFT", "TRAITVAR")
#' @param TDpredtype if "impulse" (default) TDpredictors input a single shock.  
#' If "level" they alter the base level of process, in some sense a variable CINT.
#' @param objective 'auto' selects either 'Kalman', if fitting to single subject data, 
#' or 'mxFIML' for multiple subjects. For single subject data, 'Kalman' uses the \code{mxExpectationStateSpace }
#' function from OpenMx to implement the Kalman filter. 
#' For more than one subject, 'mxFIML' specifies a wide format SEM with a row of data per subject.
#' 'mxRAM' is also an option, equivalent to 'mxFIML', but slower, using the standard RAM 
#' matrix operations. 
#' See \code{\link{ctMultigroupFit}} for the possibility to apply the Kalman filter over multiple subjects)
#' @param stationary Character vector of T0 matrix names to constrain to stationarity.  
#' Defaults to c('T0TRAITEFFECT', 'T0TIPREDEFFECT'), constraining only the between subject difference effects. 
#' Can be set to NULL to force all T0 matrices to be estimated, can be 
#' set to 'all' to constrain all T0 matrices to stationarity. 
#' @param reasonable if TRUE, constrain variance parameters to positive values. 
#' (speeds optimization, can help to optimize to a global optimum, but can also hinder)
#' @param optimizer character string, defaults to the open-source 'CSOLNP' optimizer that is distributed
#' in all versions of OpenMx. However, 'NPSOL' generally performs a little better for these problems 
#' and is the recommended option. This requires that you have installed OpenMx manually, by running:
#' \code{source('http://openmx.psyc.virginia.edu/getOpenMx.R')} 
#' @param showInits if TRUE, prints the list of user specified and auto generated 
#' starting values for free parameters.
#' @param retryattempts Number of times to retry the start value randomisation and fit procedure, if non-convergance or uncertain fits occur.
#' @param iterationSummary if TRUE, outputs limited fit details after every fit attempt.
#' @param fit2beat Passes this argument to mxTryHard function, defaults to Inf.  Specifies maximum acceptable likelihood (useful
#' if nested submodels already fit, forces optimiser to re-run until improved estimate found).
#' @param carefulFit if TRUE, first fits the specified model with a penalised likelihood function 
#' to force MANIFESTVAR, DRIFT, TRAITVAR, MANIFESTTRAITVAR parameters to remain close to 0, then
#' fits the specified model normally, using these estimates as starting values. 
#' Can help with optimization when extreme parameter estimates are returned, 
#' though results in user specified start values being ignored for the final fit.
#' @param plotOptimization If TRUE, uses checkpointing for OpenMx function \code{mxRun}, set to checkpoint every iteration, 
#' output checkpoint file to working directory, then creates a plot for each parameter's values over iterations.
#' @param meanintervals Use average time intervals for each column for calculation 
#' (both faster and inaccurate to the extent that intervals vary across individuals).
#' @param discreteTimeModel Estimate a discrete time model - ignores timing information, parameter
#' estimates will correspond to those of classical vector autoregression models, 
#' OpenMx fit object will be directly output, thus ctsem summary and plot functionality will be unavailable.
#' Time dependent predictor type also becomes irrelevant.
#'
#'  @details
#'  DATA STRUCTURE:
#'  Single row per subject. Manifest variables first, grouped by measurement occasion (with later measurements to the right), 
#'  then 1st time dependent predictor (all observations 1:(Tpoints-1)), further time dependent predictors, 
#'  time intervals between observations, time independent predictors.  
#' 
#' @examples
#' \dontrun{
#' 
#' mfrowOld<-par()$mfrow
#' par(mfrow=c(2, 3))
#' 
#' ### example from Driver, Oud, Voelkle (2015), 
#' ### simulated happiness and leisure time with unobserved heterogeneity.
#' data(ctExample1)
#' traitmodel <- ctModel(n.manifest=2, n.latent=2, Tpoints=6, LAMBDA=diag(2), 
#'   manifestNames=c('LeisureTime', 'Happiness'), 
#'   latentNames=c('LeisureTime', 'Happiness'), TRAITVAR="auto")
#' traitfit <- ctFit(datawide=ctExample1, ctmodelobj=traitmodel)
#' summary(traitfit)
#' plot(traitfit, wait=FALSE)
#' 
#' ###Example from Voelkle, Oud, Davidov, and Schmidt (2012) - anomia and authoritarianism.  
#' data(AnomAuth) 
#' AnomAuthmodel <- ctModel(LAMBDA = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2), 
#' Tpoints = 5, n.latent = 2, n.manifest = 2, MANIFESTVAR=diag(0, 2), TRAITVAR = NULL) 
#' AnomAuthfit <- ctFit(AnomAuth, AnomAuthmodel)
#' summary(AnomAuthfit)
#' plot(AnomAuthfit, wait=FALSE)
#' 
#' par(mfrow=mfrowOld)
#' 
#' 
#' ### Single subject time series - using Kalman filter (OpenMx statespace expectation)
#' data('ctExample3')
#' model <- ctModel(n.latent = 1, n.manifest = 3, Tpoints = 100, 
#'  LAMBDA = matrix(c(1, 'lambda2', 'lambda3'), nrow = 3, ncol = 1), 
#'  MANIFESTMEANS = matrix(c(0, 'manifestmean2', 'manifestmean3'), nrow = 3, 
#'    ncol = 1), T0VAR = diag(1))
#' fit <- ctFit(data = ctExample3, ctmodelobj = model)
#' 
#' ###Oscillating model from Voelkle & Oud (2013). 
#' data(Oscillating)
#' inits <- setNames(c(-38,-.5,1,1,38,.9),c('cross','auto',
#' 'diffusion22','T0var11','T0var22','m2')
#' oscillatingm<-ctModel(n.latent = 2, n.manifest=1, Tpoints=11, 
#'   MANIFESTVAR=matrix(c(0), nrow=1, ncol=1), 
#'   LAMBDA=matrix(c(1, 0), nrow=1, ncol=2), 
#'   DRIFT=matrix(c(0, "cross", 1, "auto"), nrow=2, ncol=2), 
#'   CINT=matrix(0, ncol=1, nrow=2, ), 
#'   DIFFUSION=matrix(c(0, 0, 0, "diffusion22"), nrow=2, ncol=2), 
#'   startValues = inits)
#' oscillatingf<-ctFit(Oscillating, oscillatingm)
#' summary(oscillatingf)
#' plot(oscillatingf, wait=FALSE)
#' }
#' 
#' @export

ctFit  <- function(datawide, ctmodelobj, confidenceintervals = NULL, 
  TDpredtype="impulse", objective='auto', 
  stationary=c('T0TRAITEFFECT', 'T0TIPREDEFFECT'), 
  reasonable = TRUE, optimizer='CSOLNP', 
  retryattempts=12, iterationSummary=TRUE, fit2beat=Inf, carefulFit=TRUE,  
  showInits=FALSE, 
  meanintervals=FALSE, plotOptimization=F, 
  nofit = FALSE, discreteTimeModel=FALSE){
  
  
  checkOpenMx('ctFit')
  
  n.latent<-ctmodelobj$n.latent
  n.manifest<-ctmodelobj$n.manifest
  Tpoints<-ctmodelobj$Tpoints
  n.TDpred<-ctmodelobj$n.TDpred
  n.TIpred<-ctmodelobj$n.TIpred
  startValues<-ctmodelobj$startValues
  
  manifestNames<-ctmodelobj$manifestNames
  latentNames<-ctmodelobj$latentNames
  TDpredNames<-ctmodelobj$TDpredNames
  TIpredNames<-ctmodelobj$TIpredNames
  
  #determine objective based on number of rows
  if(objective=='auto' & nrow(datawide) == 1) objective<-'Kalman'
  if(objective=='auto' & nrow(datawide) > 1 & n.TIpred + n.TDpred ==0 ) objective<-'mxFIML'
  if(objective=='auto' & nrow(datawide) > 1 & n.TIpred + n.TDpred >0 ) objective<-'mxFIML'
  
  #capture arguments
  ctfitargs<- list(confidenceintervals, TDpredtype, stationary, reasonable, optimizer, 
    retryattempts, fit2beat, carefulFit, showInits, meanintervals, nofit, objective)
  names(ctfitargs)<-c('confidenceintervals', 'TDpredtype', 'stationary', 'reasonable', 'optimizer', 
    'retryattempts', 'fit2beat', 'carefulFit', 'showInits', 'meanintervals', 'nofit', 'objective')
  
  #ensure data is a matrix
  datawide<-as.matrix(datawide)
  
  ####check data contains correct number of columns (ignore if discreteTimeModel specified)
  neededColumns<-Tpoints * n.manifest+(Tpoints - 1)+n.TDpred * (Tpoints - 1)+n.TIpred
  if(ncol(datawide)!=neededColumns & discreteTimeModel != TRUE) stop("Number of columns in data (", paste0(ncol(datawide)), ") do not match model (", paste0(neededColumns), ")")
  
  if(discreteTimeModel == TRUE){
    message('discreteTimeModel==TRUE set -- any T0 stationarity constraints will be removed, carefulFit disabled,
timing information ignored. Parameter estimates will *not* correspond with those from continous time models.')
    carefulFit<-FALSE
    stationary<-NULL    
  }
  
  
  ####0 variance predictor fix
  if(n.TDpred>0 & objective != 'Kalman'){ #check for 0 variance predictors for random predictors implementation (not needed for Kalman because fixed predictors)
    if(any(diag(var(datawide[, paste0(TDpredNames, '_T', rep(0:(Tpoints-2), each=n.TDpred))]))==0) & 
        all(is.na(suppressWarnings(as.numeric(ctmodelobj$TDPREDVAR)))) ) {
      
      ctmodelobj$TDPREDVAR <- diag(.1,n.TDpred*(Tpoints-1))
#       tdpreds<-datawide[, paste0(TDpredNames, '_T', rep(0:(Tpoints-2), each=n.TDpred))]
#       problemcols<-colnames(tdpreds)[round(diag(var(tdpreds)), digits=3)==0]
#       
#       problempreds<-unlist(lapply(TDpredNames, function(x) if(any(unlist(regexec(x, problemcols))==1)) return (x)))

      message(paste0('Time dependent predictors with 0 variance and free TDPREDVAR matrix detected - fixing TDPREDVAR matrix to diagonal 0.1 matrix to 
allow estimation, but better to fix relevant portions manually.')) 
#         paste(problempreds, collapse=', ')))
#       
#       datawide[, paste0(rep(problempreds, each=(Tpoints-1)), '_T', 0:(Tpoints-2))] <- 
#         datawide[, paste0(rep(problempreds, each=(Tpoints-1)), '_T', 0:(Tpoints-2))] +
#         rnorm(length(datawide[, paste0(rep(problempreds, each=(Tpoints-1)), '_T', 0:(Tpoints-2))]), 0, .1)
    }
  }
  
  
  ### check single subject model adequately constrained and warn
  if(nrow(datawide)==1 & 'T0VAR' %in% stationary ==FALSE & 'T0MEANS' %in% stationary == FALSE & 
      all(is.na(suppressWarnings(as.numeric(ctmodelobj$T0VAR)))) & 
      all(is.na(suppressWarnings(as.numeric(ctmodelobj$T0MEANS)))) & nofit==FALSE) stop('Cannot estimate model for single individuals unless either 
        T0VAR or T0MEANS matrices are fixed, or set to stationary')
  
  if(nrow(datawide)==1){
    message('Single subject dataset specified - ignoring any specified between subject 
      variance matrices (T0TRAITEFFECT, TRAITVAR, MANIFESTTRAITVAR, TIPREDEFFECT, TIPREDVAR, TDTIPREDCOV, T0TIPREDEFFECT)')
    if(objective != "Kalman") message('Estimation could be much faster if objective="Kalman" was specified!')  
    n.TIpred<-0
  }
  
  
  
  
  ###check which extensions (e.g. traits) need to be included 
  if(any(ctmodelobj$TRAITVAR!=0) & nrow(datawide) > 0 )   traitExtension <- TRUE
  if(all(ctmodelobj$TRAITVAR==0) | nrow(datawide)==1)    traitExtension <- FALSE
  
  if(any(ctmodelobj$MANIFESTTRAITVAR != 0)) manifestTraitvarExtension <- TRUE
  if(all(ctmodelobj$MANIFESTTRAITVAR == 0)) manifestTraitvarExtension <- FALSE
  
  ## if Kalman objective, rearrange data to long format and set Tpoints to 2 (so only single discrete algebras are generated)
  if(objective=='Kalman') {
    if(nrow(datawide) > 1) stop('To use Kalman filter implementation with multiple subjects, see function ctMultigroupFit')
    
    if(n.TDpred >0){
      if(any(is.na(datawide[, paste0(TDpredNames, '_T', 0:(Tpoints-2))] ))) stop('NA predictors are not possible with Kalman objective')
    }
    if(n.TIpred >0) message('Time independent predictors are not possible with single subject data, ignoring')  
    
    
    datawide<-ctWideToLong(datawide, #if using Kalman, convert data to long format
      manifestNames=manifestNames, TDpredNames=TDpredNames, TIpredNames=TIpredNames, 
      n.manifest=n.manifest, 
      Tpoints=Tpoints, 
      n.TIpred=n.TIpred, 
      n.TDpred=n.TDpred)
    
    colnames(datawide)[which(colnames(datawide)=='dT')]<-'dT1'
    if(n.TDpred >0){
      datawide[2:(Tpoints), TDpredNames]<-datawide[1:(Tpoints-1), TDpredNames]
      datawide[1, TDpredNames]<-0
    }
    
    
    
    Tpoints<-2
    
    
    
  }
  
  
  #function to remove FFF's from value matrices and coerce to numeric (FFF originally added to label to specify fixed value)
  removeF <- function(matrices){
    
    for(i in 1:length(matrices)){
      x <-get(matrices[i])
      x$values<-matrix(as.numeric(gsub("FFF", "", x$values)), nrow = nrow(x$values))
      assign(matrices[i], x, pos = sys.frame( - 1))
    }
  }
  
  
  
  
  
  
  
  
  ####section to define continuous time matrices from ctModel input
#   continuousMatrices<-function(){  
    
    
    #function to process ctModel specification:  seperate labels and values, fixed and free, and generate start values
    processInputMatrix <- function(x, symmetric = FALSE, diagadd = 0, randomscale=0.01, addvalues=FALSE, addCharacters = NULL, ...){
    name<-names(x)
 
      inputm<-x[[1]]
      inputmFixed <- suppressWarnings(matrix(as.numeric(inputm), nrow = nrow(inputm), ncol = ncol(inputm))) #calc fixed only matrix  
      
      #label matrix
      inputmCharacter <- inputm    
      inputmCharacter[!is.na(inputmFixed)] <- NA #extract free label matrix
      labels <- ctLabel(TDpredNames=TDpredNames, TIpredNames=TIpredNames, manifestNames=manifestNames, latentNames=latentNames, matrixname=name, n.latent=n.latent, 
        n.manifest=n.manifest, n.TDpred=n.TDpred, n.TIpred=n.TIpred, Tpoints=Tpoints)
      if(any(!is.na(inputmCharacter)))  labels[!is.na(inputmCharacter)] <- inputmCharacter[!is.na(inputmCharacter)]
      
      #starting values matrix
      values <- matrix(round(rnorm(length(inputmFixed), 0, randomscale), 3), nrow = nrow(inputm), ncol = ncol(inputm)) #generate start values according to randomscale
      if(diagadd!= 0)  {values <- values+diag(diagadd, ncol(inputm))} #increase diagonals if specified
      if(addvalues[1]!=FALSE) values <- values+addvalues #add values if specified
      if(symmetric == TRUE)  {values[row(values) > col(values)] <- t(values)[col(values) < row(values)]} #set symmetric if necessary
      values[!is.na(inputmFixed)] <- inputmFixed[!is.na(inputmFixed)] #overwite with any fixed values
      
      if(!is.null(startValues)){ #if there are some startValues specified, check each part of x
        for( i in 1:length(startValues)){
          values[grepl(names(startValues)[i], inputmCharacter)] <- startValues[i]
        }
      }

      inits<-values #output inits matrix without fixed specification (for calculation, eg discrete drift inits)
      
      values[!is.na(inputmFixed)] <- paste0(addCharacters, inputm[!is.na(inputmFixed)]) #append addcharacters to fixed values    
    
      #free matrix
      free <- suppressWarnings(!is.na(matrix(as.numeric(values), nrow = nrow(values))))
      
      output<-list(values, labels, free, inits)
      names(output)<-c('values', 'labels', 'free', 'inits')
      return(output)
    }
    
    T0VAR <- processInputMatrix(ctmodelobj['T0VAR'], rows = n.latent, cols = n.latent, symmetric = TRUE, randomscale=0.01, diagadd = 10, addCharacters = "FFF")
    T0MEANS <- processInputMatrix(ctmodelobj["T0MEANS"], rows = n.latent, cols = 1, symmetric = FALSE, randomscale=1, diagadd = 0, addCharacters = "FFF")
    MANIFESTMEANS <- processInputMatrix(ctmodelobj["MANIFESTMEANS"], rows = n.manifest, cols = 1, symmetric = FALSE, randomscale=1, diagadd = 0, addCharacters = "FFF")
    LAMBDA <- processInputMatrix(ctmodelobj["LAMBDA"], rows = n.manifest, cols = n.latent, symmetric = FALSE, randomscale=.1, addvalues=1, diagadd = 0, addCharacters = "FFF")
    MANIFESTVAR <- processInputMatrix(ctmodelobj["MANIFESTVAR"], rows = n.manifest, cols = n.manifest, symmetric = TRUE, randomscale=0, diagadd = .15, addCharacters = "FFF")    
    DRIFT <- processInputMatrix(ctmodelobj["DRIFT"], rows = n.latent, cols = n.latent, symmetric = FALSE, randomscale=0.01, diagadd=-.5, addCharacters = "FFF")
    DIFFUSION <- processInputMatrix(ctmodelobj["DIFFUSION"], rows = n.latent, cols = n.latent, symmetric = TRUE, randomscale=0.01, diagadd = .1, addCharacters = "FFF")      
    CINT <- processInputMatrix(ctmodelobj["CINT"], rows = 1, cols = n.latent, randomscale=1, addCharacters = "FFF")    
    
    
    if(traitExtension == TRUE){ #if needed, process and include traits in matrices
      T0TRAITEFFECT <- processInputMatrix(ctmodelobj["T0TRAITEFFECT"], rows = n.latent, cols = n.latent, symmetric = FALSE, randomscale=.1, diagadd = 1, addCharacters = "FFF")
      TRAITVAR <- processInputMatrix(ctmodelobj["TRAITVAR"], rows = n.latent, cols = n.latent, symmetric = TRUE, diagadd = 1, randomscale=.1, addCharacters = "FFF")
    }
    
    
    if(manifestTraitvarExtension == TRUE){
      MANIFESTTRAITVAR <- processInputMatrix(ctmodelobj["MANIFESTTRAITVAR"], rows = n.manifest, cols = n.manifest, symmetric = TRUE, randomscale=.01, diagadd = .15, addCharacters = "FFF")
    }
    
    
    if(traitExtension == TRUE && !is.null(ctmodelobj$TRAITTDPREDCOV)) {
      TRAITTDPREDCOV <- processInputMatrix(ctmodelobj["TRAITTDPREDCOV"], symmetric = FALSE, randomscale=0.01, diagadd = 0, addCharacters = "FFF")
    }
    
    if (n.TDpred > 0){
      TDPREDMEANS <- processInputMatrix(ctmodelobj["TDPREDMEANS"], symmetric = FALSE, diagadd = 0, addCharacters = "FFF")
      TDPREDEFFECT <- processInputMatrix(ctmodelobj["TDPREDEFFECT"], symmetric = FALSE, diagadd = 0, randomscale=0.1, addCharacters = "FFF")
      T0TDPREDCOV <- processInputMatrix(ctmodelobj["T0TDPREDCOV"], symmetric = FALSE, diagadd = 0, randomscale=0.01, addCharacters = "FFF") 
      TDPREDVAR <- processInputMatrix(ctmodelobj["TDPREDVAR"], symmetric = TRUE, diagadd = 1, randomscale=0.01, addCharacters = "FFF")    
    }
    
    if (n.TIpred > 0){ 
      TIPREDMEANS <- processInputMatrix(ctmodelobj["TIPREDMEANS"], symmetric = FALSE, diagadd = 0, addCharacters = "FFF")
      TIPREDEFFECT <- processInputMatrix(ctmodelobj["TIPREDEFFECT"], symmetric = FALSE, diagadd = 0, randomscale=1, addCharacters = "FFF")
      T0TIPREDEFFECT <- processInputMatrix(ctmodelobj["T0TIPREDEFFECT"], symmetric = FALSE, diagadd = 0, randomscale=1, addCharacters = "FFF")
      TIPREDVAR <- processInputMatrix(ctmodelobj["TIPREDVAR"], symmetric = TRUE, diagadd = 1, randomscale=0.01, addCharacters = "FFF")    
    }
    
    if(n.TIpred > 0 & n.TDpred > 0) TDTIPREDCOV <- processInputMatrix(ctmodelobj["TDTIPREDCOV"], symmetric = FALSE, randomscale=0.01, addCharacters = "FFF")    
    
#     
#     returnAllLocals()
  #### end continuous matrix section
  
  
  
  
  
  
  ###section to define base RAM matrices
if(objective!='Kalman') { #configure matrices 
#   defineRAM <- function(){  
    
    #basic indexes to reuse, and modify if changed
    latentstart <- 1
    latentend <- n.latent * Tpoints
    latentExtent<-latentend #includes *all* latents - traits and manifest traits
  traitend <- latentend #because no traits yet
  latenttraitend <- traitend
    manifeststart <- n.latent * Tpoints+1
    manifestend <- latentend+n.manifest * Tpoints
    manifestExtent <-manifestend # includes *all* manifests - predictors also
    intervalstart <- manifestend+1
    intervalend <- manifestend+Tpoints - 1
    
    
    
    #create A and S matrices   
    A<-list()
    S<-list()
    A$values <- matrix('FFF0', nrow = manifestend, ncol = manifestend)
    A$labels <- matrix(, nrow = manifestend, ncol = manifestend)
    S$values <- matrix("FFF0", nrow = manifestend, ncol = manifestend)
    S$labels <- matrix(NA, nrow = manifestend, ncol = manifestend)
    
    
    
    
    #DRIFT constraints
    paste0("FFF", OpenMx::omxExponential(DRIFT$inits)) -> # discrete drift values goes into A$values 
      A$values[(row(A$values) - 1 - n.latent)%/%n.latent ==  # when the rows and columns, grouped by n.latent, 
          (col(A$values) - 1)%/%n.latent & row(A$values) <= latentend] #are equal, and not greater than the total number of latent variables
    
    #create expd(discrete drift) labels
    expdLabels <- cbind(paste0("discreteDRIFT_T", 
      rep(1:(Tpoints - 1), each = n.latent^2), 
      "[", seq(1, n.latent, 1), 
      ",", 
      rep(1:n.latent, each = n.latent), "]"))
    
 
    #insert discreteDRIFT_T's into A$labels    
    expdLabels-> A$labels  [(row(A$labels) - 1 - n.latent)%/%n.latent == #this groups rows by multiples of n.latent,offset by n.latent
        (col(A$labels) - 1)%/%n.latent & #this groups columns by multiples of n.latent.  when row and column groups are equal,
        row(A$labels) <= latentend] #and less than n.latent * Tpoints, then the groups are replaced by expdlabels
    
    
    
    #LAMBDA matrix
    paste0('FFF', LAMBDA$inits) ->  A$values[(row(A$values) - 1 - latentend)%/%n.manifest == #insert LAMBDA loadings into A$values
        (col(A$values) - 1)%/%n.latent & col(A$values) <= latentend]
    
    LAMBDA$ref <- matrix(paste0('LAMBDA[', rep(1:n.manifest, times=n.latent), ',', rep(1:n.latent, each=n.manifest), ']'), nrow=n.manifest)    
    LAMBDA$ref ->  A$labels[(row(A$labels) - 1 - latentend) %/% n.manifest == #this groups rows by multiples of n.manifest, offset by n.latent * Tpoints
        (col(A$labels) - 1)%/%n.latent & #this groups columns by multiples of n.latent.  when row and column groups are equal, 
        row(A$labels) <= manifestend] #and less than the total number of variables, then the groups are replaced by
    
    
    #measurement residuals
    paste0('FFF', MANIFESTVAR$inits) ->  S$values[(row(S$values) - latentend - 1)%/%n.manifest ==  
        (col(S$values) - latentend - 1) %/% n.manifest &
        col(S$values) >  latentend &
        row(S$values) >  latentend]
    
    MANIFESTVAR$ref <- matrix(paste0('MANIFESTVAR[', 
      rep(1:n.manifest, times=n.manifest), 
      ',', 
      rep(1:n.manifest, each=n.manifest), 
      ']'), 
      nrow=n.manifest) 
    
    MANIFESTVAR$ref ->   S$labels[(row(S$labels) - latentend - 1) %/% n.manifest ==  
        (col(S$labels) - latentend - 1) %/% n.manifest &
        col(S$labels) >  latentend]
    
    
    #latent (dynamic) residuals    
    paste0("FFF", DIFFUSION$inits) ->  S$values[(row(S$values) - 1) %/%n.latent == #initial DIFFUSION values with fixed coding
        (col(S$values) - 1) %/% n.latent &
        col(S$values) <= latentend]
    
    qdlabels <- paste0("Qd", rep(1:Tpoints - 1, each = n.latent^2), 
      "[", 
      seq(1, n.latent^2, 1), 
      ",1]") 
    
    
    qdlabels-> S$labels[(row(S$labels) - 1)%/%n.latent ==  # qdlabels goes into S$labels where the row and column groups
        (col(S$labels) - 1) %/% n.latent & # which are based on number of manifest variables, are equal, 
        col(S$labels) <= latentend] # and less than total number of latent.
    
    
    #initial latent residuals
    paste0('FFF', T0VAR$inits) ->  S$values[1:n.latent, 1:n.latent] 
    
    T0VAR$ref <- paste0("T0VAR[", 1:n.latent, ",", rep(1:n.latent, each=n.latent), "]")
    
    T0VAR$ref -> S$labels[1:n.latent, 1:n.latent] #input initial variance refs to Slabel matrix
    
    
    ### 3. Filter matrix
    FILTER<-list()
    FILTER$values    <- cbind(matrix(0, nrow = n.manifest * Tpoints, ncol = n.latent * Tpoints), diag(1, nrow = n.manifest * Tpoints))
    
    FILTERnamesy    <- c(paste0(rep(latentNames, Tpoints), "_T", rep(0:(Tpoints-1), each=n.latent)), 
      paste0(rep(manifestNames, Tpoints), "_T", rep(0:(Tpoints-1), each=n.manifest)))
    
    
    FILTERnamesx     <- paste0(rep(manifestNames, Tpoints), "_T", rep(0:(Tpoints-1), each=n.manifest))
    
    ### 4. M matrix
    M<-list()
    M$values    <- matrix(c(paste0("FFF", T0MEANS$inits), #insert appropriate T0MEANS values
      rep(paste0("FFF", T0MEANS$inits), each = Tpoints-1), #follow with fixed params for later time points
      rep(MANIFESTMEANS$values, times = Tpoints)), #then insert manifest intercepts
      nrow = manifestend, ncol = 1)
    
    latentMlabels <- matrix(paste0("intd", #this name 'intd' is referenced below to check which parameters to free
      0:(latentend - 1)%/%n.latent, 
      "[", 1:n.latent, ",1]"), ncol = 1)#construct continuous latent mean labels
    
    
    latentMlabels[1:n.latent, ]  <-  paste0("T0MEANS[", 1:n.latent, ",1]") #refer to T0MEANS matrix
    
    M$labels <- matrix(c(latentMlabels, rep(MANIFESTMEANS$labels, Tpoints)), nrow = manifestend) #insert manifest mean labels
    
#     returnAllLocals() #return objects from this base matrices function to parent
   }#close base matrices definition function
  
#### end RAM matrix section
  
  
  
  
  
  
  
  ###section to append RAM matrices with process traits
  if(objective!='Kalman' & traitExtension == TRUE){ #if needed, process and include traits in matrices
#   traitMatrices <- function(){
    
    #update indices
    manifeststart <- manifeststart+n.latent
    manifestend <- manifestend+n.latent
    traitstart <- latentend+1
    traitend <- latentend+n.latent
    latenttraitend <- traitend
    latentExtent <- traitend
    manifestExtent <- manifestend
    
    #function to insert n.latent rows and columns of single specified value into matrices 
    insertProcessTraitsToMatrix <- function(x, value){      
      x <- rbind(x[latentstart:latentend, ], #insert rows and columns for traits
        matrix(value, nrow = n.latent, ncol = ncol(x)), 
        x[(latentend+1):nrow(x), ])
      x <- cbind(x[, latentstart:latentend], 
        matrix(value, ncol = n.latent, nrow = nrow(x)), 
        x[, (latentend+1):ncol(x)])
      return(x)
    }
    
    #insert extra rows and columns to RAM matrices
    A$values <- insertProcessTraitsToMatrix(A$values, "FFF0") 
    A$labels <- insertProcessTraitsToMatrix(A$labels, NA)  
    S$values <- insertProcessTraitsToMatrix(S$values, "FFF0")
    S$labels <- insertProcessTraitsToMatrix(S$labels, NA)
    
    
    A$labels[(latentstart+n.latent):latentend, 
      (latentend+1):(latentend+n.latent)] <- paste0("TRAITd", #insert labels for the traits
        rep(1:(Tpoints - 1), each = n.latent), 
        "[", 1:n.latent, ",", 
        rep(1:n.latent, each = n.latent * (Tpoints - 1)), 
        "]")
    
    A$values[(latentstart+n.latent):latentend, 
      (latentend+1):(latentend+n.latent)] <- paste0('FFF', diag(.1, n.latent))
    
    #trait influence on T0
    
    A$values[latentstart:(latentstart+n.latent - 1), (latentend+1):(latentend+n.latent)] <- paste0('FFF', T0TRAITEFFECT$inits)
    #     S$values[(latentend+1):(latentend+n.latent), latentstart:(latentstart+n.latent - 1)] <- paste0('FFF', t(T0TRAITEFFECT$inits) #kept in case we return to cov rather than regression
    
    T0TRAITEFFECT$ref <- indexMatrix(symmetrical = FALSE, dimension = n.latent, starttext = "T0TRAITEFFECT[", sep = ",", endtext = "]")
    A$labels[latentstart:(latentstart+n.latent - 1), (latentend+1):(latentend+n.latent)] <- T0TRAITEFFECT$ref    #     S$labels(latentend+1):(latentend+n.latent), latentstart:(latentstart+n.latent - 1)] <- t(T0TRAITEFFECT$ref
    
    
    #trait variance
    S$values[(latentend+1):(latentend+n.latent), (latentend+1):(latentend+n.latent)] <- paste0('FFF', TRAITVAR$inits)
    TRAITVAR$ref <- matrix(paste0("TRAITVAR[", indexMatrix(symmetrical = TRUE, dimension = n.latent, sep = ","), "]"), nrow = n.latent)
    S$labels[(latentend+1):(latentend+n.latent), (latentend+1):(latentend+n.latent)] <- TRAITVAR$ref    
    #M matrices
    M$values <- matrix(c(M$values[1:(n.latent * Tpoints), ], 
      rep("FFF0", each = (n.latent)), 
      M$values[(n.latent * Tpoints+1):(n.latent * Tpoints+n.manifest * Tpoints), ]), ncol = 1)
    
    M$labels <- matrix(c(M$labels[1:(n.latent * Tpoints), ], 
      rep(NA, each = (n.latent)), 
      M$labels[(n.latent * Tpoints+1):(n.latent * Tpoints+n.manifest * Tpoints), ]), ncol = 1)
    
    
    #FILTER matrices
    FILTER$values    <- cbind(matrix(0, nrow = n.manifest * Tpoints, ncol = n.latent * Tpoints+n.latent), diag(1, nrow = n.manifest * Tpoints))
    
    FILTERnamesy    <- c(paste0(rep(latentNames, Tpoints), "_T", rep(0:(Tpoints-1), each=n.latent)), paste0(latentNames, "Trait"), 
      paste0(rep(manifestNames, Tpoints), "_T", rep(0:(Tpoints-1), each=n.manifest)))  #subtracting manifestend from latentend gets names of manifest vars from data (because index refers to matrices)
    
    #     FILTERnamesx already created in defineRAM
    
#     returnAllLocals() #return objects from this trait function to parent function  
  }#close trait matrices function
  
  
  
  
  
  
  
  ###section to append matrices with manifest trait latents
  if( objective!='Kalman' & manifestTraitvarExtension == TRUE){
    
    #update indices
    manifeststart <- manifeststart+n.manifest #adding n.manifest latent variables to matrix indices, retaining trait indices notation rather than splitting to process and manifest
    manifestend <- manifestend+n.manifest
    manifestExtent<-manifestend
    traitstart <- latentend+1
    if(traitExtension == TRUE) traitend <- traitend + n.manifest
    if(traitExtension == FALSE) traitend <- latentend + n.manifest
    manifesttraitstart <- traitend - n.manifest + 1
    
    #function to insert n.manifest rows and columns of single specified value into matrices 
    insertManifestTraitsToMatrix <- function(x, value){      
      x <- rbind(x[latentstart:(manifesttraitstart-1), ], #insert rows and columns for manifest traits
        matrix(value, nrow = n.manifest, ncol = ncol(x)), 
        x[(manifesttraitstart):nrow(x), ])
      x <- cbind(x[, latentstart:(manifesttraitstart-1)], 
        matrix(value, ncol = n.manifest, nrow = nrow(x)), 
        x[, (manifesttraitstart):ncol(x)])
      return(x)
    }
    
    #insert new columns in matrices for manifest traits
    A$labels <- insertManifestTraitsToMatrix(A$labels, NA)  
    S$values <- insertManifestTraitsToMatrix(S$values, "FFF0")
    S$labels <- insertManifestTraitsToMatrix(S$labels, NA)
    A$values <- insertManifestTraitsToMatrix(A$values, "FFF0") 
    
    #manifest trait variance
    MANIFESTTRAITVAR$ref<- indexMatrix(dimension=n.manifest, symmetrical=TRUE, 
      sep=',', starttext='MANIFESTTRAITVAR[', endtext=']')
    
    S$values[manifesttraitstart:traitend, manifesttraitstart:traitend] <- paste0('FFF', MANIFESTTRAITVAR$inits)
    S$labels[manifesttraitstart:traitend, manifesttraitstart:traitend] <- MANIFESTTRAITVAR$ref    
    #manifest trait effect on manifests
    for(i in 0:(Tpoints-1)){
      A$values[(manifeststart+(n.manifest*i)):(manifeststart+(n.manifest* (i + 1)) - 1), 
        (manifesttraitstart):(traitend)] <- paste0("FFF", diag(n.manifest))
    }
    
    #M matrices
    M$values <- matrix(c(M$values[1:(manifesttraitstart-1), ], #include M$values before MANIFESTTRAITVARs
      rep("FFF0", each = (n.manifest)), #fix MANIFESTTRAITVAR means to 0
      M$values[(traitend+1-n.manifest):(manifestend-n.manifest), ]), ncol = 1)#include M$values after MANIFESTTRAITVARs
    
    M$labels <- matrix(c(M$labels[1:(manifesttraitstart-1), ], 
      rep(NA, each = (n.manifest)), 
      M$labels[(traitend+1-n.manifest):(manifestend-n.manifest), ]), ncol = 1)
    
    
    #F matrices
    FILTER$values    <- cbind(matrix(0, nrow = n.manifest * Tpoints, ncol = traitend), diag(1, nrow = n.manifest * Tpoints))
    
    FILTERnamesy    <- c(FILTERnamesy[1:(manifesttraitstart-1)], 
      paste0(manifestNames, 'Trait'), 
      FILTERnamesy[manifesttraitstart:length(FILTERnamesy)])  #subtracting manifestend from latentend gets names of manifest vars from data (because index refers to matrices)
    
    #FILTERnamesx already created
    
#     returnAllLocals() #return objects from this manifest trait function to parent function
  }#close manifest trait matrices function
  
  
  
  
  
  
  
  ####TDpred random matrix section
  if(n.TDpred > 0 & objective!='Kalman') {
 
    
    #function to insert rows and columns of single specified value into matrices
    insertTDpredsToMatrix <- function(target, value){
      target <- rbind(target, #insert rows and columns for traits
        matrix(value, nrow = n.TDpred * (Tpoints - 1), ncol = ncol(target)))
      target <- cbind(target, 
        matrix(value, nrow = nrow(target), ncol = n.TDpred * (Tpoints - 1)))
      return(target)
    }
    
    #update indices
    
    predictorTDstart <- manifestend+1
    predictorTDend <- manifestend+n.TDpred * (Tpoints - 1)
    predictorstart<-predictorTDstart
    predictorend<-predictorTDend
    
    
    #add rows and columns to matrices
    A$values <- insertTDpredsToMatrix(A$values, "FFF0")   
    A$labels <- insertTDpredsToMatrix(A$labels, NA)  
    S$values <- insertTDpredsToMatrix(S$values, "FFF0")
    S$labels <- insertTDpredsToMatrix(S$labels, NA)
    
    
    #create time dependent predictor effects on processes
    TDPREDEFFECT$ref <- paste0(
      "TDPREDEFFECT", 
      "_T", 
      rep(1:(Tpoints - 1), each=n.latent), 
      "[", 
      1:n.latent, 
      ",", 
      rep(1:n.TDpred, each = (Tpoints - 1)*n.latent), 
      "]")
    
    
    A$values[cbind(rep( (1+n.latent):latentend, times=n.TDpred), 
      rep(predictorTDstart:predictorTDend, each=n.latent))] <- paste0("FFF", .5) #insert starting values specifying fixed for algebras
    A$labels[cbind(rep( (1+n.latent):latentend, times=n.TDpred), 
      rep(predictorTDstart:predictorTDend, each=n.latent))] <- TDPREDEFFECT$ref #insert TDPREDEFFECT algebra references to A$labels    
    
    #add cov of time dependent predictors with T0
    T0TDPREDCOV$ref <- paste0('T0TDPREDCOV[', 1:n.latent, ',', rep( 1:(n.TDpred*(Tpoints-1)), each=n.latent ), ']')
    S$values[1:n.latent, predictorTDstart:predictorTDend] <- paste0('FFF', T0TDPREDCOV$inits) #add starting values 
    S$values[predictorTDstart:predictorTDend, 1:n.latent] <- t(paste0('FFF', T0TDPREDCOV$inits)) #add starting values 
    
    S$labels[1:n.latent, predictorTDstart:predictorTDend]  <-  T0TDPREDCOV$ref #insert combined labels to S matrix
    S$labels[predictorTDstart:predictorTDend, 1:n.latent]  <-  t(T0TDPREDCOV$ref) #insert combined labels to S matrix
    
    #add cov between all td predictors 
    #     if(fastPredictors==TRUE){ #then don't optimize, just use estimates from cov
    #       temp<-matrix(paste0('FFF', cov(datawide[, paste0(rep(TDpredNames, each=(Tpoints-1)), '_T', 0:(Tpoints-2))], use='pairwise.complete.obs')), 
    #         nrow=nrow(TDPREDVAR$values))
    #       browser()
    #         TDpreds<-datawide[, paste0(rep(TDpredNames, each=(Tpoints-1)), '_T', 0:(Tpoints-2))]
    #         
    #         covmodel<-OpenMx::mxModel(mxData(cov(TDpreds), means=apply(TDpreds, 2, mean, na.rm=T), type='cov', numObs=nrow(TDpreds)), 
    #           mxMatrix(name='A', values=0, free=F, type='Full', nrow=ncol(TDpreds), ncol=ncol(TDpreds)), 
    #           mxMatrix(name='S', values=cov(TDpreds), free=T, type='Full', nrow=ncol(TDpreds), ncol=ncol(TDpreds)), 
    #           mxMatrix(name='F', values=diag(ncol(TDpreds)), free=F, type='Full', nrow=ncol(TDpreds), ncol=ncol(TDpreds)), 
    #           mxMatrix(name='M', values=apply(TDpreds, 2, mean, na.rm=T), free=T, type='Full', nrow=1, ncol=ncol(TDpreds)), 
    #           mxExpectationRAM(M='M', dimnames=colnames(TDpreds)), 
    #           mxFitFunctionML()
    #         )
    #         tempcov<-OpenMx::mxRun(covmodel, silent=TRUE)
    #         
    #   
    #       
    #       if(all(!is.na(temp))){
    #         TDPREDVAR$values<-temp
    #         TDPREDVAR$free<-FALSE  
    #       }
    #       if(any(is.na(temp))) message('Too much missingness in TD predictors to use calculated covariance for TDPREDVAR, so estimating...')
    #     }
    TDPREDVAR$ref<-paste0('TDPREDVAR[', 1:(n.TDpred*(Tpoints-1)), ',', rep( 1:(n.TDpred*(Tpoints-1)), each=n.TDpred*(Tpoints-1) ), ']')
    S$values[predictorTDstart:predictorTDend, predictorTDstart:predictorTDend]  <- paste0('FFF', TDPREDVAR$inits)  	#insert values    
    S$labels[predictorTDstart:predictorTDend, predictorTDstart:predictorTDend ] <- TDPREDVAR$ref #insert combined labels into S matrix
    
    #introduce covariance between TDpreds and traits    
    if(traitExtension == TRUE && !is.null(ctmodelobj$TRAITTDPREDCOV)){
      TRAITTDPREDCOV$ref<-paste0('TRAITTDPREDCOV[', 1:n.latent, ',', rep( 1:(n.TDpred*(Tpoints-1)), each=n.latent ), ']')
      S$values[traitstart:(latentend+n.latent), predictorTDstart:predictorTDend ]  <- paste0('FFF', TRAITTDPREDCOV$inits) #insert starting values
      S$values[predictorTDstart:predictorTDend, traitstart:(latentend+n.latent) ]  <- t(paste0('FFF', TRAITTDPREDCOV$inits))#insert symmetric starting values
      
      S$labels[traitstart:(latentend+n.latent), predictorTDstart:predictorTDend ]  <- TRAITTDPREDCOV$ref #insert combined labels to s matrix      
      S$labels[predictorTDstart:predictorTDend, traitstart:(latentend+n.latent) ]  <- t(TRAITTDPREDCOV$ref) #insert symmetric labels      
    }
    
    #Means
    M$values  <- rbind(M$values, TDPREDMEANS$values)
    M$labels  <- rbind(M$labels, TDPREDMEANS$labels)
    
    #filter matrix    
    FILTER$values    <- cbind(matrix(0, nrow = (n.manifest * Tpoints+n.TDpred * (Tpoints - 1)), ncol = manifeststart - 1), 
      diag(1, nrow = n.manifest * Tpoints+n.TDpred * (Tpoints - 1)))
    
    FILTERnamesy <- c(FILTERnamesy, #already specified FILTERnames
      paste0(rep(TDpredNames, each=(Tpoints-1)), "_T", 0:(Tpoints-2))) #TDpred names
    
    FILTERnamesx     <- c(FILTERnamesx, 
      paste0(rep(TDpredNames, each=(Tpoints-1)), "_T", 0:(Tpoints-2)))
    
#     returnAllLocals() #return objects from this function to parent function
  }#close TD predictor matrices function
  
  
  
  
  
  
  
  ######Time independent predictors random matrix section
  if(n.TIpred > 0 & objective != 'Kalman') {
    
    #function to insert rows and columns of single specified value into matrices
    insertTIpredsToMatrix <- function(target, value){
      target <- rbind(target, 
        matrix(value, nrow = n.TIpred, ncol = ncol(target)))
      target <- cbind(target, 
        matrix(value, nrow = nrow(target), ncol = n.TIpred))
      return(target)
    }
    
    #update indices
    predictorTIstart <- manifestend + n.TDpred*(Tpoints-1) + 1
    predictorTIend <- predictorTIstart+n.TIpred-1
    predictorend<-predictorTIend
    if(n.TDpred == 0) predictorstart<-predictorTIstart
    
    #insert rows and columns to matrices
    A$values <- insertTIpredsToMatrix(A$values, "FFF0")
    A$labels <- insertTIpredsToMatrix(A$labels, NA)  
    S$values <- insertTIpredsToMatrix(S$values, "FFF0")
    S$labels <- insertTIpredsToMatrix(S$labels, NA)    
    
    #Effect of TIpreds on processes
    A$values[(n.latent+1):latentend, predictorTIstart:predictorTIend] <- paste0("FFF", .1) #add rough starting values, fixed for algebras
    
    TIPREDEFFECT$ref <- paste0( #create time Tindependent predictor algebra reference labels for A matrix 
      "TIPREDEFFECT_T", 
      rep(1:(Tpoints - 1), each = n.latent), 
      "[", 
      1:n.latent, 
      ",",    
      rep(1:n.TIpred, each = (Tpoints - 1) * n.latent), 
      "]")
    
    A$labels[(n.latent+1):latentend, predictorTIstart:predictorTIend] <- TIPREDEFFECT$ref #add time Tindependent predictor algebra references to A matrix
    
    #add effect of time independent predictors on t1
    T0TIPREDEFFECT$ref <- paste0( #create time Tindependent predictor algebra reference labels for A matrix 
      "T0TIPREDEFFECT", 
      "[", 
      1:n.latent, 
      ",",    
      rep(1:n.TIpred, each=n.latent), 
      "]")
    A$values[1:n.latent, predictorTIstart:predictorTIend] <- paste0('FFF', T0TIPREDEFFECT$inits)
    A$labels[1:n.latent, predictorTIstart:predictorTIend ]  <-  T0TIPREDEFFECT$ref 
    
    
    #add cov between TIpreds
    S$values[predictorTIstart:predictorTIend, predictorTIstart:predictorTIend]  <- TIPREDVAR$values
    S$labels[predictorTIstart:predictorTIend, predictorTIstart:predictorTIend ] <- TIPREDVAR$labels 
    
    #add cov between TDpreds and TIpreds
    if(n.TDpred > 0 & n.TIpred > 0){
      #       #fast predictor estimates if fastPredictors is set
      #       if(fastPredictors==TRUE){
      #         
      #         temp<-matrix(paste0('FFF', cov(y=datawide[, TIpredNames], 
      #          x=datawide[, paste0(rep(TDpredNames, each=(Tpoints-1)), '_T', 0:(Tpoints-2))], use='pairwise.complete.obs')), 
      #           nrow=nrow(TDTIPREDCOV$values))
      #         
      #         if(all(!is.na(temp))){
      #           TDTIPREDCOV$values<-temp
      #           TDTIPREDCOV$free<-FALSE  
      #         }
      #         if(any(is.na(temp))) message('Too much missingness in TI and TD predictors to use calculated covariance for TDTIPREDCOV, so estimating...')
      #       }
      
      
      S$values[predictorTDstart:predictorTDend, predictorTIstart:predictorTIend]  <- TDTIPREDCOV$values        
      S$labels[predictorTDstart:predictorTDend, predictorTIstart:predictorTIend] <- TDTIPREDCOV$labels 
      
      S$values[predictorTIstart:predictorTIend, predictorTDstart:predictorTDend]  <- t(TDTIPREDCOV$values)
      S$labels[predictorTIstart:predictorTIend, predictorTDstart:predictorTDend] <- t(TDTIPREDCOV$labels) 
    }
    
    #TIpred means
    M$values  <- rbind(M$values, TIPREDMEANS$values)
    M$labels  <- rbind(M$labels, TIPREDMEANS$labels)
    
    #Filter matrix
    FILTER$values    <- cbind(matrix(0, nrow = (n.manifest * Tpoints+n.TDpred * (Tpoints - 1)+n.TIpred), ncol = manifeststart - 1), 
      diag(1, nrow = n.manifest * Tpoints+n.TDpred * (Tpoints - 1)+n.TIpred))
    
    FILTERnamesy <- c(FILTERnamesy, TIpredNames)      
    FILTERnamesx     <- c(FILTERnamesx, TIpredNames)
    
#     returnAllLocals() #return objects from this function to parent function
  }#close TI predictor matrices function
  
  
  
  if(objective!='Kalman') { #do final RAM matrix adjustments - remove FFF's appropriately
    A$free  <- matrix(grepl("FFF", A$values)!= TRUE, nrow = nrow(A$values))     
    S$free <- matrix(grepl("FFF", S$values)!= TRUE, nrow = nrow(S$values))
    M$free  <- suppressWarnings(!is.na(matrix(as.numeric(M$values), nrow = nrow(M$values))))
    M$free[which(grepl("FFF", M$values) == TRUE)] <- F #these are set only after other matrices are fully specified
    removeF(c("A", "S", "M"))
  }
  
  removeF(c("DRIFT", "DIFFUSION", "CINT", 'MANIFESTMEANS', 
    "LAMBDA", "T0VAR", "MANIFESTVAR", "T0MEANS"))#removes leading FFF's which indicate fixed params and coerces matrices to numeric
  if(traitExtension == TRUE) removeF(c("T0TRAITEFFECT", "TRAITVAR"))
  if(manifestTraitvarExtension==TRUE) removeF("MANIFESTTRAITVAR")
  
  if (n.TDpred > 0) removeF(c("TDPREDEFFECT", "TDPREDMEANS", "T0TDPREDCOV", "TDPREDVAR"))
  if (n.TIpred > 0) removeF(c("TIPREDEFFECT", "TIPREDMEANS", "T0TIPREDEFFECT", "TIPREDVAR"))
  if(n.TDpred >0 & n.TIpred >0) removeF(c("TDTIPREDCOV"))
  if(!is.null(ctmodelobj$TRAITTDPREDCOV) && traitExtension == TRUE)  removeF("TRAITTDPREDCOV")
  
  
  

  
  
  
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Define OpenMx RAM Algebras for the continuous time drift matrix (A), intercept (INT), and error covariance (DIFFUSION)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    
    #set up definition variables
    defcall <- paste0("data.dT", 1:(Tpoints-1)) 
    if(meanintervals==TRUE) defcall <- apply(datawide[, paste0("dT", 1:(Tpoints-1))], 2, mean)
    
    EXPalgs <- list()
    for(i in 1:(Tpoints - 1)){
      if(discreteTimeModel==FALSE) fullAlgString <- paste0("omxExponential(DRIFT %x%", defcall[i], ")")
      
      if(discreteTimeModel==TRUE) fullAlgString <- paste0("DRIFT")
      
      EXPalgs[i] <- eval(substitute(OpenMx::mxAlgebra(theExpression, name = paste0("discreteDRIFT_T", i)), 
        list(theExpression = parse(text = fullAlgString)[[1]])))
    }
    
    INTalgs <- list()
    for(i in 1:(Tpoints - 1)){
      
      if(discreteTimeModel==FALSE){
        #         if(asymptotes==FALSE) 
        fullAlgString <- paste0('invDRIFT %*%   (discreteDRIFT_T', i, ' - II) %*% CINT') #optimize based on continuous CINT
        
        #         if(asymptotes==TRUE) fullAlgString <- paste0("(II - discreteDRIFT_T", i, ") %*% CINT") #discrete cint based on asymptotic CINT
      }
      
      if(discreteTimeModel==TRUE) fullAlgString <- paste0("CINT")
      
      INTalgs[i] <- eval(substitute(OpenMx::mxAlgebra(theExpression, name = paste0("intd", i)), 
        list(theExpression = parse(text = fullAlgString)[[1]])))
    }
    
    Qdalgs <- list()
    for(i in 1:(Tpoints - 1)){
      
      if(discreteTimeModel==FALSE){
        #         if(asymptotes==FALSE) 
        fullAlgString <- paste0("solve(DRIFTHATCH) %*% 
          ((omxExponential(DRIFTHATCH %x% ", defcall[i], ")) - (II%x%II) )%*% rvectorize(DIFFUSION)") #optimize over continuous diffusion variance
        
        #         if(asymptotes==TRUE) fullAlgString <- paste0(" ( II %x% II  - 
        #           (discreteDRIFT_T", i, ") %x% (discreteDRIFT_T", i, ") )  %*%  cvectorize(DIFFUSION) ") #from asymptotic diffusion variance
      }   
      
      if(discreteTimeModel==TRUE) fullAlgString <- paste0("rvectorize(DIFFUSION)")
      
      Qdalgs[i] <- eval(substitute(OpenMx::mxAlgebra(theExpression, name = paste0("Qd", i)), 
        list(theExpression = parse(text = fullAlgString)[[1]])))
    }
    
#end base algebra definition function
  
  
  
  
  
  
  
  ####TRAIT EXTENSION ALGEBRAS
  if(traitExtension == TRUE) {  
    
    traitalgs <- list()
    for(i in 1:(Tpoints - 1)){
      
      if(discreteTimeModel==FALSE){
        #         if(asymptotes==FALSE) 
        fullAlgString <- paste0("invDRIFT %*%   (omxExponential(DRIFT %x%", defcall[i], ") - II)")   #optimize using continuous traitvar        
        #         if(asymptotes==TRUE) 
        #         fullAlgString <- paste0("II - discreteDRIFT_T", i) #using asymptotic trait variance 
      }
      
      if(discreteTimeModel==TRUE) fullAlgString <- paste0("II")
      
      traitalgs[i] <- eval(substitute(OpenMx::mxAlgebra(theExpression, name = paste0("TRAITd", i)  ), 
        list(theExpression = parse(text = fullAlgString)[[1]])))  	
    }
    
#     returnAllLocals()
  }#end Trait algebra definition function
  
  
  
  #### Predictors algebra setup
  if(n.TDpred + n.TIpred > 0){
    
    if( n.TDpred > 0 ) { #if there are TD predictors     
      TDPREDEFFECTalgs <- list()
      
      for(j in 1:(Tpoints - 1)){
        
        if(TDpredtype=="impulse"){
          if(discreteTimeModel==FALSE) fullAlgString <- paste0("discreteDRIFT_T", j, " %*% TDPREDEFFECT")
          if(discreteTimeModel==TRUE) fullAlgString <- paste0("TDPREDEFFECT")
        }
        
        if(TDpredtype=="level") {
          if(discreteTimeModel==FALSE) fullAlgString <- paste0("invDRIFT%*% (discreteDRIFT_T", j, " - II)%*%TDPREDEFFECT")
          if(discreteTimeModel==TRUE) fullAlgString <- paste0("TDPREDEFFECT")
        }
          
        
        #         if(TDpredtype=="level" & asymptotes==TRUE) fullAlgString <- paste0(        #check this after algebra negative corrections
        #           "(II - discreteDRIFT_T", j, ") %*% t(TDPREDEFFECT)")        
        
        TDPREDEFFECTalgs[j] <- eval(substitute(OpenMx::mxAlgebra(theExpression, name = paste0("TDPREDEFFECT", "_T", j)), 
          list(theExpression = parse(text = fullAlgString)[[1]]))) 
      }
    }
    
    
    if (n.TIpred > 0){ #if there are fixed time independent predictors
      TIPREDEFFECTalgs <- list()
      for(j in 1:(Tpoints - 1)){
        #         
        #         if(asymptotes==FALSE) fullAlgString <- paste0("invDRIFT %*% (omxExponential(DRIFT %x% ", defcall[i], ") - II) %*% TIPREDEFFECT")
        #         if(asymptotes==TRUE) 
        
        if(discreteTimeModel==FALSE) fullAlgString <- paste0("invDRIFT %*% (discreteDRIFT_T", j, " - II) %*% TIPREDEFFECT") 
        if(discreteTimeModel==TRUE) fullAlgString <- paste0("TIPREDEFFECT") 
        
        TIPREDEFFECTalgs[j] <- eval(substitute(OpenMx::mxAlgebra(theExpression, name = paste0("TIPREDEFFECT", "_T", j)), 
          list(theExpression = parse(text = fullAlgString)[[1]])))  
      }
    }
    
    
    
    
    #     if(randomPredictors==FALSE & n.TDpred > 0){ #if we want to insert TD predictors via mean inputs controlled by definition variables
    # 
    #       TDpreddefs<-list()
    #       TDpreddefsfull <- mxMatrix(name = 'TDpreddefsfull', #create a definition variable matrix of all td predictors
    #         labels = paste0("data.", rep(TDpredNames, each=(Tpoints-1)), "_T", 0:(Tpoints-2) ), nrow=1, ncol=n.TDpred*(Tpoints-1))#       T0TDPREDCOVmatrices<-list()
    #       for(i in 1:(Tpoints - 1)){ #for every time point except the first, create intercept algebra that includes predictors
    #         TDpreddefs <- c(TDpreddefs, mxMatrix(name = paste0('TDpreddefs_T', i), #create a definition variable matrix
    #           labels = paste0("data.", TDpredNames, "_T", i-1), nrow=1, ncol=n.TDpred)) 
    #         
    # #         T0TDPREDCOVmatrices<-c(T0TDPREDCOVmatrices, mxMatrix(name=paste0(T0TDPREDCOV_T, i-1), 
    # #           labels=T0TDPREDCOVlabels[, seq(i, (Tpoints-1)*n.TDpred, (Tpoints-1))], 
    # #           values=T0TDPREDCOV$values[, seq(i, (Tpoints-1)*n.TDpred, (Tpoints-1))], 
    # #           free=T, nrow=n.latent, ncol=n.TDpred))
    #         
    #         fullAlgString <-  paste(deparse(INTalgs[[i]]$formula), collapse='') #set base intercept algebra        
    #         
    #         if(TDpredtype=="impulse")   secondAlgString <- paste0( " + 
    #           omxExponential(DRIFT %x%", defcall[i], ") %*% (TDPREDEFFECT %*% t(TDpreddefs_T", i, ")) " )
    #         
    #         if(TDpredtype=="level") secondAlgString <-paste0( " + 
    #           (invDRIFT%*% (omxExponential(DRIFT %x%", defcall[i], ") - II) %*% 
    # TDPREDEFFECT %*% t(TDpreddefs_T", i, "))") #add the string to the intercept algebra
    #         
    #         fullAlgString <- paste0(fullAlgString, secondAlgString) #combine the strings
    #         
    #         INTalgs[i] <- eval(substitute(mxAlgebra(theExpression, name = paste0("intd", i)), 
    #           list(theExpression = parse(text = fullAlgString)[[1]])))
    #       }
    #       
    # 
    #       
    #       
    #     }# end fixed TD predictors
    #     
    #     
    #     
    #     
    #     if(randomPredictors==FALSE & n.TIpred >0){ #if we have TI predictors
    #       
    #       for (i in 1:(Tpoints-1)){
    #         fullAlgString <-  paste(deparse(INTalgs[[i]]$formula), collapse='') #set base intercept algebra 
    #         
    #         if(asymptotes==FALSE) secondAlgString <- paste0( " +  (invDRIFT%*%(omxExponential(DRIFT %x%", defcall[i], ") - II) %*% TIPREDEFFECT %*% t(TIpreddefs))") #create TI predictor algebra
    #       
    #         if(asymptotes==TRUE) secondAlgString <- paste0(" + (II - omxExponential(DRIFT %x%", defcall[i], ")) %*% TIPREDEFFECT %*% t(TIpreddefs)")        
    #         
    #         fullAlgString <- paste0(fullAlgString, secondAlgString) #combine this string with the full intercept string
    #         
    #         INTalgs[i] <- eval(substitute(mxAlgebra(theExpression, name = paste0("intd", i)), 
    #           list(theExpression = parse(text = fullAlgString)[[1]])))
    #       }
    #     } # end fixed TI predictor / intercept algebra
    
#     returnAllLocals()
  } # end predictors model section
  
  
  
  
  
  
  
  

  
  
   #base model specification
    
    zerodiagnlatent <- matrix(NA, n.latent, n.latent) #set a matrix with 0's on diag for upper / lower bounds
    diag(zerodiagnlatent)<-0.000001
    zerodiagnmanifest <- matrix(NA, n.manifest, n.manifest) #set a matrix with 0's on diag for upper / lower bounds
    diag(zerodiagnmanifest)<-0.000001
    
    if('T0VAR' %in% stationary) {
      T0VAR$labels<-paste0('asymDIFFUSION[', 1:(n.latent^2), ',1]')
      T0VAR$free<-FALSE
      asymDIFFUSIONalg<- OpenMx::mxAlgebra(name='asymDIFFUSION', -solve(DRIFTHATCH) %*% cvectorize(DIFFUSION))
    }
    if('T0MEANS' %in% stationary){
      T0MEANS$labels<-paste0('asymCINT[', 1:n.latent, ',1]')
      T0MEANS$free<-FALSE
      asymCINTalg<- OpenMx::mxAlgebra(name='asymCINT', -invDRIFT %*% CINT )      
    }
    
    
    model  <-  OpenMx::mxModel("ctsem", #type="RAM", #begin specifying the mxModel
      mxData(observed = datawide, type = "raw"), 
      
      mxMatrix(type = "Iden", nrow = n.latent, ncol = n.latent, free = FALSE, name = "II"), #identity matrix
      
      mxMatrix(name='LAMBDA', free=LAMBDA$free, values=LAMBDA$values, dimnames=list(manifestNames, latentNames), 
        labels=LAMBDA$labels, nrow=n.manifest, ncol=n.latent), 
      
      mxMatrix(name='MANIFESTVAR', free=MANIFESTVAR$free, values=MANIFESTVAR$values, 
        labels=MANIFESTVAR$labels, nrow=n.manifest, ncol=n.manifest, lbound=if(reasonable==TRUE){zerodiagnmanifest}else{NA}), 
      
      mxMatrix(name = "T0MEANS", free=T0MEANS$free, labels=T0MEANS$labels, #T0MEANS matrix
        values=T0MEANS$values, nrow=n.latent, ncol=1), 
      
      mxMatrix(name = "T0VAR", values=T0VAR$values, labels=T0VAR$labels, 
        lbound=if(reasonable==TRUE){zerodiagnlatent}else{NA}, ncol=n.latent, nrow=n.latent, free=T0VAR$free), 
      
      mxMatrix(type = "Full", labels = DIFFUSION$labels, values = DIFFUSION$values, #DIFFUSION matrix of dynamic innovations
        byrow = TRUE, free = DIFFUSION$free, name = "DIFFUSION", nrow = n.latent, ncol = n.latent, 
        lbound = if(reasonable==TRUE){zerodiagnlatent}else{NA}), 
      
      mxMatrix(type = "Full", labels = DRIFT$labels, , values = DRIFT$values, byrow = TRUE, #continuous DRIFT matrix
        free = DRIFT$free, name = "DRIFT"), 
      
      mxMatrix(type = "Full", labels = CINT$labels, values = CINT$values, 
        free = CINT$free, nrow=nrow(CINT$values), ncol=ncol(CINT$values), name = "CINT"),  #continuous intercept matrix
      INTalgs, 	#continuous intercept discrete translation algebras
      Qdalgs, # error covariance (DIFFUSION) discrete translation algebras
      mxAlgebra(name='invDRIFT', solve(DRIFT)), 
      EXPalgs #include matrix exponential algebras to equate discrete matrix to DRIFT matrix             
    )
    
    if('T0VAR' %in% stationary) model<-OpenMx::mxModel(model, asymDIFFUSIONalg)
    if('T0MEANS' %in% stationary) model<-OpenMx::mxModel(model, asymCINTalg)
    
    #     if(asymptotes==FALSE) 
    model<-OpenMx::mxModel(model, 
      mxAlgebra(DRIFT%x%II + II%x%DRIFT, name = "DRIFTHATCH"), #used in continuous DIFFUSION algebras
      mxMatrix("Full", values = (matrix(1, n.latent^2, n.latent^2) - diag(n.latent^2)), name = "tempb") #used in continuous DIFFUSION algebras  
    )
    
    if(objective!='Kalman') model<-OpenMx::mxModel(model, #include RAM matrices
      
      mxMatrix(values = A$values, free = A$free, labels = A$labels, dimnames = list(FILTERnamesy, FILTERnamesy), name = "A"),   #directed effect matrix   
      
      mxMatrix(values = S$values, free = S$free, labels = S$labels, dimnames = list(FILTERnamesy, FILTERnamesy), name = "S"),   #symmetric effect matrix
      
      mxMatrix(values = FILTER$values, free = FALSE, dimnames = list(FILTERnamesx, FILTERnamesy), name = "F"),  #filter matrix
      
      mxMatrix(free = t(M$free), values = t(M$values), labels = t(M$labels), dimnames = list(1, FILTERnamesy), name = "M") #mean matrix
    )
    
    
    ###Trait variance 
    if(traitExtension==TRUE){
      
      if('T0TRAITEFFECT' %in% stationary){
        T0TRAITEFFECT$labels <-paste0('T0TRAITEFFECTalg[', 1:n.latent, ',', rep(1:n.latent, each=n.latent), ']')
        T0TRAITEFFECT$free <-FALSE
        T0TRAITEFFECTalg<- OpenMx::mxAlgebra(name='T0TRAITEFFECTalg', -invDRIFT)     
      }
      
      model <- OpenMx::mxModel(model, 
        #       mxAlgebra(name = "traitloadings", (invDRIFT %*% (discreteDRIFT - II))), #trait loading matrix
        mxMatrix(type = "Full", labels = TRAITVAR$labels, values = TRAITVAR$values, byrow = TRUE, 
          free = TRAITVAR$free, lbound = if(reasonable==TRUE){zerodiagnlatent}else{NA}, name = "TRAITVAR"), 
        
        mxMatrix(name="T0TRAITEFFECT", labels=T0TRAITEFFECT$labels, #T0TRAITEFFECT matrix
          values=T0TRAITEFFECT$values, free=T0TRAITEFFECT$free, lbound=if(reasonable==TRUE){zerodiagnlatent}else{NA}, 
          nrow=n.latent, ncol=n.latent), 
        traitalgs #include trait RAM algebras     
      )
      
      if('T0TRAITEFFECT' %in% stationary) model<-OpenMx::mxModel(model, T0TRAITEFFECTalg)
      
    }
    
    ###Manifest Trait variance 
    if(manifestTraitvarExtension==TRUE){
      model <- OpenMx::mxModel(model, 
        mxMatrix(type = "Full", 
          labels = MANIFESTTRAITVAR$labels, 
          values = MANIFESTTRAITVAR$values, byrow = TRUE, 
          free = MANIFESTTRAITVAR$free, 
          lbound = if(reasonable==TRUE){zerodiagnmanifest}else{NA}, 
          name = "MANIFESTTRAITVAR")    
      )
    }
    
    
    
    
    
    
    
    #model options
    if(optimizer==TRUE) {
      message("Setting NPSOL optimizer for OpenMx temporarily") 
      mxOption(NULL, "Default optimizer", "NPSOL")
    }
    if(optimizer==FALSE) {
      message("Setting CSOLNP optimizer for OpenMx temporarily") 
      mxOption(NULL, "Default optimizer", "CSOLNP")
    }
    
    #     model <- mxOption(model, "Standard Errors", "No")
    #     model <- mxOption(model, "Calculate Hessian", "No")
    #     model <- mxOption(model, "No Sort Data", "ctsem")
    
    #     mxOption(model, "Derivative level", 0) #0
    #           mxOption(model, "Function precision", 1e-35) #1e-14
    #       #     mxOption(model, "Infinite bound size", 1e+15) #1.0e+15
#     mxOption(model, "Line search tolerance", .99) #.3
    #           mxOption(model, "Feasibility tolerance", 1.0e-14) #1.0e-05
    #       #     mxOption(model, "mvnMaxPointsA", 0) #0
    #       #     mxOption(model, "mvnMaxPointsB", 0) #0
    #           mxOption(model, "mvnMaxPointsC", 2500000) #5000
    #           mxOption(model, "mvnAbsEps", .000000001) #.001
    #           mxOption(model, "mvnRelEps", 0.0000000001) #0
    #     mxOption(model, "Verify level", 3) #-1
    #     mxOption(model, "Minor print level", 5) #-1
    #     mxOption(model, "Print level", 5) #-1
    #     mxOption(model, "Print file", "test") #-1
    
 #end base model spec
  
  
  
  
  
  
  
  
  ###predictor extension model specification
  if(n.TDpred + n.TIpred > 0){
    
    if(n.TDpred > 0 ) { #if there are fixed TD predictors 
            TDpredlbounds<-matrix(NA, nrow=(Tpoints-1)*n.TDpred, ncol=(Tpoints-1)*n.TDpred)
            if(reasonable==TRUE) diag(TDpredlbounds)<-.000001
      
      if('T0TDPREDCOV' %in% stationary){
        T0TDPREDCOV$values<-0
        T0TDPREDCOV$free<-FALSE     
      }
      
      model <- OpenMx::mxModel(model, 
        TDPREDEFFECTalgs, 
        mxMatrix(type = "Full", 
          labels = TDPREDEFFECT$labels, values = TDPREDEFFECT$values, free = TDPREDEFFECT$free, name = "TDPREDEFFECT"))
      
      if(objective!='Kalman'){
        
        model <- OpenMx::mxModel(model, 
          mxMatrix(type = "Full", 
            labels = TDPREDVAR$labels, values = TDPREDVAR$values, free = TDPREDVAR$free, name = "TDPREDVAR", 
            lbound = TDpredlbounds), 
          
          mxMatrix(type = "Full", 
            labels = T0TDPREDCOV$labels, values = T0TDPREDCOV$values, free = T0TDPREDCOV$free, name = "T0TDPREDCOV")
        )
        
        if(traitExtension==TRUE) model<-OpenMx::mxModel(model,        
          mxMatrix(labels = TRAITTDPREDCOV$labels, values = TRAITTDPREDCOV$values, free = TRAITTDPREDCOV$free, name = "TRAITTDPREDCOV")
        )
      }
    }
    
    
    if (n.TIpred > 0){ #if there are fixed time independent predictors
      
      if('T0TIPREDEFFECT' %in% stationary){
        T0TIPREDEFFECT$labels <-paste0('asymTIPREDEFFECT[', 1:n.latent, ',', rep(1:n.latent, each=n.latent), ']')
        T0TIPREDEFFECT$free<-FALSE
        asymTIPREDEFFECTalg<- OpenMx::mxAlgebra(name='asymTIPREDEFFECT', -invDRIFT %*% TIPREDEFFECT)  
        model<-OpenMx::mxModel(model, asymTIPREDEFFECTalg)
      }
      
      model <- OpenMx::mxModel(model, 
        mxMatrix(type = "Full", labels = TIPREDEFFECT$labels, values = TIPREDEFFECT$values, free = TIPREDEFFECT$free, name = "TIPREDEFFECT"), 
        mxMatrix(name='T0TIPREDEFFECT', values = T0TIPREDEFFECT$values, nrow=nrow(T0TIPREDEFFECT$values), ncol=ncol(T0TIPREDEFFECT$values), 
          free=T0TIPREDEFFECT$free, labels=T0TIPREDEFFECT$labels), 
        TIPREDEFFECTalgs
      )
    }
    
    
    
    
    #     if(randomPredictors==FALSE & n.TDpred > 0){ #if we want to insert TD predictors via mean inputs controlled by definition variables
    #       
    #       model <- OpenMx::mxModel(model, #add TDPREDEFFECT matrix to model
    #         mxMatrix(type = "Full", 
    #           labels = TDPREDEFFECTlabels, 
    #           values = TDPREDEFFECT$values, 
    #           free = TDPREDEFFECT$free, 
    #           name = "TDPREDEFFECT"), 
    # #         T0TDPREDCOVmatrices, 
    #         TDpreddefs, TDpreddefsfull, INTalgs 
    #         ) #add TDpreddefs and intalgs and T0cov to model
    #       
    #       
    #             model <- OpenMx::mxModel(model, remove=TRUE, 'T0MEANS')        
    #             model <- OpenMx::mxModel(model, 
    #               mxAlgebra(name = "T0MEANS", T0MEANBASE + T0TDPREDCOV %*% t(TDpreddefsfull)), 
    #               mxMatrix(name='T0TDPREDCOV', nrow=nrow(T0TDPREDCOV), ncol=ncol(T0TDPREDCOV), 
    #                 free=T0TDPREDCOV$free, labels=T0TDPREDCOVlabels), 
    #               mxMatrix(name = "T0MEANBASE", free=T0MEANS$free, labels=T0MEANS$labels, 
    #                 values=T0MEANS$valuesnrow=n.latent, ncol=1)
    #             )
    #     }# end fixed TD predictors
    #     
    #     
    #     
    #     
    #     if(randomPredictors==FALSE & n.TIpred >0){ #if we have TI predictors
    #       
    #       model <- OpenMx::mxModel(model, remove=TRUE, 'T0MEANS')        
    #       model <- OpenMx::mxModel(model, 
    #         mxMatrix(type = "Full", labels = TIPREDEFFECTlabels, values = TIPREDEFFECT$values free = TIPREDEFFECT$free, name = "TIPREDEFFECT"), 
    #         INTalgs, 
    #         mxMatrix(name="TIpreddefs", labels=paste0("data.", TIpredNames), nrow=1, ncol=n.TIpred), #create a definition variable matrix
    #         mxAlgebra(name = "T0MEANS", T0MEANBASE + T0TIPREDEFFECT %*% t(TIpreddefs)), 
    #         mxMatrix(name='T0TIPREDEFFECT', nrow=nrow(T0TIPREDEFFECT), ncol=ncol(T0TIPREDEFFECT), 
    #           free=T0TIPREDEFFECT$free, labels=T0TIPREDEFFECTlabels), 
    #         mxMatrix(name = "T0MEANBASE", free=T0MEANS$free, labels=T0MEANS$labels
    #           values=T0MEANS$valuesnrow=n.latent, ncol=1)
    #       )
    #     } # end fixed TI predictor / intercept algebra
    
#     returnAllLocals()
  } # end predictors model section
  
  
  
  
  
  
  
#   setobjective<-function(){
    ######### objective functions
    
    if(objective == "mxRAM") {  
      model <- OpenMx::mxModel(model, 
        mxAlgebra(F%*%solve(bigI - A)%*%S%*%t(solve(bigI - A))%*%t(F), name = "expCov"), 
        mxAlgebra(t(F%*%(solve(bigI - A))%*%t(M)), name = "expMean"), 
        mxMatrix(type = "Iden", nrow = nrow(A$labels), ncol = ncol(A$labels), name = "bigI"), 
        mxExpectationRAM(M = "M"), 
        mxFitFunctionML(vector=FALSE)
    )
    }
    
    
    if(objective == "mxFIML") { #split RAM style matrices into components for faster processing
      
      manifestExtent<-manifestend #update manifest extents - perhaps requires predictors
      latentExtent<-latentend
      if(n.TDpred + n.TIpred > 0) manifestExtent <- manifestExtent + n.TDpred *(Tpoints-1) + n.TIpred
      if(n.TDpred + n.TIpred > 0) latentExtent <- latentExtent + n.TDpred *(Tpoints-1) + n.TIpred
      if(traitExtension==TRUE) latentExtent<-latentExtent+n.latent
      if(manifestTraitvarExtension==TRUE) latentExtent <- latentExtent + n.manifest
      
      latents<-1:(n.latent*Tpoints)
      if(traitExtension==TRUE) latents<-c(latents,(n.latent*Tpoints+1):(n.latent*Tpoints+n.latent))
      if(manifestTraitvarExtension==TRUE) latents<-c(latents,(manifesttraitstart):(traitend))
      if(n.TDpred > 0) latents<-c(latents,predictorTDstart:predictorTDend)
      if(n.TIpred >0) latents <- c(latents,predictorTIstart:predictorTIend)
    
      Amanifestvalues <- A$values[manifeststart:manifestExtent, latents]
      Smanifestvalues<-S$values[manifeststart:manifestExtent, manifeststart:manifestExtent]
      Smanifestlabels<-S$labels[manifeststart:manifestExtent, manifeststart:manifestExtent]
      Smanifestfree<-S$free[manifeststart:manifestExtent, manifeststart:manifestExtent]
      
      Amanifestcovvalues <- A$values[manifeststart:manifestExtent, latents] #need this to incorporate predictor and manifest covariance
      if(n.TDpred+n.TIpred > 0) {
        Amanifestcovvalues[(n.manifest*Tpoints+1):(n.manifest*Tpoints+n.TIpred+n.TDpred*(Tpoints-1)),
        (traitend+1):(traitend+n.TIpred+n.TDpred*(Tpoints-1))] [col(diag(n.TIpred+n.TDpred*(Tpoints-1))) == 
            row(diag(n.TIpred+n.TDpred*(Tpoints-1)))] <- 1
        
        Smanifestvalues[(n.manifest*Tpoints+1):(n.manifest*Tpoints+n.TIpred+n.TDpred*(Tpoints-1)),
          (n.manifest*Tpoints+1):(n.manifest*Tpoints+n.TIpred+n.TDpred*(Tpoints-1))]  <-0
        Smanifestlabels[(n.manifest*Tpoints+1):(n.manifest*Tpoints+n.TIpred+n.TDpred*(Tpoints-1)),
          (n.manifest*Tpoints+1):(n.manifest*Tpoints+n.TIpred+n.TDpred*(Tpoints-1))]  <-NA
        Smanifestfree[(n.manifest*Tpoints+1):(n.manifest*Tpoints+n.TIpred+n.TDpred*(Tpoints-1)),
          (n.manifest*Tpoints+1):(n.manifest*Tpoints+n.TIpred+n.TDpred*(Tpoints-1))]  <-F
      }
  
      #base model components
      model <- OpenMx::mxModel(model, 
        
        mxMatrix(name='Alatent', values=A$values[latents, latents], 
          labels=A$labels[latents, latents], 
          free=A$free[latents, latents], 
          nrow=latentExtent, ncol=latentExtent), 
        
        mxMatrix(name='Slatent', values=S$values[latents, latents], 
          labels=S$labels[latents, latents], 
          free=S$free[latents, latents], 
          nrow=latentExtent, ncol=latentExtent), 
        
        mxMatrix(name='Mlatent', values=M$values[latents, 1], 
          labels=M$labels[latents, 1], 
          free=M$free[latents, 1], 
          nrow=1, ncol=latentExtent), 
        
        mxMatrix(name='Smanifest', labels=Smanifestlabels, #includes predictors!!
          values=Smanifestvalues, 
          free=Smanifestfree, 
          nrow=n.manifest*Tpoints+n.TDpred*(Tpoints-1)+n.TIpred, ncol=n.manifest*Tpoints+n.TDpred*(Tpoints-1)+n.TIpred), 
#           nrow=n.manifest*Tpoints, ncol=n.manifest*Tpoints), 
        
        mxMatrix(name='Mmanifest', labels=M$labels[manifeststart:manifestExtent], 
          values=M$values[manifeststart:manifestExtent], free=M$free[manifeststart:manifestExtent], 
          nrow=1, ncol=n.manifest*Tpoints+n.TDpred*(Tpoints-1)+n.TIpred), 
#           nrow=1, ncol=n.manifest*Tpoints), 
        
        mxMatrix(name='Amanifest', values=Amanifestvalues, 
          labels=A$labels[manifeststart:manifestExtent, latents], 
          free=A$free[manifeststart:manifestExtent, latents], 
          nrow=n.manifest*Tpoints+n.TDpred*(Tpoints-1)+n.TIpred, ncol=latentExtent), 
#           nrow=n.manifest*Tpoints, ncol=latentExtent), 
        
        mxMatrix(name='Amanifestcov', values=Amanifestcovvalues, 
          labels=A$labels[manifeststart:manifestExtent, latents], 
          free=A$free[manifeststart:manifestExtent, latents], 
          nrow=n.manifest*Tpoints+n.TDpred*(Tpoints-1)+n.TIpred, ncol=latentExtent), 
        
        mxMatrix(name='Ilatent', type='Iden', nrow=latentExtent, ncol=latentExtent), 
        
        mxAlgebra(name='invIminusAlatent', solve(Ilatent - Alatent))
      )
      
      fullSlatentlist<-'Slatent'   #begin a list of entities that sum together to generate the S matrix (ie traits, predictors, add to this)
      fullSmanifestlist<-'Smanifest'
      fullMlatentlist<-'Mlatent'
      
#      if(traitExtension==TRUE){
#        fullSlatentlist<-c(fullSlatentlist, 'traitCov')
#        model <- OpenMx::mxModel(model, 
#          
#          mxMatrix(name='Atrait', values=A$values[1:latentend, traitstart:latenttraitend], 
#            labels=A$labels[1:latentend, traitstart:latenttraitend], 
#            free=A$free[1:latentend, traitstart:latenttraitend], 
#            nrow=latentend, ncol=n.latent), 
#          
#          mxMatrix(name='Strait', values=S$values[traitstart:latenttraitend, traitstart:latenttraitend], 
#            labels=S$labels[traitstart:latenttraitend, traitstart:latenttraitend], 
#            free=S$free[traitstart:latenttraitend, traitstart:latenttraitend], 
#            nrow=n.latent, ncol=n.latent),       
#       
#          mxAlgebra(name='traitCov', Atrait %*%  Strait %*%  t(Atrait))  #latent trait covariance            
#        )
#      }
      
#       if(manifestTraitvarExtension==TRUE){
#         fullSmanifestlist<-c(fullSmanifestlist, 'manifesttraitCov')
#         model <- OpenMx::mxModel(model, 
#           
#           mxMatrix(name='Amanifesttrait', values=A$values[manifeststart:manifestExtent, manifesttraitstart:traitend], 
#             labels=A$labels[manifeststart:manifestExtent, manifesttraitstart:traitend], 
#             free=A$free[manifeststart:manifestExtent, manifesttraitstart:traitend], 
#             nrow=Tpoints*n.manifest+n.TDpred*(Tpoints-1)+n.TIpred, ncol=n.manifest), 
#           
#           mxMatrix(name='Smanifesttrait', values=S$values[manifesttraitstart:traitend, manifesttraitstart:traitend], 
#             labels=S$labels[manifesttraitstart:traitend, manifesttraitstart:traitend], 
#             free=S$free[manifesttraitstart:traitend, manifesttraitstart:traitend], 
#             nrow=n.manifest, ncol=n.manifest),       
#           
#           mxAlgebra(name='manifesttraitCov', Amanifesttrait %*%  Smanifesttrait %*%  t(Amanifesttrait))  #updated S manifest covariance with manifest traits
#         )
#       }
      
      
#       if(n.TDpred + n.TIpred > 0){
# #         fullSlatentlist<-c(fullSlatentlist, 'predCov')
# #         fullMlatentlist<-c(fullMlatentlist, 'predMean')
# #         browser()
#         model <- OpenMx::mxModel(model, 
#           
#           mxMatrix(name='Apred', values=A$values[1:latentend, predictorstart:predictorend], 
#             labels=A$labels[1:latentend, predictorstart:predictorend], 
#             free=A$free[1:latentend, predictorstart:predictorend], 
#             nrow=latentend, ncol=n.TDpred*(Tpoints-1)+n.TIpred), 
#           
#           mxMatrix(name='Spred', values=S$values[predictorstart:predictorend, predictorstart:predictorend], 
#             labels=S$labels[predictorstart:predictorend, predictorstart:predictorend], 
#             free=S$free[predictorstart:predictorend, predictorstart:predictorend], 
#             nrow=n.TDpred*(Tpoints-1)+n.TIpred, ncol=n.TDpred*(Tpoints-1)+n.TIpred),       
#           
# #           mxMatrix(name='Mpred', values=M$values[predictorstart:predictorend], 
# #             labels=M$labels[predictorstart:predictorend], 
# #             free=M$free[predictorstart:predictorend], 
# #             nrow=1, ncol=n.TDpred*(Tpoints-1) + n.TIpred) 
#           
#           mxMatrix(name='predmultiplier',nrow=latentExtent,ncol=n.TDpred*(Tpoints-1)+n.TIpred,
#             values=c(rep(0,latentend),rep(1,n.TDpred*(Tpoints-1)+n.TIpred)),free=F),
#           
#           mxAlgebra(name='Smanifestpredsmall', Amanifest %*% invIminusAlatent %*% fullSlatent %*% predmultiplier),
#           mxAlgebra(name='Smanifestpred',
#           
#           
#           
# #           mxAlgebra(name='predMean', Mpred %*% t(Apred)), 
# #           mxAlgebra(name='predCov', Apred %*%  Spred %*%  t(Apred))  #predictor covariance 
#           
#      
#         
# 
#         )
#       }
      
      fullSlatent<-eval(parse(text=paste0("mxAlgebra(name='fullSlatent', ", paste(fullSlatentlist, collapse=' + '), ")")))
      fullSmanifest<-eval(parse(text=paste0("mxAlgebra(name='fullSmanifest', ", paste(fullSmanifestlist, collapse=' + '), ")")))
      fullMlatent<-eval(parse(text=paste0("mxAlgebra(name='fullMlatent', ", paste(fullMlatentlist, collapse=' + '), ")")))
        
      model <- OpenMx::mxModel(model, 
        fullSlatent, fullSmanifest, fullMlatent, 
#         mxAlgebra(Amanifest %*% solve(bigI - Alatent) %*% Slatent %*% t(solve(bigI - Alatent)) %*% t(Amanifest) + Smanifest, name = "expCov"), 
#         mxAlgebra(Amanifest %*% invIminusAlatent %*% fullSlatent %*% t(invIminusAlatent) %*% t(Amanifest) + fullSmanifest + Smanifestpred, name = "expCov"), 
        mxAlgebra(Amanifestcov %*% invIminusAlatent %*% fullSlatent %*% t(invIminusAlatent) %*% t(Amanifestcov) +Smanifest, name = "expCov"), 
        mxAlgebra(t(Amanifest %*% (invIminusAlatent %*% t(fullMlatent))) + Mmanifest, name = "expMean"), 
#         mxMatrix(type = "Iden", nrow = nrow(A$labels, ncol = ncol(A$labels, name = "bigI"), 
        mxExpectationNormal(covariance = "expCov", means = "expMean", dimnames = FILTERnamesx), 
        mxFitFunctionML()
      )}
    
    if(objective == "means") {
      model <- OpenMx::mxModel(model, 
        mxAlgebra(F%*%solve(bigI - A)%*%S%*%t(solve(bigI - A))%*%t(F), name = "expCov"), 
        mxAlgebra(t(F%*%(solve(bigI - A)) %*% t(M)), name = "expMean"), 
        mxMatrix(type = "Iden", nrow = nrow(A$labels), ncol = ncol(A$labels), name = "bigI"), 
        #         mxExpectationNormal(covariance = "expCov", means = "expMean", dimnames = FILTERnamesx), 
        #         mxFitFunctionML(vector=FALSE)
        #         mxAlgebra( sum((expMean - datarow) * (expMean - datarow)), name='rowAlgebra'), 
        #         mxAlgebra( log2pi +log(det(expCov))+ (datarow - expMean) %*% solve(expCov) %*% t(datarow-expMean), name='rowAlgebra'), 
        mxAlgebra( (expMean - datarow) %*% t(expMean - datarow), name='rowAlgebra'), 
        
        mxAlgebra(expression=sum(existenceVector), name="numVar_i"), 
        mxMatrix("Full", 1, 1, values = log(2*pi), name = "log2pi"), 
        #         mxAlgebra(expression=omxSelectRowsAndCols(expCov, existenceVector), 
        #           name="filteredExpCov"), 
        #         mxAlgebra(solve(filteredExpCov), name='invFilteredExpCov'), 
        #         mxAlgebra(expression=log2pi + log(det(filteredExpCov)), 
        #           name ="firstHalfCalc"), 
        #         mxAlgebra(expression=(filteredDataRow - filteredExpMean) %*% solve(filteredExpCov) %*% t(filteredDataRow - filteredExpMean), 
        #           name = "secondHalfCalc"), 
        #         mxAlgebra(expression=(firstHalfCalc + secondHalfCalc), 
        #           name="rowAlgebra"), 
        
        
        mxAlgebra(omxSelectCols(expMean, existenceVector), name="filteredExpMean"), 
          
        mxMatrix(name='datarow', nrow=1, ncol=n.manifest*Tpoints, values=.5, 
          labels=paste0('data.', manifestNames, '_T', rep(0:(Tpoints-1), each=n.manifest)), free=F), 
          
        mxAlgebra(sum(rowResults), name='reduceAlgebra'),
          
        mxFitFunctionRow( rowAlgebra='rowAlgebra', 
          reduceAlgebra='reduceAlgebra', 
          dimnames=paste0(manifestNames, '_T', rep(0:(Tpoints-1), each=n.manifest)))   
        
      )}
    
    
    if(objective=='Kalman'){ 
      model<-OpenMx::mxModel(model, 
        
        mxMatrix(name='D', values=MANIFESTMEANS$values, labels=MANIFESTMEANS$labels, 
          nrow=n.manifest, ncol=1, free=MANIFESTMEANS$free), 
        
        mxMatrix(name='u', values=1, nrow=1, ncol=1, free=F), 
        
        mxMatrix(name='Qd', labels=paste0('Qd1[', 1:n.latent^2, ',1]'), 
          values=DIFFUSION$inits, nrow=n.latent, ncol=n.latent, free=F), 
        
        #         mxAlgebra(name='intd', t(intd1)), 
        
        mxExpectationStateSpace(A='discreteDRIFT_T1', B='intd1', C='LAMBDA', 
          D="D", Q='Qd', R='MANIFESTVAR', x0='T0MEANS', P0='T0VAR', u="u"), 
        
        mxFitFunctionML(vector=FALSE)
      )
      
      if(n.TDpred>0){
        intd1labels<-matrix(paste0('intd1[', 1:n.latent, ',','1]'), nrow=n.latent)
        TDPREDEFFECT_T1labels<-matrix(paste0('TDPREDEFFECT_T1[', 1:n.latent, ',', rep(1:n.TDpred, each=n.latent), ']'), nrow=n.latent)
        TDPRED$ref<-paste0('data.', TDpredNames)
        
        model<-OpenMx::mxModel(model, 
          mxMatrix(name='B', free=FALSE , nrow=n.latent, ncol=n.TDpred+1, 
            labels=cbind(intd1labels, TDPREDEFFECT_T1labels)), 
          
          mxMatrix(name='D', nrow=n.manifest, ncol=1+n.TDpred, 
            free=c(MANIFESTMEANS$free, rep(FALSE, n.TDpred*n.manifest)), 
            values=c(MANIFESTMEANS$values, rep(0, n.TDpred*n.manifest)), 
            labels=c(MANIFESTMEANS$labels, rep(NA, n.TDpred*n.manifest))), 
          
          mxMatrix(name='u', ncol=1, nrow=n.TDpred+1, free=FALSE, 
            values=c(1, rep(0, n.TDpred)), 
            labels=c(NA, TDPRED$ref)),         
          
          mxExpectationStateSpace(A='discreteDRIFT_T1', B='B', C='LAMBDA', D="D", Q='Qd', R='MANIFESTVAR', x0='T0MEANS', P0='T0VAR', u="u")
        )
      }
    }
    
    
    
    if(carefulFit==TRUE) {
      originalmodel<-model
      
      if(traitExtension==TRUE) penalties <- OpenMx::mxAlgebra(name='penalties', 
        sum(T0VAR*T0VAR) + sum(DRIFT*DRIFT) + sum(DIFFUSION*DIFFUSION) + sum(MANIFESTVAR*MANIFESTVAR) +
          sum(TRAITVAR * TRAITVAR))
      
      if(manifestTraitvarExtension==TRUE & traitExtension==FALSE) penalties <- OpenMx::mxAlgebra(name='penalties', 
        sum(T0VAR*T0VAR) + sum(DRIFT*DRIFT) + sum(DIFFUSION*DIFFUSION) + sum(MANIFESTVAR*MANIFESTVAR) +
          sum(MANIFESTTRAITVAR * MANIFESTTRAITVAR) )
      
      if(traitExtension==FALSE & manifestTraitvarExtension==FALSE) penalties <- OpenMx::mxAlgebra(name='penalties', 
        sum(T0VAR*T0VAR) + sum(DRIFT*DRIFT) + sum(DIFFUSION*DIFFUSION) + sum(MANIFESTVAR*MANIFESTVAR) )
      
      modelwithpenalties <- OpenMx::mxModel(model, 
        #             mxExpectationNormal(covariance = "expCov", means = "expMean", dimnames = FILTERnamesx), 
        #             mxAlgebra(F%*%solve(bigI - A)%*%S%*%t(solve(bigI - A))%*%t(F), name = "expCov"), 
        #             mxAlgebra(t(F%*%(solve(bigI - A))%*%t(M)), name = "expMean"), 
        #             mxMatrix(type = "Iden", nrow = nrow(Alabels), ncol = ncol(Alabels), name = "bigI"), 
        penalties, 
        mxFitFunctionML(vector=FALSE)
      )
      
      model<-OpenMx::mxModel('ctsemCarefulFit', 
        modelwithpenalties, 
        #             mxMatrix(type='Full', name='FIMLpenaltyweight', nrow=1, ncol=1, values=FIMLpenaltyweight, free=F), 
        mxAlgebra(name='FIMLpenaltyweight', ctsem.fitfunction / 1000 ), 
        mxAlgebra(sum(ctsem.fitfunction)+ctsem.penalties*FIMLpenaltyweight, name='m2ll'), #
        mxFitFunctionAlgebra('m2ll')
      )
    }
    
    
    #     if(objective == "mxFIMLpenalised") {  ## attempt for multigroup
    #       
    #       if(traitExtension==TRUE) penalties <- mxAlgebra(name='penalties', sum(DRIFT*DRIFT + TRAITVAR*TRAITVAR + 
    #           MANIFESTVAR*MANIFESTVAR + tr(DRIFT)/2) * FIMLpenaltyweight)
    #       if(traitExtension==FALSE) penalties <- mxAlgebra(name='penalties', sum(DRIFT*DRIFT + 
    #           MANIFESTVAR*MANIFESTVAR + tr(DRIFT)/2) * FIMLpenaltyweight)
    #       
    #       
    #       model <- OpenMx::mxModel(model, 
    #         #                    fitmodel, 
    #         mxExpectationNormal(covariance = "expCov", means = "expMean", dimnames = FILTERnamesx), 
    #         mxAlgebra(F%*%solve(bigI - A)%*%S%*%t(solve(bigI - A))%*%t(F), name = "expCov"), 
    #         mxAlgebra(t(F%*%(solve(bigI - A))%*%t(M)), name = "expMean"), 
    #         mxMatrix(type = "Iden", nrow = nrow(Alabels), ncol = ncol(Alabels), name = "bigI"), 
    #         #                         mxData(observed = datawide, type = "raw"), 
    #         mxMatrix(type='Full', name='FIMLpenaltyweight', nrow=1, ncol=1, values=FIMLpenaltyweight, free=F), 
    #         mxFitFunctionML(vector=TRUE), 
    #         penalties, 
    #         mxAlgebra(sum(fitfunction)+penalties, name='m2ll'), #
    #         mxFitFunctionAlgebra('m2ll')
    #       )
    #     }
    
    
    #       
    #       if(objective == "mxRowFIML") {
    #         model <- addmxrowobjective(model, dimlabels = FILTERnamesx)
    #         model <- OpenMx::mxModel(model, 
    #           mxRAMObjective("A", "S", "F", "M"), 
    #           mxAlgebra(F%*%solve(bigI - A)%*%S%*%t(solve(bigI - A))%*%t(F), name = "expCov"), 
    #           mxAlgebra(t(F%*%(solve(bigI - A))%*%t(M)), name = "expMean"), 
    #           mxMatrix(type = "Iden", nrow = nrow(Alabels), ncol = ncol(Alabels), name = "bigI")
    #         )}
    #       
    #       if(objective == "CondLogLik") {
    #         model <- OpenMx::mxModel(model, type = "default", mxFitFunctionR(CondLogLik), 
    #           mxAlgebra(F%*%solve(bigI - A)%*%S%*%t(solve(bigI - A))%*%t(F), name = "expCov"), 
    #           mxAlgebra(t(F%*%(solve(bigI - A))%*%t(M)), name = "expMean"), 
    #           mxMatrix(type = "Iden", nrow = nrow(Alabels), ncol = ncol(Alabels), name = "bigI")
    #         )}
    
    
    
#     returnAllLocals()
#   } #end setobjective function
  
  
  

  
  
  #   objective <- 'mxRAM' #only option for the moment so not included in initial arguments
  #   if(carefulFit==TRUE){
  #     targetObjective <- objective
  #     objective <- 'mxFIMLpenalised'
  #     
  #   }
  
  
  model <- OpenMx::omxAssignFirstParameters(model, indep = FALSE) #randomly selects inits for labels with 2 or more starting values
  if(showInits==TRUE & carefulFit==TRUE) { #
    message('Starting values: ')
    newstarts <- OpenMx::omxGetParameters(model)
    message(paste(names(newstarts), ": ", newstarts, "\n"))
  }
  
  if(plotOptimization==T){
    model<- OpenMx::mxOption(model, 'Always Checkpoint', 'Yes')
    model<- OpenMx::mxOption(model, 'Checkpoint Units', 'iterations')
    model<- OpenMx::mxOption(model, 'Checkpoint Count', 1)    
  }
  
  ###optimization options
  if(mxOption(NULL, "Default optimizer")=='CSOLNP') originaloptimizer<-'NPSOL'
  if(mxOption(NULL, "Default optimizer")=='NPSOL') originaloptimizer<-'CSOLNP'
  #   model<- mxOption(model, "mvnMaxPointsC", 150000) #default 5000, but 150000 improved solution finding...
  
  ###fit model
  
  if(carefulFit==TRUE) {
    carefulFit<-FALSE
    #     browser()
    mxobj<-try(suppressWarnings(OpenMx::mxRun(model))) #fit with the penalised likelihood
#         mxobj<-OpenMx::mxRun(model) #fit with the penalised likelihood
    newstarts <- try(OpenMx::omxGetParameters(mxobj)) #get the params
    if(showInits==TRUE) {
      message('Generated start values from carefulFit=TRUE')
      message(paste(names(newstarts), ": ", newstarts, "\n"))
    }
    model<-originalmodel #revert to our single layer model without the penalties fit function
    #     model<-OpenMx::mxModel(model, 'penalties', remove=TRUE) #and remove the penalties object
    if(class(newstarts)!="try-error") model<-OpenMx::omxSetParameters(model, labels=names(newstarts), values=newstarts) #set the params of it
    #     objective<-targetObjective #revert our objective to whatever was specified
    #     setobjective() #and set it
  }
  
  if(nofit == TRUE) return(model) #if we're not fitting the model, just return the unfitted openmx model
  
  if(nofit == FALSE){ #but otherwise...  
    
    if(!is.null(confidenceintervals))  {
      #     engine<-ifelse(npsol==TRUE, 'NPSOL', 'CSOLNP')
      #     mxOption(NULL, "Default optimizer", "NPSOL")
      model <- OpenMx::mxModel(model, 
        mxCI(confidenceintervals, interval = 0.95, type = "both")) #if 95% confidence intervals are to be calculated
      #     model <- OpenMx::mxModel(model, 
      #       mxComputeSequence(steps=list(
      #         mxComputeGradientDescent(engine=engine), 
      #         mxComputeConfidenceInterval(engine="NPSOL"), 
      # #         MxComputeNumericDeriv(), 
      #         MxComputeStandardError()
      # #         MxComputeHessianQuality(), 
      # #         MxComputeReportDeriv()        
      #       )))
      
    }
    
    #     if(retryattempts > 0){ 
    mxobj <- ctmxTryHard(model, 
      fit2beat = fit2beat, showInits=showInits, checkHess=FALSE, greenOK=TRUE, 
      #         intervals = ifelse(!is.null(confidenceintervals), TRUE, FALSE), 
      #         confidenceintervals=confidenceintervals, 
      iterationSummary=iterationSummary, bestInitsOutput=FALSE, 
      extraTries=retryattempts, loc=1, scale=.3, paste=FALSE)
    #     }
    
    #     if(retryattempts < 1 ) {
    #       mxobj <- OpenMx::mxRun(model, useOptimizer = useOptimizer, 
    #         intervals = ifelse(!is.null(confidenceintervals), TRUE, FALSE))
    #     }   
  }
  
  
  if(plotOptimization==TRUE){
    
    checkpoints<-read.table(file='ctsem.omx', header=T, sep='\t')
    if(carefulFit==TRUE) checkpoints<-rbind(checkpoints, NA, NA, NA, NA, NA, read.table(file='ctsemCarefulfit.omx', header=T, sep='\t'))
    mfrow<-par()$mfrow
    par(mfrow=c(3, 3))
    for(i in 6:ncol(checkpoints)) {
      plot(checkpoints[, i], main=colnames(checkpoints)[i])
    }
    par(mfrow=mfrow)
    deleteCheckpoints <- readline('Remove created checkpoint file, ctsem.omx? y/n \n')
    if(deleteCheckpoints=='y') file.remove(file='ctsem.omx')
  }
  
  mxOption(NULL, "Default optimizer", originaloptimizer) #reset optimizer
  
  if(discreteTimeModel==FALSE) {
    out <- list(mxobj, ctmodelobj, ctfitargs, datawide) #roll unfitted and fitted model and ctmodelobj into one list item
    names(out) <- c("mxobj", "ctmodelobj", "ctfitargs", "datawide")
    class(out) <- "ctsemFit" #and give it the class of a ctsemFit object
  }
  if(discreteTimeModel==TRUE) out <- mxobj

  return(out)
}