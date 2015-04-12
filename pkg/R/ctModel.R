#' Define a ctsem model
#' 
#' This function is used to specify a continuous time structural equation model, 
#' which can then be fit to data with function \code{\link{ctFit}}.  
#' @param n.manifest Number of manifest indicators per individual at each measurement occasion / time point.  
#' Manifest variables are included as the first element of the wide data matrix, with all the 1:n.manifest manifest variables 
#' at time 1 followed by those of time 2, and so on.
#' 
#' @param n.latent Number of latent processes.
#' 
#' @param Tpoints Number of time points, or measurement occasions, in the data.  This will generally be the maximum 
#' number of time points for a single individual, but may be one extra if sample relative time intervals are used, 
#' see \code{\link{ctIntervalise}}.  
#' 
#' @param LAMBDA n.manifest*n.latent loading matrix relating latent to manifest variables, 
#' with latent processes 1:n.latent along the columns, and manifest variables
#' 1:n.manifest in the rows.
#' 
#' @param manifestNames n.manifest length vector of manifest variable names as they appear in the data structure, 
#' without the _Tx time point suffix.  Defaults to Y1, Y2, etc.
#' 
#' @param latentNames n.latent length vector of latent variable names 
#' (used for naming parameters, defaults to eta1, eta2, etc).
#' 
#' @param T0VAR symmetric n.latent*n.latent matrix of latent process initial variance / covariance. 
#' "auto" freely estimates all parameters and generates starting values.
#' 
#' @param T0MEANS n.latent*1 matrix of latent process means at first time point, T0. 
#' "auto" freely estimates all parameters and generates starting values.
#' 
#' @param MANIFESTMEANS n.manifest*1 matrix of manifest means.  
#' "auto" fixes all parameters to 0.
#' 
#' @param MANIFESTVAR symmetric n.manifest*n.manifest matrix of variance / covariance 
#' between manifests at each measurement occasion (i.e. measurement error / residual).  
#' "auto" freely estimates variance parameters and generates starting values, 
#' and fixes covariances between manifests to 0.
#' 
#' @param DRIFT n.latent*n.latent DRIFT matrix of continuous auto and cross effects, 
#' relating the processes over time. 
#' "auto" freely estimates all parameters and generates starting values.
#' 
#' @param CINT n.latent * 1 matrix of continuous intercepts, allowing for non 0 
#' stationary levels of the processes.  
#' "auto" freely estimates all parameters and generates starting values.
#' 
#' @param DIFFUSION symmetric n.latent*n.latent matrix of diffusion process 
#' variance and covariance (generates latent error / innovation).
#' "auto" freely estimates all parameters and generates starting values.
#' 
#' @param TRAITVAR Either NULL, if no trait / individual heterogeneity effect, 
#' or symmetric n.latent*n.latent matrix of trait variance / covariance.
#' "auto" freely estimates all parameters and generates starting values.
#' 
#' @param T0TRAITEFFECT n.latent*n.latent matrix of effect of traits on latents at T0. 
#' If TRAITVAR is NULL (default), T0TRAITEFFECT is also set to NULL, 
#' and ignored in any case.
#' 
#' @param MANIFESTTRAITVAR either NULL (default) if no trait variance / individual heterogeneity in the level of
#' the manifest indicators, otherwise a symmetric n.manifest * n.manifest variance / covariance matrix. 
#' Set to "auto" to include and free all parameters.
#' 
#' @param n.TDpred Number of time dependent predictors in the dataset.  
#' Each time dependent predictor must have Tpoints-1 columns 
#' (A value of a predictor at the final time point would not have any time to show an effect), 
#' which are included in the data matrix after manifest variables,
#' with the Tpoints-1 observations of the first predictor followed by those of the second and so on.
#' 
#' @param TDpredNames n.TDpred length vector of time dependent predictor variable names, 
#' as they appear in the data structure, without the _Tx time point suffix.  
#' Default names are TD1, TD2, etc.
#' 
#' @param TDPREDMEANS (n.TDpred*(Tpoints-1)) rows * 1 column matrix of time dependent predictor means.
#' If 'auto', the means are freely estimated.  Otherwise, 
#' the means for the Tpoints-1 observations of your first time dependent predictor 
#' are followed by those of TDpred 2, and so on.
#' 
#' @param TDPREDEFFECT n.latent*n.TDpred matrix of effects from time dependent predictors to latent processes.
#' Effects from 1:n.TDpred columns TDpredictors to 1:n.latent rows of latent processes.
#' "auto" freely estimates all parameters and generates starting values.
#' 
#' @param T0TDPREDCOV n.latent rows * ((Tpoints-1) rows * n.TDpred) columns covariance matrix 
#' between latents at T0 and time dependent predictors.
#' "auto" freely estimates all parameters and generates starting values.
#' 
#' @param TDPREDVAR symmetric (n.TDpred * (Tpoints-1)) rows * (n.TDpred * (Tpoints-1)) columns variance / covariance
#' matrix for time dependent predictors.
#' "auto" (default) freely estimates all parameters and generates random starting values.
#' 
#' @param TRAITTDPREDCOV n.latent rows * (n.TDpred*(Tpoints-1)) columns covariance matrix of 
#' latent traits and time dependent predictors.  
#' The Tpoints-1 columns of the first preditor are followed by those of the second and so on.
#' Covariances with the trait variance of latent process 1 are specified in row 1, process 2 in row 2, etc.
#' "auto" (default) freely estimates covariance between time dependent predictors and traits, 
#' (if both exist, otherwise NULL).
#' 
#' @param TDTIPREDCOV (n.TDpred * (Tpoints-1)) rows * n.TIpred columns covariance
#' matrix between time dependent and time independent predictors.
#' "auto" (default) freely estimates all parameters and generates random starting values.
#' 
#' @param n.TIpred Number of time independent predictors. 
#' Each TIpredictor is inserted at the right of the data matrix, after the time intervals.
#' 
#' @param TIpredNames n.TIpred length vector of time independent predictor variable names,
#' as they appear in the data structure.  Default names are TI1, TI2, etc.
#'  
#' @param TIPREDMEANS n.TIpred * 1 matrix of time independent predictor means.
#' If 'auto', the means are freely estimated.  
#' 
#' @param TIPREDEFFECT n.latent*n.TIpred effect matrix of time independent predictors on latent processes.
#' "auto" freely estimates all parameters and generates starting values.
#' 
#' @param T0TIPREDEFFECT n.latent*n.TIpred effect matrix of time independent 
#' predictors on latents at T0.
#' "constrained" replicates labels from TIPREDEFFECT, forcing the time 0 effect to be the same as the long run effect.
#' This can be freed by setting to 'auto' or manually specifying, but may make optimization more difficult.
#' 
#' @param TIPREDVAR symmetric n.TIpred * n.TIpred variance / covariance
#' matrix for all time independent predictors.
#' "auto" (default) freely estimates all parameters and generates random starting values.
#' 
#' @param startValues a named vector, where the names of each value must match a parameter in the specified model,
#' and the value sets the starting value for that parameter during optimization.
#' If not set, random starting values representing relatively stable processes with small effects and 
#' covariances are generated by ctFit.  
#' Better starting values may improve model fit speed and the chance of an appropriate model fit.  
#' 
#' @param initsToFixed if TRUE, uses the provided list of starting values as fixed values instead of starting values
#' (convenient way to extract previous fit parameters using omxGetParameters and constrain a new model)
#' @export

ctModel<-function(n.manifest, n.latent, Tpoints, LAMBDA, 
  manifestNames='auto', latentNames='auto', 
  T0VAR="auto", T0MEANS="auto", MANIFESTMEANS="auto", MANIFESTVAR="auto", 
  DRIFT="auto", CINT="auto", DIFFUSION="auto",
  TRAITVAR=NULL, T0TRAITEFFECT="auto", MANIFESTTRAITVAR=NULL,
  n.TDpred=0, TDpredNames='auto', TDPREDMEANS="auto", TDPREDEFFECT="auto", 
  T0TDPREDCOV="auto", TDPREDVAR="auto",  
  TRAITTDPREDCOV="auto", TDTIPREDCOV='auto',
  n.TIpred=0, TIpredNames='auto', TIPREDMEANS="auto", TIPREDEFFECT="auto", 
  T0TIPREDEFFECT="auto", TIPREDVAR="auto", 
  startValues=NULL, initsToFixed=FALSE){
  
  
  
  #####FUNCTIONS
  
  
  isx<-function(x){ #if that returns FALSE for both FALSE and NA results
    if(is.na(x)){
      return(FALSE)
      break()
    }
    if(x==FALSE){
      return(FALSE)
      break()
    }
    if(x==TRUE){
      return(TRUE)
    }}
  
  checkSymmetry<-function(x){  #this checks the symmetry of matrix x, if not symmetric it stops and outputs an error
    if(isSymmetric(x)==F){
      stop(paste0(substitute(x)," must be symmetric"))
    }
  }
  


  
  

  
  
  ###### RUN SEQUENCE
  
  
  #names
  if(all(manifestNames=='auto')) manifestNames=paste0('Y',1:n.manifest)
  if(all(latentNames=='auto')) latentNames=paste0('eta',1:n.latent)
  
  if(length(manifestNames) != n.manifest) stop("Length of manifestNames does not equal n.manifest!") 
  if(length(latentNames) != n.latent) stop("Length of latentNames does not equal n.latent!") 
  
  if(n.TDpred > 0){
    if(all(TDpredNames=='auto')) TDpredNames=paste0('TD',1:n.TDpred)
  if(length(TDpredNames) != n.TDpred) stop("Length of TDpredNames does not equal n.TDpred!") 
  }
  
  if(n.TIpred > 0){
  if(all(TIpredNames=='auto')) TIpredNames=paste0('TI',1:n.TIpred)
  if(length(TIpredNames) != n.TIpred) stop("Length of TIpredNames does not equal n.TIpred!") 
  }
  
  
#matrices
  if(T0MEANS[1]=="auto") T0MEANS<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,
    manifestNames=manifestNames,latentNames=latentNames,matrixname="T0MEANS",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  
  if(ncol(T0MEANS)>1) stop("Specified T0MEANS matrix has more than one column")
  if(nrow(T0MEANS)!=n.latent) stop ("Specified T0MEANS matrix rows not equal to n.latent")
  
  
  
  
  
  if(T0VAR[1]=="auto") T0VAR<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,
    manifestNames=manifestNames,latentNames=latentNames,matrixname="T0VAR",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  
  if(nrow(T0VAR)!=n.latent | ncol(T0VAR)!=n.latent) stop("Dimensions of T0VAR matrix are not n.latent * n.latent")
  
  
  
  
  
  if(MANIFESTMEANS[1]=="auto") MANIFESTMEANS<- matrix(0,nrow=n.manifest)
  if(ncol(MANIFESTMEANS)>1) stop("Specified MANIFESTMEANS matrix has more than one column")
  if(nrow(MANIFESTMEANS)!=n.manifest) stop ("Specified MANIFESTMEANS matrix rows not equal to n.manifest")
  
  
  
  
  
  if(MANIFESTVAR[1]=="auto") MANIFESTVAR<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,
    manifestNames=manifestNames,latentNames=latentNames,matrixname="MANIFESTVAR",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  
  MANIFESTVAR[row(MANIFESTVAR) != col(MANIFESTVAR)] <-0
  
  if(nrow(MANIFESTVAR)!=n.manifest | ncol(MANIFESTVAR)!=n.manifest) stop("Dimensions of MANIFESTVAR matrix are not n.manifest * n.manifest")
  
  
  
  
  
  if(DRIFT[1]=="auto") DRIFT<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="DRIFT",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  if(nrow(DRIFT)!= n.latent | ncol(DRIFT)!= n.latent) stop("Dimensions of DRIFT matrix are not n.latent * n.latent")
  
  
  
  
  
  if(CINT[1]=="auto") CINT<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,manifestNames=manifestNames,
    latentNames=latentNames,matrixname="CINT",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  if(ncol(CINT)!=1 | nrow(CINT)!=n.latent) stop("Dimensions of CINT are not 1 * n.latent")
  
  
  
  
  
  if(DIFFUSION[1]=="auto") DIFFUSION<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="DIFFUSION",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  if(nrow(DIFFUSION)!=n.latent | ncol(DIFFUSION) != n.latent) stop("Dimensions of DIFFUSION are not n.latent * n.latent")

  
  
  
  
  
  if(!is.null(TRAITVAR) && TRAITVAR[1]=="auto") TRAITVAR<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="TRAITVAR",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  if(!is.null(TRAITVAR) && nrow(TRAITVAR)+ncol(TRAITVAR)!=n.latent*2) stop ("Dimensions of TRAITVAR are not n.latent * n.latent")
  
  
  
  
  
  if(!is.null(TRAITVAR) && T0TRAITEFFECT[1]=="constrained")  T0TRAITEFFECT<-diag(n.latent)
  if(!is.null(TRAITVAR) && T0TRAITEFFECT[1]=="auto")  T0TRAITEFFECT<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="T0TRAITEFFECT",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  if(is.null(TRAITVAR) && T0TRAITEFFECT[1]=="auto") T0TRAITEFFECT<-NULL  
  if(!is.null(TRAITVAR)) if(nrow(T0TRAITEFFECT)!= n.latent | ncol(T0TRAITEFFECT) != n.latent) stop("Dimensions of T0TRAITEFFECT are not n.latent * n.latent")
  
  
  
  
  
  if(!is.null(MANIFESTTRAITVAR) && MANIFESTTRAITVAR[1]=="auto") MANIFESTTRAITVAR<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="MANIFESTTRAITVAR",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  if(!is.null(MANIFESTTRAITVAR) && nrow(MANIFESTTRAITVAR)+ncol(MANIFESTTRAITVAR)!=n.manifest*2) stop ("Dimensions of MANIFESTTRAITVAR are not n.manifest * n.manifest")
#   if(is.null(MANIFESTTRAITVAR)) MANIFESTTRAITVAR <- matrix(0,nrow=n.manifest,ncol=n.manifest)
  
  
  
  
  
  
  #effect matrix of time dependent predictors (with time independent effect strength)
  if(TDPREDEFFECT[1]=="auto" && n.TDpred!=0)  TDPREDEFFECT<-ctLabel(TDpredNames=TDpredNames,
    TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="TDPREDEFFECT",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  if(n.TDpred==0) TDPREDEFFECT<-NULL
  if(n.TDpred!=0) if(nrow(TDPREDEFFECT) != n.latent | ncol(TDPREDEFFECT) != n.TDpred) stop("Dimensions of TDPREDEFFECT are not n.latent * n.TDpred")
  
  
  

  # means of TD predictors
  if(TDPREDMEANS[1]=="auto" & n.TDpred > 0) TDPREDMEANS <- ctLabel(TDpredNames=TDpredNames,
    TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="TDPREDMEANS",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  
  if(n.TDpred==0) TDPREDMEANS<-NULL
  
  if(n.TDpred > 0) if(nrow(TDPREDMEANS) != n.TDpred*(Tpoints-1) | ncol(TDPREDMEANS) != 1) stop(
    "Dimensions of TDPREDMEANS are not (n.TDpred*(Tpoints-1)) rows * 1 column")
  
  
  
  
  
  # means of TI predictors
  if(TIPREDMEANS[1]=="auto" & n.TIpred >0) TIPREDMEANS <- ctLabel(TDpredNames=TDpredNames,
    TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="TIPREDMEANS",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  
  if(n.TIpred==0) TIPREDMEANS<-NULL
  
  if(n.TIpred > 0) if(nrow(TIPREDMEANS) != n.TIpred | ncol(TIPREDMEANS) != 1) stop(
    "Dimensions of TIPREDMEANS are not n.TIpred * 1")
  
  
  
  
  if(T0TDPREDCOV[1]=="auto" && n.TDpred>0) T0TDPREDCOV<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="T0TDPREDCOV",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  if(n.TDpred==0) T0TDPREDCOV<-NULL

  if(n.TDpred>0) if(nrow(T0TDPREDCOV) != n.latent | ncol(T0TDPREDCOV) != n.TDpred*(Tpoints-1))  stop(
    "Dimensions of T0TDPREDCOV are not n.latent * (n.TDpred*(Tpoints-1))")
  
  
  
  
  
  if(TIPREDEFFECT[1]=="auto" && n.TIpred>0) TIPREDEFFECT<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="TIPREDEFFECT",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  if(n.TIpred==0) TIPREDEFFECT<-NULL
    if(n.TIpred>0) if(nrow(TIPREDEFFECT) != n.latent | ncol(TIPREDEFFECT) != n.TIpred) stop("Dimensions of TIPREDEFFECT are not n.latent * n.TIpred")

  
  
  
  if(T0TIPREDEFFECT[1]=="constrained" && n.TIpred > 0) T0TIPREDEFFECT<-TIPREDEFFECT
  if(T0TIPREDEFFECT[1]=="auto" && n.TIpred > 0) T0TIPREDEFFECT <- ctLabel(TDpredNames=TDpredNames,
    TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="T0TIPREDEFFECT",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  
  if(n.TIpred==0) T0TIPREDEFFECT <- NULL
  
  if(n.TIpred>0) if(nrow(T0TIPREDEFFECT) != n.latent | ncol(T0TIPREDEFFECT) != n.TIpred) stop("Dimensions of T0TIPREDEFFECT are not n.latent * n.TIpred")

  
  
  
  
  if(TIPREDVAR[1]=="auto" && n.TIpred > 0)  TIPREDVAR <- ctLabel(TDpredNames=TDpredNames,
    TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="TIPREDVAR",n.latent=n.latent,
      n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints) 
  
  if(n.TIpred==0) TIPREDVAR <- NULL
  
    if(n.TIpred > 0) if(nrow(TIPREDVAR) != n.TIpred | ncol(TIPREDVAR) != n.TIpred) stop(
      "Dimensions of TIPREDVAR are not n.TIpred * n.TIpred")
  
  
  
  
  
  if(all(TDPREDVAR=="auto") && n.TDpred >0)  TDPREDVAR <- ctLabel(TDpredNames=TDpredNames,
    TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="TDPREDVAR",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)  
  
  if(n.TDpred==0) TDPREDVAR<-NULL
  
  if(n.TDpred > 0) if(nrow(TDPREDVAR) != n.TDpred*(Tpoints-1) | ncol(TDPREDVAR) != n.TDpred*(Tpoints-1)) stop(
    "Dimensions of TDPREDVAR are not (n.TDpred*(Tpoints-1)) * (n.TDpred*(Tpoints-1))")
  
  
  
  
  if(TDTIPREDCOV=="auto" & n.TDpred > 0 & n.TIpred > 0)  TDTIPREDCOV <- ctLabel(TDpredNames=TDpredNames,
    TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="TDTIPREDCOV",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  
   
  if(n.TDpred == 0 | n.TIpred == 0) TDTIPREDCOV<-NULL
  
  if(n.TDpred > 0 & n.TIpred > 0) if(nrow(TDTIPREDCOV) != (n.TDpred*(Tpoints-1)) | ncol(TDTIPREDCOV) != n.TIpred) stop(
    "Dimensions of TDTIPREDCOV are not (n.TDpred*(Tpoints-1)) * n.TIpred")
  
  
  
  
  
  # TDpredictor and trait covariances
  if(n.TDpred > 0) { #then set covariance matrix for heterogeneity and predictors    
    if(TRAITTDPREDCOV[1] == "auto") TRAITTDPREDCOV<-ctLabel(TDpredNames=TDpredNames,
      TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="TRAITTDPREDCOV",n.latent=n.latent,
      n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)  
    
    if(nrow(TRAITTDPREDCOV) != n.latent | ncol(TRAITTDPREDCOV) != n.TDpred*(Tpoints-1)) stop(
      "Dimensions of TRAITTDPREDCOV are not n.latent * (n.TDpred*(Tpoints-1))")
  }
 
  if(n.TDpred==0)  TRAITTDPREDCOV<-NULL #if no predictors, assume no covariance between heterogeneity and predictors
  
  
  
  ###implement model checks here
  checkSymmetry(DIFFUSION)
  checkSymmetry(MANIFESTVAR)
  checkSymmetry(T0VAR)
 
  if(any(dim(LAMBDA)!=c(n.manifest,n.latent))) stop("Incorrect LAMBDA structure specified - check number or rows and columns")
  if(!is.null(TRAITVAR[1])) { checkSymmetry(TRAITVAR) }
  
  #check no phit1traits without traitvars
  if(isx(is.null(TRAITVAR))){
    if(isx(!all(TRAITVAR==0))){
      if(isx(any(!is.na(T0TRAITEFFECT)))){
        if(isx(any(T0TRAITEFFECT>0))){
          stop("Latent and Trait covariances specified (T0TRAITEFFECT) with no trait variance (TRAITVAR)")
        }}}}
  
  
  

  
  
  completemodel<-list(n.manifest,n.latent,n.TDpred,n.TIpred,Tpoints,LAMBDA,
    manifestNames,latentNames,TDpredNames,TIpredNames,
    T0VAR,T0MEANS,MANIFESTMEANS,MANIFESTVAR,DRIFT,CINT,DIFFUSION,
    TRAITVAR,T0TRAITEFFECT,MANIFESTTRAITVAR,
    TDPREDEFFECT, TDPREDMEANS, T0TDPREDCOV, TRAITTDPREDCOV,
    TIPREDEFFECT, TIPREDMEANS, T0TIPREDEFFECT,
    TIPREDVAR, TDPREDVAR,TDTIPREDCOV,   
    startValues)
  
  names(completemodel)<-c("n.manifest","n.latent","n.TDpred","n.TIpred","Tpoints","LAMBDA",
    "manifestNames","latentNames","TDpredNames","TIpredNames",
    "T0VAR","T0MEANS","MANIFESTMEANS","MANIFESTVAR","DRIFT","CINT","DIFFUSION",
    "TRAITVAR","T0TRAITEFFECT","MANIFESTTRAITVAR",
    'TDPREDEFFECT', 'TDPREDMEANS', 'T0TDPREDCOV', 'TRAITTDPREDCOV',
    'TIPREDEFFECT', 'TIPREDMEANS', 'T0TIPREDEFFECT',
    'TIPREDVAR', 'TDPREDVAR','TDTIPREDCOV',  
    "startValues")
  
  
  
  if(initsToFixed==TRUE){ #if we want to take init params as fixed instead - can be a simple way to vary models
    
    for( i in 1:length(completemodel)){
      
      if(is.matrix(completemodel[[i]])){
        
        fixed<-completemodel[[i]]

        
        if(!is.null(startValues)){ #if there are some startValues specified, check each part of x
          for( j in 1:length(startValues)){
            
              fixed[grepl(names(startValues)[j], fixed)] <- startValues[j]
          }
        }
        
        completemodel[[i]]<-fixed     
      }
  }
  }
  
  
  class(completemodel)<-"ctsemInit"
  return(completemodel)
}