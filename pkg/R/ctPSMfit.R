#' ctPSMfit
#' 
#' Fits a specified ctsem model using the PSM package. 
#'
#' @param datawide a dataset in the wide format used by ctsem.
#' @param ctmodelobj A ctsem model specified using \code{\link{ctModel}}. As between subject effects and predictors are 
#' not implemented currently, ensure the model contains none of these or the fit may fail.
#'
#' @return PSM fit data
#' @export

ctPSMfit<-function(datawide,ctmodelobj){
  
  if (!requireNamespace("PSM", quietly = TRUE)) {
    stop("PSM package needed for this function to work. Please install it.",
      call. = FALSE)
  }
  
  psmdata<-lapply(1:nrow(datawide), function(x) {
    
    long<-ctWideToLong(datawide[x,,drop=F], Tpoints=ctmodelobj$Tpoints, 
      n.manifest=ctmodelobj$n.manifest, manifestNames=ctmodelobj$manifestNames)
    
    long<-suppressMessages(ctDeintervalise(long))
    colnames(long)[ncol(long)] <-'Time'
    out<-list(Time=long[,'Time'], Y=t(long[,ctmodelobj$manifestNames]))
    out$U=matrix(1,nrow=1,ncol=ncol(out$Y))
    return(out)
  })
  
  
  flatcm<-unlist(ctmodelobj[-which(names(ctmodelobj) %in% c('latentNames', 'manifestNames', 'T0VAR'))])
  flatcm<-flatcm[-which(!is.na(as.numeric(flatcm)))]
  names(flatcm)<-flatcm
  inits<-unlist(lapply(flatcm, function(x) rnorm(1,.3,.01)))
  
  MyModel <- vector(mode="list")
  MyModel$Matrices=function(phi) {
    
    matA<-ctmodelobj$DRIFT

    matA[which(matA %in% names(phi))] <- 
      as.numeric(unlist(phi[match(matA[which(matA %in% names(phi))],names(phi))]))
    mode(matA)<-'numeric'
    diag(matA)<- -exp(diag(matA))
    
    matB<-ctmodelobj$CINT
    matB[which(matB %in% names(phi))] <- 
      as.numeric(unlist(phi[match(matB[which(matB%in% names(phi))],names(phi))]))
    mode(matB)<-'numeric'
    
    matC<-ctmodelobj$LAMBDA
    matC[which(matC %in% names(phi))] <- 
      as.numeric(unlist(phi[match(matC[which(matC %in% names(phi))],names(phi))]))
    mode(matC)<-'numeric'
    
    matD<-ctmodelobj$MANIFESTMEANS
    matD[which(matD %in% names(phi))] <- 
      as.numeric(unlist(phi[match(matD[which(matD%in% names(phi))],names(phi))]))
    mode(matD)<-'numeric'
    
    out<-list(matA=matA,
      matB=matB,
      matC=matC,
      matD=matD
    )
    
    out
  }
  
  MyModel$h = function(eta,theta,covar) {
    
    
    if(!is.null(ctmodelobj$TRAITVAR)){
      
      names(theta)<-flatcm[-which(flatcm %in% ctmodelobj$TRAITVAR )]
      
      matB<-ctmodelobj$CINT
      matB[which(matB %in% names(theta))] <- 
        as.numeric(unlist(theta[match(matB[which(matB %in% names(theta))],names(theta))]))
      mode(matB)<-'numeric'
      
      matB<-matB+eta

      phi <- theta
      phi[which(names(phi) %in% ctmodelobj$CINT)] <- matB #this assignment may not be reliable - relies on ordering of phi matching matB
      return(phi)
    }
    
    
    if(is.null(ctmodelobj$TRAITVAR)){
    phi <- theta
    return(phi)
    }
    
  }
  
  MyModel$S = function(phi) {
    S<-ctmodelobj$MANIFESTVAR
    
    S[which(S %in% names(phi))] <- 
      as.numeric(unlist(phi[match(S[which(S %in% names(phi))],names(phi))]))
    mode(S)<-'numeric'
    diag(S)<-exp(diag(S))
    S<-t(S) %*% S
    
  }
  MyModel$SIG = function(phi) {
    SIG<-ctmodelobj$DIFFUSION
    
    SIG[which(SIG %in% names(phi))] <- 
      as.numeric(unlist(phi[match(SIG[which(SIG %in% names(phi))],names(phi))]))
    mode(SIG)<-'numeric'
    diag(SIG)<-exp(diag(SIG))
    SIG<-t(SIG) %*% SIG
  }
  MyModel$X0 = function(Time,phi,U) {
    X0<-ctmodelobj$T0MEANS
    
    X0[which(X0 %in% names(phi))] <- 
      as.numeric(unlist(phi[match(X0[which(X0 %in% names(phi))],names(phi))]))
    mode(X0)<-'numeric'
    X0
  }
  MyModel$ModelPar = function(THETA) {
    
    names(THETA)<-flatcm
    
    if( !is.null(ctmodelobj$TRAITVAR)){
     
    OMEGA<-ctmodelobj$TRAITVAR
    out<-list(theta=lapply(THETA[-which(names(THETA) %in% OMEGA )], function(x) x))
    names(out$theta)<-flatcm[-which(flatcm %in% OMEGA )]
    
    OMEGA[which(OMEGA %in% names(THETA))] <- 
      as.numeric(unlist(THETA[match(OMEGA[which(OMEGA %in% names(THETA))],names(THETA))]))
    mode(OMEGA)<-'numeric'
    diag(OMEGA)<-exp(diag(OMEGA))
    OMEGA<-t(OMEGA) %*% OMEGA
    
    out$OMEGA<-OMEGA
    
    }
    
    
    if( is.null(ctmodelobj$TRAITVAR)){
    out<-list(theta=lapply(THETA, function(x) x))
    names(out$theta)<-flatcm
    }
    
    out
  }
  
  
  fit<-PSM::PSM.estimate(MyModel, psmdata, Par=list(Init=inits), control=list(optimizer='ucminf', grad='forward'),fast=F)
  
  return(list(PSMfit=fit, PSMmodel=MyModel, PSMinitpars=list(Init=inits)))
}
