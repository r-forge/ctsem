#' Simulate continuous time data
#' 
#' This function generates data according to the specified ctsem model object, 
#' which must contain fixed values for parameters.
#' @param ctmodelobj ctsem model object from \code{\link{ctModel}}.
#' @param n.subjects Number of subjects to output.
#' @param burnin Number of initial time points to discard (to simulate stationary data)
#' @param dT Time interval (delta T) to use, defaults to 1.
#' @param TDpredtype Time dependent predictor type to generate data for, if specified.
#' @param asymptotes Are the parameters provided asymptotic paramters, or the regular continuous time parameters?
#' @export

ctGenerate<-function(ctmodelobj,n.subjects=1000,burnin=0,TDpredtype='impulse',dT=1,asymptotes=FALSE){
  
  ###read in model
  for(i in 1:length(ctmodelobj)){ #this loop reads in the specified continuous time model
    assign(names(ctmodelobj[i]),ctmodelobj[[i]])

    if(is.matrix(ctmodelobj[[i]])){ #if this element is a matrix, continue on...
     
    if(any(is.na(suppressWarnings(as.numeric(get(names(ctmodelobj[i]))))))){ #if it contains character labels
      assign(names(ctmodelobj[i]),matrix(0,nrow=nrow(get(names(ctmodelobj[i]))), #set the values to 0 instead
                                      ncol=ncol(get(names(ctmodelobj[i])))))
      message(paste0(names(ctmodelobj[i])," contained character labels - setting matrix to 0"))
    }
    
     #set any matrices to numeric elements
      assign(names(ctmodelobj[i]), 
        matrix(as.numeric(get(names(ctmodelobj[i]))),nrow=nrow(get(names(ctmodelobj[i]))),
                                      ncol=ncol(get(names(ctmodelobj[i])))))
    }
  }
  
  #set up extra matrices
  DRIFTHATCH <- DRIFT %x% diag(n.latent) + diag(n.latent) %x% DRIFT #generate drifthatch
  if(asymptotes==FALSE) dynresidualcov <- matrix(solve(DRIFTHATCH)%*%((OpenMx::expm(DRIFTHATCH %x% dT)) - #generate dynamic error cov from continuous value
                                                  diag(1,n.latent^2))%*%rvectorize(DIFFUSION),nrow=n.latent)
  if(asymptotes==TRUE) dynresidualcov <- matrix((diag(n.latent^2) - OpenMx::expm(DRIFT %x% dT) %x% OpenMx::expm(DRIFT %x% dT)) %*% c(DIFFUSION),nrow=n.latent)
  
  
  if(!all(is.numeric(T0VAR))) { #if T0VAR does not have all values fixed
    print("No T0VAR specified - generated process will be at equilibrium")
    T0VAR<-diag(100,n.latent)+1 #arbitrarily set it
    if(burnin < 200){ #and if burnin has not been set
      burnin <- 200 #ensure burnin is high enough that arbitrary phit1 doesn't matter
    }
  }
  
  if(!all(is.numeric(T0MEANS))) T0MEANS<-matrix(0,ncol=n.latent) #if T0MEANS has not been fixed arbitrarily set it

  
  
  
  
  ####traits
  traiteffect<-matrix(0,nrow=n.subjects,ncol=n.latent) #create traiteffect matrix with 0 effect
  if(!is.null(TRAITVAR[1])) { #if traits are specified
    if(!all(is.numeric(T0TRAITEFFECT))){#if not all of T0TRAITEFFECT is not fixed
      T0TRAITEFFECT<-diag(100,n.latent) #arbitrarily set it
      if(burnin==0) burnin<-200 #and ensure burnin is adequate (as for T0VAR)
    }
if(asymptotes==FALSE) traitloading<- solve(DRIFT) %*% (OpenMx::expm(DRIFT %x% dT) - diag(1,n.latent)) #calc trait loadings to latents for small(continuous) traitvar
if(asymptotes==TRUE) traitloading<-  diag(n.latent)- OpenMx::expm(DRIFT %x% dT)  #calc trait loadings to latents for asymptotic traitvar
    trait <- MASS::mvrnorm(n=n.subjects,mu=rep(0,n.latent),Sigma=TRAITVAR,tol=1) #generate trait effects    
    traiteffect<- trait %*% t(traitloading)
 
#     withinphi <- T0VAR - T0TRAITEFFECT %*% TRAITVAR %*% T0TRAITEFFECT #calculate non-trait related variance for T0
  }
  
  
  #####predictors
  
  TDpreds<-matrix(NA,nrow=n.subjects,ncol=n.TDpred*(Tpoints-1))
  TDpredeffects <- matrix(0,nrow=n.subjects,ncol=n.latent*(Tpoints-1)) #create 0 TDpredeffects in case no TDpreds
  if (n.TDpred>0) { #but if TDpreds exist
    if(TDpredtype=='level') TDpredparam <- solve(DRIFT) %*% (OpenMx::expm(DRIFT %x% dT) - diag(1, n.latent)) %*%  TDPREDEFFECT #calculate effect size
    if(TDpredtype=='impulse') TDpredparam <- OpenMx::expm(DRIFT %x% dT) %*%  TDPREDEFFECT #calculate effect size
    TDpreds <- MASS::mvrnorm(n=n.subjects,mu=TDPREDMEANS, #generate TDpred variables from TDPREDMEANS and TDPREDVAR
                       Sigma=TDPREDVAR ,tol=1)
    
    TDpreds <- TDpreds + traiteffect %*% TRAITTDPREDCOV
    
    for(l in 1:n.latent){ #for each latent process
      for(i in 1:(Tpoints-1)){#for each latent state variable except the final ones
        #     TDpredeffects[,(i-1)*n.latent+l] <-  TDpreds[,seq(i,n.TDpred*(Tpoints-1),(Tpoints-1))] %*% t(TDpredparam[l, ]) #create TDpred effects
        TDpredeffects[,(i-1)*n.latent+l] <-  TDpreds[,seq(i,n.TDpred*(Tpoints-1),(Tpoints-1))] %*% t(TDpredparam[l, ,drop=F]) #create TDpred effects
      }
    }
  }    
  
  TIpreds<-matrix(NA,nrow=n.subjects,ncol=n.TIpred)
  TIpredeffects<-matrix(0,nrow=n.subjects,ncol=n.latent) #set 0 effects in case no TIpreds
  if (n.TIpred>0) { #if TIpreds exist
    if(asymptotes==FALSE) TIPREDEFFECTdiscrete <-  solve(DRIFT) %*% (OpenMx::expm(DRIFT %x% dT) - diag(1, n.latent)) %*%  TIPREDEFFECT #calculte their effect
    if(asymptotes==TRUE) TIPREDEFFECTdiscrete <-  (diag(1, n.latent) - OpenMx::expm(DRIFT %x% dT) ) %*%  TIPREDEFFECT #calculate their effect
    
    TIpreds <- matrix(
      MASS::mvrnorm(n=n.subjects,
              mu=TIPREDMEANS, #generate TIpreds
              Sigma=TIPREDVAR, tol=1)
      ,nrow=n.subjects)
      #     TIpredeffects<-matrix(TIpreds*rep(TIpredparam,n.subjects*n.latent),nrow=n.subjects,ncol=n.TIpred) #generate effects on the latents
      TIpredeffects<-matrix(TIpreds %*% t(TIPREDEFFECTdiscrete^2),nrow=n.subjects,ncol=n.latent) #generate effects on the latents
  }
  
  ####cint
  
  cinteffect <- matrix(solve(DRIFT) %*% (OpenMx::expm(DRIFT %x% dT) - diag(1, n.latent)) %*% 
                         CINT,byrow=T,nrow=n.subjects,ncol=n.latent) #continuous intercept effect
  
  
  
  Tpoints<-Tpoints+burnin #add burnin to Tpoints (after we checked if extra burnin was needed for T0 cov, and after predictor generation)
  T0VAReffect<-MASS::mvrnorm(n=n.subjects,mu=rep(0,n.latent),Sigma=(T0VAR),tol=1) #create effect of non-trait variance at T0
  
  
  
  
  
  
  ########create latents  
  latents <- matrix(,nrow=n.subjects,ncol=n.latent*Tpoints) #create latent matrix
  latents[,1:n.latent] <- traiteffect + matrix(T0MEANS,nrow=n.subjects,ncol=n.latent,byrow=T)+ T0VAReffect #create T0 latents including all possible effects
  
  #burnin
  if(burnin>0){
    for(i in seq(n.latent+1,n.latent*(burnin+1),n.latent)){ #for every time block of latents after the first
      drifteffect <- latents[,(i-n.latent):(i-1)] %*% t(OpenMx::expm(DRIFT %x% dT))#effect of past time points    
      qeffect <- MASS::mvrnorm(n=n.subjects,mu=rep(0,n.latent),Sigma=dynresidualcov,tol=1) #effect of q noise
      
      latents[,i:(i+n.latent-1)]<-drifteffect+cinteffect+qeffect+traiteffect+TIpredeffects#sum of all constant effects at t to create latents
    }
    latents<-as.matrix(latents[,-1:-(burnin*n.latent),drop=FALSE]) #remove burnin from latents
    Tpoints<-Tpoints-burnin #remove burnin from Tpoints
  }
  
  
  #post burnin latents (here TDpredictors are added)
  for(i in seq(n.latent+1,n.latent*Tpoints,n.latent)){ #for every time block of latents after the first
    drifteffect <- latents[,(i-n.latent):(i-1)] %*% t(OpenMx::expm(DRIFT %x% dT))#effect of past time points    
    qeffect <- MASS::mvrnorm(n=n.subjects,mu=rep(0,n.latent),Sigma=dynresidualcov,tol=1) #effect of q noise
    latents[,i:(i+n.latent-1)]<-drifteffect+cinteffect+qeffect+traiteffect+TIpredeffects+
      TDpredeffects[,(i-n.latent) : (i-1)]  #sum of all effects at t to create latents
  }
  
  #generate indicators from latents

  manifests<-matrix(,nrow=n.subjects,ncol=n.manifest*Tpoints) #create latent matrix
  
  #   
  mantraiteffect<-matrix(0,nrow=n.subjects,ncol=n.manifest) #create traiteffect matrix with 0 effect
  if(!is.null(MANIFESTTRAITVAR[1])) { #if traits are specified
    mantraiteffect <- MASS::mvrnorm(n=n.subjects,mu=rep(0,n.manifest),Sigma=MANIFESTTRAITVAR,tol=1) #generate trait effects
  }
  
  #   
  for(i in 1:Tpoints){
    manifests[,((i-1)*n.manifest+1):(i*n.manifest)] <- latents[,((i-1)*n.latent+1):(i*n.latent)] %*% t(LAMBDA)
  }
  
  for(i in 1:n.manifest){
    manifests[,seq(i,Tpoints*n.manifest,n.manifest)]<-manifests[,seq(i,Tpoints*n.manifest,n.manifest),drop=FALSE]+
      rnorm(Tpoints*n.subjects,0,sqrt(MANIFESTVAR[i,i]))+ #measurement residual
      MANIFESTMEANS[i]+rep(mantraiteffect[,i],Tpoints)
  }
  
  intervals<-matrix(dT,nrow=n.subjects,ncol=(Tpoints-1)) #add intervals to latent output
  out<-as.matrix(cbind(manifests,TDpreds,intervals,TIpreds)) #output variables
  colnames(out)<-ctWideNames(Tpoints=Tpoints,n.manifest=n.manifest,n.TIpred=n.TIpred,n.TDpred=n.TDpred)
  rownames(out)<-1:nrow(out)
  return(out)
}
