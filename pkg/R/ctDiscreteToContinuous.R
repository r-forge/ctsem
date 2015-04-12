#'ctDiscreteToContinuous
#'Converts values in ctModel matrices from discrete time to continuous time parameters.
#'@param ctmodelobj model specified via ctModel function. Must contain only fixed numeric values, 
#'parameter values should represent those from a discrete time model.
#'@param timeInterval time interval the provided discrete parameters represent. 1 may be appropriate.
#'@details Does not convert T0TRAITEFFECT or T0TIPREDEFFECT matrices yet.
#'@export
ctDiscreteToContinuous <- function(ctmodelobj,timeInterval){
  
  
  for(i in 1:length(ctmodelobj)){ #convert matrices to numeric
    if(is.matrix(ctmodelobj[[i]])){
      ctmodelobj[[i]] <- matrix(as.numeric(ctmodelobj[[i]]),nrow=nrow(ctmodelobj[[i]]))
    }
  }
  
  II<-diag(ctmodelobj$n.latent)
  
  ctmodelobj$DRIFT <- OpenMx::logm(ctmodelobj$DRIFT)
  
  DRIFTHATCH<-(ctmodelobj$DRIFT %x% diag(2) + diag(2) %x% ctmodelobj$DRIFT)
  
  ctmodelobj$DIFFUSION <- matrix(DRIFTHATCH %*% solve((OpenMx::expm(DRIFTHATCH * timeInterval)) - (II%x%II)) %*% 
    OpenMx::rvectorize(ctmodelobj$DIFFUSION),nrow(II))
  
  ctmodelobj$CINT <- ctmodelobj$DRIFT %*% solve(OpenMx::expm(ctmodelobj$DRIFT * timeInterval) - II) %*% 
      (ctmodelobj$CINT)
  
  smalltraitloadings <- solve(DRIFT) %*% (-OpenMx::expm(ctmodelobj$DRIFT * timeInterval)) - solve(DRIFT)
  
  if(!is.null(ctmodelobj$TRAITVAR)) ctmodelobj$TRAITVAR<-matrix(
    solve(smalltraitloadings %x% smalltraitloadings) %*% c(ctmodelobj$TRAITVAR),nrow=2)
  
  if(ctmodelobj$n.TIpred > 0) ctmodelobj$TIPREDEFFECT <- -ctmodelobj$DRIFT %*% 
    solve(diag(2) - OpenMx::expm(ctmodelobj$DRIFT * timeInterval) )    %*%  (ctmodelobj$TIPREDEFFECT)
  
  if(ctmodelobj$n.TDpred > 0) ctmodelobj$TDPREDEFFECT <- 
    solve(OpenMx::expm(ctmodelobj$DRIFT * timeInterval)) %*% ctmodelobj$TDPREDEFFECT

  return(ctmodelobj)
}