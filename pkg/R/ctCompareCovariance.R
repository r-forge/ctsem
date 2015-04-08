#' ctCompareCovariance
#' Compares model implied to observed auto and cross correlations for panel data fit with ctsem. 
#' @param ctfitobj Fitted model object from OpenMx or ctsem.
#' @param outputmatrices if TRUE, output expected, observed, and residual correlation matrices as well as plots.
#' @param pause if TRUE (default) output plots interactively, one at a time.  If FALSE, output without stopping.
#' @param varlist if "all" include all variables in dataset.  Otherwise, specify numeric vector of variables to include.
#' @param ... additional arguments passed to plot. 
#' @export
ctCompareCovariance<-function(ctfitobj,
  outputmatrices=FALSE,pause=TRUE,varlist="all",...){ 

  checkOpenMx('ctCompareCovariance')
  
  if(class(ctfitobj)!="ctsemFit"){ 
    stop('Not a ctsemFit object')   
  }
  

  n.latent<-ctfitobj$ctmodelobj$n.latent
  n.TIpred<-ctfitobj$ctmodelobj$n.TIpred
  n.TDpred<-ctfitobj$ctmodelobj$n.TDpred
  n.manifest<-ctfitobj$ctmodelobj$n.manifest
  Tpoints<-ctfitobj$ctmodelobj$Tpoints
  
  if(varlist=="all") varlist <- 1:n.manifest

  datawide<-ctfitobj$mxobj$data$observed[,1:(n.manifest*Tpoints)]
  manifestindices<-n.latent*Tpoints+(1:(n.manifest*Tpoints))
#   modelexp<-cov2cor(ctfitobj$mxobj$expectation$UnfilteredExpCov[manifestindices,manifestindices])
  modelexp<-cov2cor(ctfitobj$mxobj$expCov$result[1:(n.manifest*Tpoints),1:(n.manifest*Tpoints)])
  
  OpenMx::mxEval(expCov,ctfitobj$mxobj,compute=T,defvar.row=nrow(datawide))

   original<-round(cor(datawide,use="pairwise.complete.obs"),digits=3) #subobtimal
  residuals <-  original-modelexp

  out<-list(original,modelexp,residuals)
  names(out)<-c("original","expected","residual")

  for(i in varlist){
    for(j in varlist){
    
      plot(rep(seq(n.manifest,Tpoints*n.manifest,n.manifest)/n.manifest,each=Tpoints)+ #x section
          seq(-.3,.3,length.out=Tpoints),#placement within x section
        (out$original[cbind(rep(seq(j,Tpoints*n.manifest,n.manifest),times=Tpoints),
          rep(seq(i,Tpoints*n.manifest,n.manifest),each=Tpoints))]),
        ylim=ylim, main=paste0("Cor variable ",i,"*",j),
        ylab="Correlation",xlab="Measurement occasion",col="blue",pch=16,...)
      abline(v=seq(1.5,Tpoints-.5,1),col="grey",lty=2)
      
      
      points(rep(seq(n.manifest,Tpoints*n.manifest,n.manifest)/n.manifest,each=Tpoints)+
          seq(-.3,.3,length.out=Tpoints),
        out$expected[cbind(rep(seq(j,Tpoints*n.manifest,n.manifest),times=Tpoints),
          rep(seq(i,Tpoints*n.manifest,n.manifest),each=Tpoints))],col="red")
      
      legend("bottomleft",legend=list("original","expected"),text.col=c("blue","red"),bty="n")
      
      if(i*j<n.manifest^2){
        if(pause==T){
          message("Press [enter] to display next graph")
          readline()
        }
      }
    }
  }
  if(outputmatrices==TRUE) return(out)
}