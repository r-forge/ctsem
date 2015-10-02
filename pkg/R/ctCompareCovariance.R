#' ctCompareCovariance
#' Compares model implied to observed auto and cross correlations for panel data fit with ctsem. 
#' @param fitobj Fitted model object from OpenMx or ctsem.
#' @param outputmatrices if TRUE, output expected, observed, and residual correlation matrices as well as plots.
#' @param pause if TRUE (default) output plots interactively, one at a time.  If FALSE, output without stopping.
#' @param varlist if "all" include all variables in dataset.  Otherwise, specify numeric vector of variables to include.
#' @param ylim vector of min and max Y axis limits for plot.
#' @param ... additional arguments passed to plot. 
#' @export
ctCompareCovariance<-function(fitobj,
  outputmatrices=FALSE,pause=TRUE,varlist="all",ylim=c(-1,1),...){ 

  if(class(fitobj)=="ctsemFit"){ 
    n.manifest<-fitobj$ctmodelobj$n.manifest
    Tpoints<-fitobj$ctmodelobj$Tpoints
    mxobj<-fitobj$mxobj    
  }
  
  if(class(fitobj)=="MxModel" &
      any( c(!exists('n.manifest'), !exists('Tpoints')))) message(
        'OpenMx object specified - additional arguments required: n.manifest, Tpoints') 
  

  if(varlist=="all") varlist <- 1:n.manifest
  
  datawide<-mxobj@data@observed
  if(dim(datawide)[1] == dim(datawide)[2]) original<-cov2cor(datawide) else {
    datawide<-datawide[,paste0(fitobj$ctmodelobj$manifestNames, '_T',rep(0:(Tpoints-1),each=n.manifest))]
    original<-round(cor(datawide,use="pairwise.complete.obs"),digits=3) #subobtimal
  } 

  modelexp<-cov2cor(mxobj$fitfunction$info$expCov[1:(n.manifest*Tpoints),1:(n.manifest*Tpoints)])

  residuals <-  original-modelexp

  out<-list(original,modelexp,residuals)
  names(out)<-c("original","expected","residual")
  


  for(i in varlist){
    for(j in varlist){
    
      suppressWarnings(plot(rep(seq(n.manifest,Tpoints*n.manifest,n.manifest)/n.manifest,each=Tpoints)+ #x section
          seq(-.3,.3,length.out=Tpoints),#placement within x section
        (out$original[cbind(rep(seq(j,Tpoints*n.manifest,n.manifest),times=Tpoints),
          rep(seq(i,Tpoints*n.manifest,n.manifest),each=Tpoints))]),
        main=paste0("Cor variable ",i,"*",j),
        ylab="Correlation",xlab="Measurement occasion",col="blue",pch=16,ylim=c(-1,1),...))
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