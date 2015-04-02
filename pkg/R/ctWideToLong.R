#' Convert ctsem wide to long format
#' @param dataset ctsem wide format data
#' @param Tpoints number of measurement occasions in data
#' @param n.manifest number of manifest variables
#' @param n.TDpred number of time dependent predictors
#' @param n.TIpred number of time independent predictors
#' @param manifestNames Character vector of manifest variable names.
#' @param TDpredNames Character vector of time dependent predictor names.
#' @param TIpredNames Character vector of time independent predictor names.
#' @export

ctWideToLong<-function(dataset,Tpoints,n.manifest,n.TDpred=0, n.TIpred=0, 
  manifestNames='auto',TDpredNames='auto',TIpredNames='auto'){ 
  
  #names
  if(all(manifestNames=='auto')) manifestNames=paste0('Y',1:n.manifest)
  if(length(manifestNames) != n.manifest) stop("Length of manifestNames does not equal n.manifest!") 

  if(n.TDpred > 0){
    if(all(TDpredNames=='auto')) TDpredNames=paste0('TD',1:n.TDpred)
    if(length(TDpredNames) != n.TDpred) stop("Length of TDpredNames does not equal n.TDpred!") 
  }
  
  if(n.TIpred > 0){
    if(all(TIpredNames=='auto')) TIpredNames=paste0('TI',1:n.TIpred)
    if(length(TIpredNames) != n.TIpred) stop("Length of TIpredNames does not equal n.TIpred!") 
  }
  
  datawide<-as.matrix(dataset,nrow=nrow(dataset),ncol=ncol(dataset)) #set to matrix for manipulations
  n.subjects<-nrow(datawide) #calculate number of subjects in dataset
  
  
  datalong<-matrix(NA,nrow=n.subjects*Tpoints,ncol=1+n.manifest+1+n.TDpred+n.TIpred) #create blank datalong matrix
  colnames(datalong)<-paste0("name",1:ncol(datalong)) #create default colnames
  colnames(datalong)[1:(1+n.manifest)]<-c("subject",manifestNames) #set column names for manifests
  colnames(datalong)[ncol(datalong)-n.TIpred]<-"dT" #set colnames for time intervals
  
  if(n.TDpred>0)  colnames(datalong)[(1+n.manifest+1):(1+n.manifest+n.TDpred)]<-TDpredNames #set TDpredictor colnames    
  if(n.TIpred>0)  colnames(datalong)[(ncol(datalong)-(n.TIpred-1)):ncol(datalong)]<-TIpredNames #set TIpredictor colnames 
 
  for(j in 1:n.subjects){ #for each subject
    for(i in 1:Tpoints){ #at each time point
      
      datalong[(j-1)*Tpoints+i, 1:(1+n.manifest)]<-
        c(j, #add a subject indicator /id to datalong
          datawide[j,(i-1)*n.manifest+(1:n.manifest)]) #add manifests to datalong
      
      if(i==1) datalong[(j-1)*Tpoints+i,(ncol(datalong)-n.TIpred)]<-0 #if first obs, add interval of 0
      
      if(i!=1) datalong[(j-1)*Tpoints+i,(ncol(datalong)-n.TIpred)]<- #if not first obs, set time interval to
        datawide[j,Tpoints*n.manifest+(Tpoints-1)*n.TDpred+(i-1)] #observed interval
      
      if(n.TDpred>0) for(p in 1:n.TDpred){ #if time dependent predictors exist, for each:
        datalong[(j-1)*Tpoints+i,(1+n.manifest+p)]<-c(
          ifelse(i==Tpoints,NA,datawide[j,(n.manifest*Tpoints+(p-1)*(Tpoints-1)+(i))])) #add time dependent predictor to datalong if not last obs
      }
      
      if(n.TIpred>0) { #if time independent predictors exist
        datalong[(j-1)*Tpoints+i,(1+n.manifest+n.TDpred+1+1:n.TIpred)]<-
          datawide[j,(n.manifest*Tpoints+(Tpoints-1)+n.TDpred*(Tpoints-1)+1:n.TIpred),drop=FALSE] #add time independent predictor to all obs
      }        
    }
  }
  return(datalong)
}