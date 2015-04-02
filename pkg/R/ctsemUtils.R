# helper function to generate an index matrix, or return unique elements of a matrix
indexMatrix<-function(dimension,symmetrical=FALSE,upper=FALSE,sep=NULL,starttext=NULL,endtext=NULL,
  unique=FALSE,rowoffset=0,coloffset=0,indices=FALSE,diagonal=TRUE,namesvector=NULL){
  if(is.null(namesvector)) namesvector=1:9999
  if(indices==T) sep<-c(",")
  tempmatrix<-matrix(paste0(starttext,namesvector[1:dimension+rowoffset],sep,rep(namesvector[1:dimension+coloffset],each=dimension),endtext),nrow=dimension,ncol=dimension)
  if(upper==TRUE) tempmatrix<-t(tempmatrix)
  if(symmetrical==TRUE) tempmatrix[col(tempmatrix)>row(tempmatrix)] <-t(tempmatrix)[col(tempmatrix)>row(tempmatrix)]
  if(unique==TRUE && symmetrical==TRUE) tempmatrix<-tempmatrix[lower.tri(tempmatrix,diag=diagonal)]
  if(indices==T){
    tempmatrix<-matrix(c(unlist(strsplit(tempmatrix,","))[seq(1,length(tempmatrix)*2,2)],
      unlist(strsplit(tempmatrix,","))[seq(2,length(tempmatrix)*2,2)]),ncol=2)
  }
  return(tempmatrix)
}

#plots means of ctsem wide panel data
meanplot<-function(data,n.manifest,Tpoints,...){
  plot(1:Tpoints,rep(0,times=Tpoints),ylim=c(min(colMeans(data[,1:(n.manifest*Tpoints)],na.rm=T)),
    max(colMeans(data[,1:(n.manifest*Tpoints)],na.rm=T))),...)
  for(i in 1:n.manifest){
    points(1:Tpoints,colMeans(data[,seq(i,n.manifest*Tpoints,n.manifest)],na.rm=T),col=(i+1),type="b")
  }}

#Allows easy creation of matrices with a diagonal of character vectors.

chardiag<-function(x,dim=x){
  out<-diag(dim)
  diag(out)<-x
  return(out)
}

#' ctGetInits
#' 
#' Extracts estimates from a fitted ctsem model and returns in ctsem init matrix layout
#' @param ctfitobject ctsem fit object to extract new starting values from
#' @export
#' @import OpenMx

ctGetInits<-function(ctfitobject){
  if(class(ctfitobject)!='ctsemFit') stop('Specified fitobject is not of class ctsemFit')
  inits<-matrix(omxGetParameters(ctfitobject$mxobj),ncol=1)
  inits<-cbind(names(omxGetParameters(ctfitobject$mxobj)),inits)
  return(inits)
}



#Remove observations at random from a wide ctsem format dataset
# @param retainpercent Percent of each row to retain. Must be between 1 and 100, retainpercent/100*Tpoints must generate a whole number.
ctRemoveObservations<-function(datawide,Tpoints,n.manifest,n.TDpred=0,n.TIpred=0,retainpercent=100,
  manifestNames="auto",TDpredNames="auto",TIpredNames='auto'){
  
  datalong<-ctWideToLong(datawide,
    Tpoints=Tpoints,
    n.manifest=n.manifest,
    n.TIpred=n.TIpred,
    n.TDpred=n.TDpred) #convert to long
  datalong<-ctDeintervalise(datalong,
    n.manifest=n.manifest,
    n.TIpred=n.TIpred,
    n.TDpred=n.TDpred) #convert to absolute time
  
  samplelist<-c()
  for(i in 1:nrow(datawide)){
    samplelist<-c(samplelist,
      sample(  ((i-1)*Tpoints+1) : 
          (i*Tpoints), 
        ceiling(retainpercent/100 * Tpoints))) #keep retainpercent rows of each individual at random
  }
  
  datalong<-datalong[samplelist,,drop=F] #retain only samplelist rows
  
  datalong<-datalong[order(datalong[,"subject"],datalong[,"AbsTime"]),] #fix ordering of rows by time
  
  if(manifestNames=='auto') manifestNames <- paste0("Y",1:n.manifest)
  if(n.TDpred >0 && TDpredNames=='auto') TDpredNames <- paste0("TD",1:n.TDpred) 
  if(n.TDpred ==0) TDpredNames <- NULL
  if(n.TIpred >0 && TIpredNames=='auto') TIpredNames <- paste0("TI",1:n.TIpred) 
  if(n.TIpred == 0) TIpredNames <- NULL
  
  datawide<-ctLongToWide(datalong,
    id="subject",
    manifestNames=manifestNames,
    TDpredNames=TDpredNames,
    TIpredNames=TIpredNames,
    time="AbsTime") 
  
  datawide<-ctIntervalise(datawide,
    manifestNames=manifestNames,
    TDpredNames=TDpredNames,
    TIpredNames=TIpredNames,
    Tpoints=Tpoints*retainpercent/100,
    n.manifest=n.manifest,
    n.TDpred=n.TDpred,
    n.TIpred=n.TIpred,
    mininterval=.001,
    individualRelativeTime=TRUE)
  rownames(datawide)<-1:nrow(datawide)
  return(datawide)
}



# sets default column names for wide ctsem datasets

ctWideNames<-function(n.manifest,n.TDpred=0,Tpoints,n.TIpred=0,manifestNames='auto',TDpredNames='auto',TIpredNames='auto'){
  
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
  
  manifestnames<-paste0(manifestNames,"_T",rep(0:(Tpoints-1),each=n.manifest))
  if(n.TDpred > 0) TDprednames<-paste0(TDpredNames,"_T",rep(0:(Tpoints-2),each=n.TDpred)) else TDprednames<-NULL
  intervalnames<-paste0("dT",1:(Tpoints-1))
  if(n.TIpred>0) TIprednames <- paste0(TIpredNames) else TIprednames <- NULL
  return(c(manifestnames,TDprednames,intervalnames,TIprednames))
}