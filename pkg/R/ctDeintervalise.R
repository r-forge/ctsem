#' ctDeintervalise
#' 
#' Converts intervals in ctsem long format data to absolute time
#' @param datalong data to use, in ctsem long format (attained via function ctWideToLong)
#' @param n.manifest Number of manifest variables
#' @param n.TDpred Number of time dependent predictor variables
#' @param n.TIpred Number of time independent predictor variables
#' @param startoffset Number of units of time to offset by when converting.
#' @export
ctDeintervalise<-function(datalong,n.manifest,n.TDpred,n.TIpred,startoffset=0){
#   `[` <- function(..., drop=FALSE) base::`[`(...,drop=drop) #to set DROP=FALSE on all bracket subset operations
  message(paste0("Converting intervals to absolute time:  Any missing intervals on 1st row of each subject are assumed to occur at earliest measurement time (", startoffset ,"), any other missing intervals render subsequent intervals for the subject unusable so time variables are set NA"))

  initialmissingcount <- ifelse(is.na(datalong[1,1+n.manifest+n.TDpred+1]),1,0)
  othermissingcount<-0
  datalong[1,1+n.manifest+n.TDpred+1]<-sum(c(datalong[1,1+n.manifest+n.TDpred+1],startoffset),na.rm=T) #datalong row 1 equals first interval and offset
  
  for(i in 2:nrow(datalong)){ #for subsequent rows
          if(datalong[i,1]==datalong[i-1,1]){ #check if the subject is the same as the row above
            othermissingcount <- ifelse(is.na(datalong[i,1]==datalong[i-1,1]),othermissingcount+1, othermissingcount)            
            datalong[i,1+n.manifest+n.TDpred+1]<-sum(datalong[(i-1):i,1+n.manifest+n.TDpred+1],na.rm=FALSE)+startoffset #if same subject, sum the new interval with the prev total time
      } else {
        initialmissingcount <- ifelse(is.na(datalong[i,1+n.manifest+n.TDpred+1]),initialmissingcount+1,initialmissingcount)
        datalong[i,1+n.manifest+n.TDpred+1]<-sum(c(datalong[i,1+n.manifest+n.TDpred+1],startoffset),na.rm=T) #otherwise create new total time with new interval and offset
      }
  }
  colnames(datalong)[1+n.manifest+n.TDpred+1]<-'AbsTime'
  return(datalong)
}