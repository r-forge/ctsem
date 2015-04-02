#' Fits a multiple group continuous time model.
#' 
#' Fits a single continuous time structural equation models to multiple groups (where each group contains 1 or more subjects),
#' by default, all parameters are free across groups.  Can also be used to easily estimate seperate models for each group.
#' 
#' @param datawide Wide format data, as used in \code{\link{ctFit}}.  See \code{\link{ctLongToWide}} to
#' easily convert long format data.
#' @param groupings Vector of character labels designating group membership for each row of datawide.  
#' These will be prefixed on relevant parameter estimates in the summary.
#' @param ctmodelobj Continuous time model to fit, specified via \code{\link{ctModel}} function.
#' @param fixedmodel Modified version of ctmodelobj, wherein any parameters you wish to keep 
#' fixed over groups should be given the value 'groupfixed'.  
#' If specified, all other parameters will be free across groups.
#' @param freemodel Modified version of ctmodelobj, wherein any parameters you wish to free across groups
#' should be given the label 'groupfree'.  
#' If specified, all other parameters will be fixed across groups.  
#' If left NULL, the default, all parameters are free across groups.
#' @param confidenceintervals Character vector of parameter labels to estimate confidence intervals for.  
#' Unlike with \code{\link{ctFit}}, entire matrices cannot be specified.
#' @param showInits Displays start values prior to optimization
#' @param fitasgroup If FALSE, OpenMx models for each group are fit individually and output in a list.  
#' By default, groups are specified and fit as sub-models of a single large OpenMx model.
#' @param plots Only relevant if fitasgroup=FALSE.  If TRUE, plots each model fit 
#' with mean trajectory, autoregression and cross regression plots.
#' @param summaries Same as plots but for summary output.
#' @param carefulFit if TRUE, first fits the specified model with a penalised likelihood function 
#' to force MANIFESTVAR, DRIFT, TRAITVAR, MANIFESTTRAITVAR parameters to remain close to 0, then
#' fits the specified model normally, using these estimates as starting values. 
#' Can help with optimization, though results in user specified inits being ignored for the final fit.
#' @param freemodelvarlimits Experimental
#' @param varzeroweight Experimental
#' @param retryattempts Number of fit retries to make.
#' @param plotOptimization plots graphs of optimization progress after fitting, uses OpenMx checkpointing.
#' @param ... additional arguments to pass to \code{\link{ctFit}}.
#' @return Returns either an OpenMx fit object (if fitasgroup=TRUE), or a list containing 
#' individual ctsem fit objects.
#' @details Additional \code{\link{ctFit}} parameters may be specified as required.
#' @examples #Two group model, all parameters except LAMBDA[3,1] constrained across groups.
#' data(ctExample4)
#' basemodel<-ctModel(n.latent=1, n.manifest=3, Tpoints=20,
#'                    LAMBDA=matrix(c(1, 'lambda2', 'lambda3'), nrow=3, ncol=1),
#'                    MANIFESTMEANS=matrix(c(0, 'manifestmean2', 'manifestmean3'), 
#'                    nrow=3, ncol=1), TRAITVAR = 'auto')
#' 
#' freemodel<-basemodel
#' freemodel$LAMBDA[3,1]<-'groupfree'
#' groups<-paste0('g',rep(1:2, each=10),'_')
#' 
#' multif<-ctMultigroupFit(datawide=ctExample4, groupings=groups,
#'                        ctmodelobj=basemodel, freemodel=freemodel)
#' summary(multif)
#' 
#' 
#' \dontrun{
#' #fixed model approach
#' fixedmodel<-basemodel
#' fixedmodel$LAMBDA[2,1]<-'groupfixed'
#' groups<-paste0('g',rep(1:2, each=10),'_')
#' 
#' multif<-ctMultigroupFit(datawide=ctExample4, groupings=groups,
#'                        ctmodelobj=basemodel, fixedmodel=fixedmodel)
#' summary(multif) 
#'}
#' 
#' 
#' @seealso \code{\link{ctFit}} and \code{\link{ctModel}}
#' @export
#' @import OpenMx


ctMultigroupFit<-function(datawide,groupings,ctmodelobj,fixedmodel=NA,freemodel=NA,
  freemodelvarlimits=NA,varzeroweight=0,fitasgroup=TRUE, carefulFit=FALSE,
  retryattempts=5,showInits=TRUE,confidenceintervals=NULL,
  plots=FALSE,summaries=FALSE,plotOptimization=FALSE,...){

  
  if(any(suppressWarnings(!is.na(as.numeric(groupings))))) stop("grouping variable must not contain purely numeric items")
  if(length(groupings)!= nrow(datawide)) stop('length of groupings does not equal number of rows of datawide')
  
  if(all(is.na(fixedmodel))) fixedmodel<-ctmodelobj #so that it is not null or na
  
  startparams<-c() #to fill as needed
  
  
  
  omxmodels<-list() #blank list preparing for model input
  for(i in unique(groupings)){ #for every specified group
    #     
    singlegroup<-datawide[which(groupings == i),,drop=F] #data for the group
    
    singlectspec<-ctmodelobj
    varparams<-c() #blank placeholder for adding any variance constrained free parameters
    
    if(all(is.na(freemodel))) freemodel <- lapply(ctmodelobj,function(x) { x<-rep('groupfree',length(x))}) #if no freemodel specified then free all params at this point
    
    for(m in 1:length(ctmodelobj)) { #for every element of the ctmodelobj list
      if(is.matrix(ctmodelobj[[m]])){ #if the element is a matrix
      for(j in 1:length(ctmodelobj[[m]])){ #for every slot in the matrix
        if(freemodel[[m]][j]=="groupfree"){ #if the slot is free in freemodel and not fixed in fixedmodel
        jnum<-suppressWarnings(as.numeric(ctmodelobj[[m]][j])) #check if it is numeric
        if(!is.na(ctmodelobj[[m]][j]) && is.na(jnum)) { #if the slot is neither NA or fixed to a value, then
          singlectspec[[m]][j] <- paste0(i,ctmodelobj[[m]][j]) #give the label a group specific prefix
          
#           if(!is.na(suppressWarnings(as.numeric(freemodelvarlimits[[m]][j])))) { #if there are variance limits imposed on this free param
#             newvarparam<-
#             
#             varparams<-cbind(varparams,     
#           }
          
        }
      }
        if(any(!is.na(fixedmodel))){
        if(fixedmodel[[m]][j]=="groupfixed"){
          jnum<-suppressWarnings(as.numeric(ctmodelobj[[m]][j])) #check if it is numeric
          if(!is.null(ctmodelobj[[m]][j]) && is.na(jnum)) { #if the slot is neither null or fixed to a value, then
            singlectspec[[m]][j] <- paste0(ctmodelobj[[m]][j]) #use the global label
          }
        }
      }
    }
      }
    }
   

    
    
    
    if(carefulFit==TRUE) message('Begin carefulFit start value estimation for group ', i)
    omxmodel<-ctFit(singlegroup,singlectspec,nofit=TRUE, carefulFit=carefulFit,...) #omxmodel for group i
    
    if(carefulFit==TRUE) {
      startparams<-c( startparams[ !( names(startparams) %in%  #get inits
          names(omxGetParameters(omxmodel))) ], #that are not found in the new fits inits
        omxGetParameters(omxmodel) ) #and combine the two vectors
    }
    
    omxmodel<-mxRename(omxmodel, newname=i) #change name of omxmodel for group i
      
    
    if(fitasgroup==TRUE) omxmodels[[i]]<-omxmodel #if fitting single multigroup model, add omxmodel for group i to list of omxmodels for all groups
    
    if(fitasgroup==FALSE) { #if fitting groups seperately
    omxfit<-ctFit(singlegroup,singlectspec,...) #generate omxfit for group i
      omxmodels[[i]]<-omxfit #add to list of outputs
      if(plots==TRUE) plot(omxfit,wait=F)
      if(summaries==TRUE) print(summary(omxfit))
      
    }
    
  } #end loop over groups
  
  
  
  if(fitasgroup==TRUE) {
    

    if(any(!is.na(freemodelvarlimits))){ #extracts global parameter names for params with variance limits
   
      varparglobalnames <- suppressWarnings(unlist(ctmodelobj)[unlist(freemodel)=='groupfree' & 
         unlist(fixedmodel)!='groupfixed' &  
          !is.na(suppressWarnings(as.numeric(unlist(freemodelvarlimits))))])
      varparglobalnames<-varparglobalnames[duplicated(varparglobalnames)==FALSE]

       varpargroupnames<-suppressWarnings(matrix(c(  #extracts group specific parameter names for params with variance limits
         paste0(rep(groupings, times=length(varparglobalnames)),
           rep(varparglobalnames,each=length(unique(groupings))))
       ),ncol=length(varparglobalnames)))
      

      L3varlabels<- matrix(paste0(  #setup labels for var/cov matrix of parameters with constrained variance
        'L3var_',
        varparglobalnames,
        '_X_',
        rep(varparglobalnames,each=length(varparglobalnames))), 
        nrow=length(varparglobalnames),
        ncol=length(varparglobalnames))
      
      L3varlabels[upper.tri(L3varlabels)]<-t(L3varlabels)[upper.tri(L3varlabels)]
      
#       test<-test1
#       test[upper.tri(test)]<-t(test)[upper.tri(test)]
      
      
      
      level3model<-OpenMx::mxModel('L3VarianceModel',
        
        mxMatrix(name='L3VarianceParameters',
          nrow=nrow(varpargroupnames), ncol=length(varparglobalnames),
          labels=varpargroupnames, free=TRUE),
        
#         mxMatrix( type="Symm", nrow=length(varparglobalnames), 
#           ncol=length(varparglobalnames),
#           labels=L3varlabels,
#           values=diag(.9,length(varparglobalnames)),
#           free=TRUE, name="expCov" ),
        
#         mxMatrix(name="expMean", nrow=1, ncol=length(varparglobalnames), free=TRUE, 
#           labels=paste0('L3mean_',varparglobalnames)),
        
#         mxMatrix(type='Full',values=length(varparglobalnames), #should this be number of variances or number of both var / cov?
#           ncol=1,nrow=1,free=F,name='L3variablecount'),
#         
        mxMatrix(type='Full',values=length(unique(groupings)), #should this be number of variances or number of both var / cov?
          ncol=1,nrow=1,free=F,name='groupcount'),
      
#         mxMatrix("Full", 1, 1, values = log(2*pi), name = "log2pi"),
        #                           mxAlgebra(
        #                             expression=log2pi %*% 2 + log(det(expCov) + driftvar - expMean) %&% solve(expCov),
        #                             name ="llvector",
        #                           ),
        
        
        
#         mxAlgebra(expression=  1/2 * ( log(det(expCov)) + tr(S %*% solve(expCov)) - log(det(S))  ) + # -(varcount)
#             1/2 * ( t(m - expMean) %*% solve(expCov) %*% (m-expMean)  ),  name ="llvector"),
          
#         mxAlgebra(  name ="llvector", expression= sum( (expCov - S ) * (expCov - S ) ) +  #working but slow
#             sum( (m-expMean) * (m-expMean) )),
        
#         
#         mxAlgebra(  name ="llvector", expression= (groupcount-1) %*% 
#             (log(det(expCov)) - log(det(S)) + tr(S %*% solve(expCov)) - L3variablecount +
#                 (groupcount / (groupcount - 1)) %*% (m - expMean) %*% solve(expCov) %*% t(m - expMean) +1) ),
        
#         mxMatrix(name='llvector',nrow=1,ncol=2,values=.3,free=F),
        
        
          #                             expression=  -2 * (  (-20/2) * log(2*pi) - (20/2) *(log(expCov^2)) - (1/(2*expCov^2)) * sum(driftvar - expMean)^2   ),
          #                             expression = (-varcount / 2) * log(2*pi) - varcount/2*log(det(expCov)) - 
          #                               (1/2)*sum(  t(driftvar - expMean) %*% solve(expCov) %*% (driftvar - expMean)
          
        
        
        
        mxMatrix(name='sumMatrix',values=1,free=F,nrow=1,ncol=length(unique(groupings))),
        
#        
#         l3VarParDeviance <- 
#           mxEval(L3VarianceParameters - bigMeans, fullmodel$submodels$L3VarianceModel,compute=T)
        
        mxAlgebra(name='m', sumMatrix %*% (L3VarianceParameters) %x% (1/groupcount),dimnames=list(NULL,varparglobalnames)), #observed means vector
        mxMatrix(name='bigMeans',
          labels=paste0(
          'm[1,',
          rep(1:length(varparglobalnames), each= length(unique(groupings))),
          ']'),
          ncol=length(varparglobalnames),
          nrow=length(unique(groupings))),
        
        mxAlgebra(name='L3VarParDeviance', L3VarianceParameters - bigMeans),
        mxAlgebra(name='S',t(L3VarParDeviance) %*% L3VarParDeviance,dimnames=list(varparglobalnames,varparglobalnames)), #observed covariance matrix
        mxMatrix(name='varzeroweight',values=varzeroweight, nrow=1, ncol=1),
#         mxAlgebra((sum(llvector)+sum(expCov * expCov)) %x% varzeroweight, name='fitAlgebra'),
        mxAlgebra(sum(tr(S) * tr(S)) * varzeroweight, name='fitAlgebra'),
        mxFitFunctionAlgebra('fitAlgebra')
        
      )      #end L3 model spec
      
      fullmodel <- OpenMx::mxModel('ctsemMultigroup', #output multigroup omxmodel
        mxFitFunctionMultigroup(c(paste0(unique(groupings)),'L3VarianceModel.fitfunction')),
        omxmodels, level3model)
    }  #end variance constraint section
    
    
    if(all(is.na(freemodelvarlimits))) fullmodel <- OpenMx::mxModel('ctsem multigroup', #output multigroup omxmodel
      mxFitFunctionMultigroup(c(paste0(unique(groupings)))),
      omxmodels)
    
    if(carefulFit==TRUE) fullmodel<-omxSetParameters(fullmodel,labels=names(startparams),values=startparams)
    
    fullmodel<-omxAssignFirstParameters(fullmodel)

    if(!is.null(confidenceintervals)) fullmodel <- OpenMx::mxModel(fullmodel, mxCI(confidenceintervals,interval = 0.95,type = "both")) #if 95% confidence intervals are to be calculated

    if(plotOptimization==T){

      fullmodel<-mxOption(fullmodel,'Always Checkpoint', 'Yes')
      fullmodel<-mxOption(fullmodel,'Checkpoint Units', 'iterations')
      fullmodel<-mxOption(fullmodel,'Checkpoint Count', 1)    
    }
    
    multiout<-OpenMx::mxTryHard(fullmodel,
      showInits=showInits,
      #         intervals = ifelse(!is.null(confidenceintervals),TRUE,FALSE),
#       confidenceintervals=confidenceintervals,
      bestInitsOutput=FALSE,
      extraTries=retryattempts,loc=1,scale=.2,paste=FALSE,iterationSummary=TRUE) 
    
    
    if(plotOptimization==TRUE){
      
      checkpoints<-read.table(file='ctsemMultigroup.omx',header=T,sep='\t')
      mfrow<-par()$mfrow
      par(mfrow=c(3,3))
      for(i in 6:ncol(checkpoints)) {
        plot(checkpoints[,i],main=colnames(checkpoints)[i])
      }
      par(mfrow=mfrow)
      deleteCheckpoints <- readline('Remove created checkpoint file, ctsemMultigroup.omx? y/n \n')
      if(deleteCheckpoints) file.remove(file='ctsemMultigroup.omx')
    }
    
    
    return(multiout)
  }
  
  if(fitasgroup==FALSE) return(omxmodels)
  
}


