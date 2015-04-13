# ctmxTryHard 
# Wrapper to mxRun that makes multiple attempts to reach an acceptable solution.
# Temporarily included updated version of mxTryHard from OpenMx within ctsem and renamed to avoid confusion. 

ctmxTryHard<-function (model, extraTries = 10, greenOK = FALSE, loc = 1, 
  scale = 0.25, checkHess = TRUE, fit2beat = Inf, paste = TRUE,
  iterationSummary=FALSE, bestInitsOutput=TRUE, showInits=FALSE,
  ...) 
{
  
  stopflag <- FALSE
  numdone <- 0
  bestfitsofar<-Inf
  inits<-omxGetParameters(model)
  while (!stopflag) {
    if(iterationSummary==TRUE) message(paste0('\nBegin fit attempt ', numdone+1, ' of at maximum ', extraTries +1, ' tries'))
    params <- omxGetParameters(model)
    if(showInits==TRUE) {
      message('Starting values:  ')
      message(paste0(names(params),' : ', params,'\n'))
    }
    numdone <- numdone + 1
    
    fit <- suppressWarnings(try(mxRun(model, suppressWarnings = T, unsafe=T, silent=T,intervals=FALSE)))
    if (class(fit) == "try-error" || fit$output$status$status== -1) {
      newparams<-omxGetParameters(model) #get recent fit
      if(exists('bestfit')) newparams<-bestfit.params #if bestfit exists use this instead
      if(numdone %% 4 == 0) newparams<-inits #sometimes, use initial start values instead
      #       if(numdone %% 5 == 0) { #sometimes, switch optimizers
      #         if(mxOption(NULL, "Default optimizer")=='CSOLNP') newoptimizer<-'NPSOL'
      #         if(mxOption(NULL, "Default optimizer")=='NPSOL') newoptimizer<-'CSOLNP'
      #         message(paste0('Switching to ',newoptimizer,' optimizer for OpenMx temporarily')) 
      #         mxOption(NULL, "Default optimizer", newoptimizer)
      #       }
      model <- omxSetParameters(model, labels = names(newparams), 
        values = newparams * runif(length(params),loc-scale,loc+scale))  #set to multiply bestfit.params instead of params
    }
    else { #if fit was not an error
      if (fit$output$minimum <= bestfitsofar) {
        bestfit <- fit
        bestfit.params <- omxGetParameters(bestfit)
      }
      
      if(fit$output$minimum < bestfitsofar) bestfitsofar <- fit$output$minimum
      
      if (length(fit$output$calculatedHessian) == 0) {
        checkHess <- FALSE
      }
      if (checkHess) {
        if (sum(is.na(fit$output$calculatedHessian)) > 
            0) {
          checkHess <- FALSE
        }
      }
      
      stopflag <- ifelse(checkHess, (fit$output$status[[1]] <= 
          greenOK) & (all(eigen(fit$output$calculatedHessian, 
            symmetric = T, only.values = T)$values > 0)) & 
          (fit$output$minimum <= fit2beat) & (fit$output$minimum <= bestfitsofar), (fit$output$status[[1]] <=  #added bestfitsofar condition
              greenOK) & (fit$output$minimum <= fit2beat) & (fit$output$minimum <= bestfitsofar) )
      if (!stopflag) {
        model <- fit
        newparams<-omxGetParameters(fit) #get recent fit
        if(exists('bestfit')) newparams<-bestfit.params #if bestfit exists use this instead
        if(numdone %% 4 == 0) newparams<-inits #sometimes, use initial start values instead
        #         if(numdone %% 5 == 0) { #sometimes, switch optimizers
        #           if(mxOption(NULL, "Default optimizer")=='CSOLNP') newoptimizer<-'NPSOL'
        #           if(mxOption(NULL, "Default optimizer")=='NPSOL') newoptimizer<-'CSOLNP'
        #           message(paste0('Switching to ',newoptimizer,' optimizer for OpenMx temporarily')) 
        #           mxOption(NULL, "Default optimizer", newoptimizer)
        #         }
        model <- omxSetParameters(model, labels = names(params), 
          values = newparams * runif(length(params),loc-scale,loc+scale))
        fit2beat <- ifelse(fit$output$minimum < fit2beat, fit$output$minimum, 
          fit2beat)
        if(iterationSummary==TRUE){
          message(paste0("Attempt ",numdone," fit:  "))
          message(paste(names(params),": ", fit$output$estimate,"\n"))
          message(paste0("-2LL = ", fit$output$Minus2LogLikelihood))
        }
      }
      
      if(stopflag){
        message('\nSolution found\n')
        fit<-bestfit
        if(length(fit$intervals)>0){ #only calculate confidence intervals once the best fit is established
          fit<-omxSetParameters(fit, labels=names(bestfit.params),values=bestfit.params)
          message("Refit using best inits and estimate confidence intervals\n") 
          #           mxOption(NULL, "Default optimizer", "NPSOL")
          cifit<-suppressWarnings(try(mxRun(fit,intervals=TRUE,suppressWarnings=T,silent=T)))
          if(class(cifit) == "try-error" || cifit$output$status$status== -1) {
            message('Confidence interval estimation generated errors\n')
          } else {
            if (length(OpenMx::summary(cifit)$npsolMessage) > 0) message('Warning messages generated from confidence interval refit\n')
            fit<-cifit
          }
          
        }
        if (length(OpenMx::summary(fit)$npsolMessage) > 0) {
          warning(OpenMx::summary(fit)$npsolMessage)
        }
        
        params <- bestfit.params
        
        if(iterationSummary==TRUE){
          message(paste(names(bestfit.params),": ", bestfit$output$estimate,"\n"))
          message(paste0("-2LL = ", bestfit$output$Minus2LogLikelihood))
        }
        
      }
    } #end 'if fit not an error' section
    if (numdone > extraTries & stopflag==FALSE) { #added stopflag==FALSE
      message('\nRetry limit reached')
      stopflag <- TRUE
      if (exists("bestfit")) {
        fit <- bestfit
        params <- bestfit.params
        if(length(fit$intervals)>0){ #calculate intervals for best fit, even though imperfect
          fit<-omxSetParameters(fit, labels=names(bestfit.params),values=bestfit.params)
          message("Refit using best inits and estimate confidence intervals\n") 
          #           mxOption(NULL, "Default optimizer", "NPSOL")
          cifit<-suppressWarnings(try(mxRun(fit,intervals=TRUE,suppressWarnings=T,silent=T)))
          if(class(cifit) == "try-error" || cifit$output$status$status== -1) {
            message('Confidence interval estimation generated errors, returning fit without confidence intervals\n')
          } else {
            fit<-cifit
          }
        }
        if (length(fit$output$status$statusMsg) > 0) { 
          warning(fit$output$status$statusMsg)
        }
        if(fit$output$status$code==6) message('\nUncertain solution found - consider parameter validity, try again, increase extraTries, change inits, change model, or check data!\n')
        if(iterationSummary==TRUE){
          message(paste(names(bestfit.params),": ", bestfit$output$estimate,"\n"))
          message(paste0("-2LL = ", bestfit$output$Minus2LogLikelihood))
        }
      }
      if (!exists("bestfit")) {
        if (length(fit$output$status$statusMsg) > 0) { 
          warning(fit$output$status$statusMsg)
        }
      }
    }
  }
  if(bestInitsOutput){
    message("\nStart values from best fit:")
    if(paste) message(paste(params, sep=",", collapse = ",")) 
    if(!paste)  message(paste(names(params),": ", params,"\n"))
  }
  return(fit)
}
