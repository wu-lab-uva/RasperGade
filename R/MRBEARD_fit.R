library(bbmle)
#
fitContinuousRateModelFromPIC = function(x,l,model,num_categories=30,start.value=NULL,method=NULL){
  if(any(is.finite(model$boundary))){
    picfunc = match.fun("dDiscreteRateRBMPIC")
  }else{
    picfunc = match.fun("dDiscreteRatePIC")
  }
  LL = function(...){
    argg = c(as.list(environment()))
    this.model = model
    this.model$fullparams[names(argg)] = lapply(1:length(argg),function(i){
      transformParameters(x = argg[[i]],range=model$param_range[[names(argg)[i]]])
    })
    rate = discretizeDistribution(n = num_categories,model = this.model)
    this.model$rate = rate$rate
    this.model$probs = rate$probs
    pp = sapply(1:length(x),function(i){
      picfunc(x = x[[i]],l=l[[i]],model = this.model)})
    #if(any(pp<=0)) print(which(pp<=0))
    ll = -sum(log(pp))
    return(ll)
  }
  formals(LL) = formals(model$LLfunc)
  if(is.null(start.value)){
    start.value = lapply(names(formals(LL)),function(i){
      transformParameters(x=model$fullparams[[i]],range = model$param_range[[i]],inverse = TRUE)
    })
    names(start.value) = names(formals(LL))
  }
  if(is.null(method)){
    if(length(start.value)>1){
      method="Nelder-Mead"
    }
    else{
      method="BFGS"
    }
  }
  result = do.call("mle2", 
                   c(list(LL, 
                          start = start.value, 
                          data = list(x = x,l=l), 
                          method = method,
                          control = list(trace=TRUE)
                   )))
  model$fullparams[names(result@coef)] = lapply(1:length(result@coef),function(i){transformParameters(x=result@coef[i],
                                                                                  range=model$param_range[[names(result@coef)[i]]])})
  rate = discretizeDistribution(n = num_categories,model = model)
  model$rate = rate$rate
  model$probs = rate$probs
  model$AIC = AIC(result)
  return(model)
}
#
addZero2DiscretizedDistribution = function(p0,model,k=1e-6){
  if(is.null(model$rate)|is.null(model$probs)) stop("The rate distribution must be discretized first.")
  model$rate = c(k*min(model$rate),model$rate)
  model$probs = c(p0,(1-p0)*model$probs)
  return(model)
}




