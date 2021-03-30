#' @title  Calculate pseudo-PIC
#' @description  Calculate pseudo-PIC that follows the standard normal distribution
#' @export
#' @rdname pseudoPIC
pseudo.pic = function(phy,x,pfunc,params){
  # get the number of Poisson processes
  numPE = length(params)/2-1
  # calculate the total evolution rate (variance of trait change per unit branch length)
  rate = sum(sapply(1:numPE,function(i){
    prod(params[-1:0 + 2*i])
  }))+params[length(params)-1]
  # reconstruct ancestral states
  epsilon = params[length(params)]
  anc = reconstructAncestralStates(phy = phy,x=x,rate= rate,epsilon=epsilon)
  # calculate the cumulative probability of trait differences between nodes under the distribution defined by pfunc
  pp = sapply(1:length(anc$contrast),function(i){
    do.call(pfunc,c(list(q=anc$contrast[i],t=anc$l[i],epsilon=anc$epsilon[i]),as.list(params[-length(params)])))
  })
  # transform trait differences into normal quantiles (pseudo-PIC)
  p.pic = qnorm(p = pp)
  return(p.pic)
}

#' @title  Get a component of the fitted model
#' @description  Get a component of the fitted model from the list of results
#' @export
#' @rdname getModel
get.params = function(x){x$params}

#' @export
#' @rdname getModel
get.AIC = function(x){x$AIC}

#' @title  Total evolution rate
#' @description  Calculate the long-term evolution rate combining BM and PE
#' @export
#' @rdname processVariance
total.process.variance = function(params){
  # number of Poisson processes
  numPE = length(params)/2-1
  # process variance of each Poisson process
  process.variance = sapply(1:numPE,function(i){
    params[2*i-1]*params[2*i]
  })
  # the total process variance
  return(sum(process.variance)+params[2*numPE+1])
}

#' @title  Prepare model parameters for optimization
#' @description  Prepare the parameter vector so it can be supplied to the optimizer
#' @export
#' @rdname initializeParams
initialize.parameters.SP = function(x,...){as.list(log(x))}

#' @title  Summarize fitted model
#' @description  Summarize fitted model for properly scaled AIC and parameters
#' @export
#' @rdname summarizeModel
summarize.model.parameter = function(params,time=2000,scale=1){
  # get the number of Poisson processes
  numPE = length(params)/2-1
  # correct for scaling
  unscale.coef = rep(scale,length(params))
  unscale.coef[2*(1:numPE)-1] = 1
  unscale.params = params/unscale.coef
  # calculating average waiting time
  wait.time = sapply(1:numPE,function(i){
    unname(time/params[2*i-1])
  })
  # calculate relative jump size
  relative.size = sapply(1:numPE,function(i){
    sqrt(unname(params[2*i]/params[2*numPE+2]))
  })
  # calculate jump rate per unit time
  jump.rate = 1/wait.time
  # calculate process variance (trait change variance per unit branch length)
  process.variance = sapply(1:numPE,function(i){
    params[2*i-1]*params[2*i]
  })
  # calculate relative contribution
  contribution = sapply(1:numPE,function(i){
    process.variance[i]/(sum(process.variance)+params[2*numPE+1])
  })
  return(list(summary=data.frame(time=wait.time,rate=jump.rate,size=relative.size,contribution = contribution),params=unscale.params))
}

#' @export
#' @rdname summarizeModel
normalize.AIC = function(AIC,scale=1,sample.size=6667){
  return(AIC-log(scale)*sample.size)
}
