# load required packages
library(ape)
library(bbmle)
library(pracma)
library(doParallel)
library(extraDistr)

### function to select only independent tip-pair contrasts
## phy is a phylo-class from ape package
## x is a named vector that contains the trait value for each tip
## remove.zero determines if zero-contrasts are removed from the results
# output value of this function is a list whose elements are vectors/lists
# of PIC, trait value, branch length,node index and tip labels
findTerminalPIC = function (phy, x, remove.zero = FALSE) {
  # add node labels if not already exist
  if (is.null(phy$node.label)) 
    phy = makeNodeLabel(phy)
  # keep only tip-values
  x = x[phy$tip.label]
  # find nodes that have only tip-children
  tip.node = setdiff(sapply(1:Ntip(phy), getLatestAncestor, 
                            phy = phy), sapply(Ntip(phy) + 1:phy$Nnode, getLatestAncestor, 
                                               phy = phy))
  # tip labels for the tip-pairs
  tip.label = lapply(tip.node, function(n) {
    phy$tip.label[phy$edge[phy$edge[, 1] == n, 2]]
  })
  # branches leading to the tips
  term.l = lapply(tip.node, function(n) {
    phy$edge.length[phy$edge[, 1] == n]
  })
  # tip values of the tip-pairs
  term.x = lapply(tip.node, function(n) {
    x[phy$edge[phy$edge[, 1] == n, 2]]
  })
  # get all phylogenetically independent contrasts
  phy.pic = pic(x = x, phy = phy)
  # get PIC of tip-pairs
  term.pic = phy.pic[tip.node - Ntip(phy)]
  # remove zero-contrasts if required
  if (remove.zero) {
    p0 = sum(term.pic == 0)/length(term.pic)
    idx = term.pic != 0
    return(list(pic = term.pic[idx], x = term.x[idx], l = term.l[idx], 
                node = tip.node[idx],tip = tip.label[idx] ,p0 = p0))
  }
  else {
    return(list(pic = term.pic, x = term.x, l = term.l, node = tip.node, tip = tip.label))
  }
}

### function to quickly reconstruct ancestral states and calculate PIC with time-independent variation
## phy is a phylo-class from ape package
## x is a named vector that contains the trait value for each tip
## rate is the total rate of evolution (variance of trait change per unit branch length)
## epsilon is the variance of time-independent variation between nodes
# output value of this function is a list whose elements are vectors of tip and ancestral states, tip and
# ancestral uncertainty, and trait change, branch length, total uncertainty, total variance between nodes
reconstructAncestralStates = function(phy,x,rate,epsilon){
  # pre-allocate ancestral states and their corresponding variance
  trait = c(x[phy$tip.label],numeric(phy$Nnode))
  trait.var = c(rep(epsilon/2,Ntip(phy)),numeric(phy$Nnode))
  # pre-allocate trait difference, branch length and time-independent variation between nodes
  dx = numeric(phy$Nnode)
  dl = dx
  de = dx
  # find the order of branches to reconstruct from tip to root
  idx = reorder.phylo(x = phy,order = "postorder",index.only = TRUE)
  # reconstruct ancestral state using Felsenstein's method
  for(i in seq(1,length(idx),2)){
    # find the index of ancestor and descendants 
    anc.idx = phy$edge[idx[i],1]
    des.idx = phy$edge[idx[c(i,i+1)],2]
    # find trait values and branch lengths
    this.x = trait[des.idx]
    this.l = phy$edge.length[idx[c(i,i+1)]]
    # find time-independent variation
    this.epsilon = trait.var[des.idx]
    # reconstruct ancestral states by Felsenstein's PIC method
    trait[anc.idx] = (this.x[1]*(this.l[2]*rate+this.epsilon[2])+
                        this.x[2]*(this.l[1]*rate+this.epsilon[1]))/
      (sum(this.l)*rate+sum(this.epsilon))
    # the uncertainty (time-independent variation) of reconstructed ancestral states
    trait.var[anc.idx] = ((this.l[2]*rate+this.epsilon[2])*
                            (this.l[1]*rate+this.epsilon[1]))/
      (sum(this.l)*rate+sum(this.epsilon))
    # get trait change, total branch length and total time-independent variation betwee nodes
    dx[anc.idx-Ntip(phy)] = this.x[1]-this.x[2]
    dl[anc.idx-Ntip(phy)] = sum(this.l)
    de[anc.idx-Ntip(phy)] = sum(this.epsilon)
  }
  return(list(trait=trait,var=trait.var,contrast=dx,l = dl,epsilon=de,lprime=dl*rate+de))
}

### function to calculate pseudo-PIC that normalizes trait changes into standard normal distribution
## phy is a phylo-class from ape package
## x is a named vector that contains the trait value for each tip
## pfunc is one of the cumulative probability function in PEpois.R
## params is the parameter of the pfunc
# output value of this function is a vector of pseudo-PIC
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

### function to get parameters from the fitting output
## x is a list with an element named params, which is the output of fitting functions in this script
get.params = function(x){x$params}

### function to get AIC from the fitting output
## x is a list with an element named AIC, which is the output of fitting functions in this script
get.AIC = function(x){x$AIC}

### function to calculate the process variance
## params is a named vector of parameters whose elements is oredered as:
## [lambda, size], sigma, epsilon
## the elements in the bracket can repeat several times
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

### function to log-transform the model parameters before fitting
## x is a named vector of parameters
## size.ratio is the constraint of jump size between each Poisson process
initialize.parameters = function(x,size.ratio){
  numPE = length(x)/2-1
  for(i in seq(numPE,2,-1)){
    x[paste0("size",i)]=x[paste0("size",i-1)]/x[paste0("size",i)]-size.ratio
  }
  return(as.list(log(x)))
}

### function to log-transform the model parameters before fitting (PE1)
## x is a named vector of parameters
initialize.parameters.SP = function(x,...){as.list(log(x))}

### function to log-transform the model parameters before fitting (PE2)
## x is a named vector of parameters
initialize.parameters.MP = function(x,...){
  x["size2"]=x["size1"]/x["size2"]-size.ratio
  return(as.list(log(x)))
}

### function to log-transform the model parameters before fitting (PE3)
## x is a named vector of parameters
initialize.parameters.TP = function(x,...){
  x["size3"]=x["size2"]/x["size3"]-size.ratio
  x["size2"]=x["size1"]/x["size2"]-size.ratio
  return(as.list(log(x)))
}

### function to fit a model through high-low tolerance phases
## fit.func is one of the model fitting functions defined in this script
## start.params is a list whose first element contains the list of starting model parameters
## params1 and params2 are the other arguments passed to the fitting function in the two phases
## AIC.func is the function to extract AIC from the results
## params.function is the function to extract parameters from the results
## ini.function is the function to transform the model parameters before fitting
## AIC.cutoff is the threshold to stop fitting
# output value of this function is a list that contains the final AIC, final parameter and results of each phase
fit.model.byNM.2step = function(fit.func="mle2",start.params,params1,params2,
                                AIC.func="AIC",params.func=function(x){x@coef},ini.func="as.list",
                                AIC.cutoff=2,...){
  # fit high-tolerance phase model
  res1 = fit.model.byNM(fit.func=fit.func,start.params = start.params,params=params1,
                        AIC.func=AIC.func,params.func=params.func,ini.func=ini.func,
                        AIC.cutoff=AIC.cutoff,...)
  # if infinite AIC is returned, stop fitting
  if(is.infinite(res1$AIC)) return(list(params = res1$params, AIC = res1$AIC,high=res1,low=NULL))
  # update starting parameters
  start.params[[1]] = do.call(ini.func,list(res1$params))
  # fit low-tolerance phase model
  res2 = fit.model.byNM(fit.func=fit.func,start.params=start.params,params=params2,
                        AIC.func=AIC.func,params.func=params.func,ini.func=ini.func,
                        AIC.cutoff=AIC.cutoff,...)
  return(list(params = res2$params, AIC = res2$AIC,high=res1,low=res2))
}

### function to fit a model by maximum likelihood with multiple steps
## fit.func is one of the model fitting functions defined in this script
## start.params is a list whose first element contains the list of starting model parameters
## params is the other arguments passed to the fitting function
## AIC.func is the function to extract AIC from the results
## params.function is the function to extract parameters from the results
## ini.function is the function to transform the model parameters before fitting
## AIC.cutoff is the threshold to stop fitting
# output value of this function is a list that contains the final AIC, final parameter and results of each step
fit.model.byNM = function(fit.func="mle2",start.params,params,
                          AIC.func="AIC",params.func=function(x){x@coef},ini.func="as.list",
                          AIC.cutoff=2,...){
  # initialize AICs
  lastAIC = Inf
  deltaAIC = Inf
  # initialize index
  i = 0
  # get initial parameters
  ini.params = start.params[[1]]
  # evaluate starting AIC
  res = list(do.call(fit.func,c(params,start.params,list(eval.only=TRUE))))
  # if starting AIC is infinite, stop fitting
  if(is.infinite(res[[1]]$AIC)) return(list(AIC = res[[1]]$AIC,params=res[[1]]$params,full=res))
  # continues to fit until AIC cutoff is reached
  while(deltaAIC>AIC.cutoff){
    t0= Sys.time()
    # update starting parameters
    start.params[[1]] = ini.params
    # fit the model again
    mdl = do.call(fit.func,c(params,start.params))
    # update the model parameters
    current.params = do.call(params.func,list(mdl))
    # transform the new model parameters
    ini.params = do.call(ini.func,list(current.params))
    # update AIC
    currentAIC = do.call(AIC.func,list(mdl))
    deltaAIC = lastAIC - currentAIC
    lastAIC = currentAIC
    # update index
    i = i+1
    t1=Sys.time()
    res[[i+1]] = list(AIC=currentAIC,params=current.params,time=list(t0=t0,t1=t1))
  }
  return(list(AIC = currentAIC,params=current.params,full=res))
}


### functions to fit BM and PE1 models without a phylogeny
## x is a vector of trait changes between nodes
## l is a vector of branch lengths between nodes
## start.value is a list of starting parameters
## numCores is the number of cores when run in parallel
## fixed is a list of fixed parameters
## laplace determines if Laplace distribution is used for time-independent variation
## eval.only determines if the function should only evaluate the AIC of starting parameters
## reltol is the relative tolerance of the optimizer
# output of these functions is a list of AIC, model parameter and detailed optimization results.
fitPE = function(x,l,start.value=NULL,numCores=1,fixed=list(),laplace=TRUE,eval.only=FALSE,reltol=sqrt(.Machine$double.eps),...){
  # set the probability function to use
  if(laplace){
    dfunc = dPEpoislaplace
  }else{
    dfunc = dPEpoisnorm
  }
  # set the likelihood function
  if(numCores>1){
    # Setup parallel backend to use many processors
    cores=detectCores()
    cat(sprintf("There's %d cores available\n",cores))
    #cl <- makeCluster(min(cores[1]-1,numCores),type="PSOCK") # Not to overload
    cl <- makeCluster(min(cores[1]-1,numCores),type = "FORK") # Not to overload
    cat(sprintf("Using %d cores\n",min(cores[1]-1,numCores)))
    registerDoParallel(cl)
    cat(sprintf("%d cluster registered\n",min(cores[1]-1,numCores)))
    LL = function(lambda,size,sigma,epsilon){
      lambda= exp(lambda)
      size = exp(size)
      sigma=exp(sigma)
      epsilon=exp(epsilon)
      pp = parSapply(cl,1:length(x),function(i){dfunc(x = x[i],t = l[i],lambda = lambda,size=size,sigma = sigma,epsilon=epsilon)})
      ll = -sum(log(pp))
      return(ll)
    }
  }else{
    LL = function(lambda,size,sigma,epsilon){
      lambda= exp(lambda)
      size = exp(size)
      sigma=exp(sigma)
      epsilon=exp(epsilon)
      pp = sapply(1:length(x),function(i){dfunc(x = x[i],t = l[i],lambda = lambda,size=size,sigma = sigma,epsilon=epsilon)})
      ll = -sum(log(pp))
      return(ll)
    }
  }
  # if there is no starting parameters, set one
  if(is.null(start.value)) start.value = list(lambda=log(1),size=log(1),sigma = log(var(x/sqrt(l))),epsilon=log(var(x)))
  # remove fixed parameters from the starting parameters
  start.value[names(fixed)]=NULL
  # set the method to do optimization
  if(length(start.value)<2){
    method  = "BFGS"
  }else{
    method = "Nelder-Mead"
  }
  # evaluate AIC only if instructed
  if(eval.only){
    res = exp(unlist(c(start.value,fixed)))
    names(res) = c(names(start.value),names(fixed))
    res = res[c("lambda","size","sigma","epsilon")]
    this.AIC = 2*do.call(LL,c(start.value,fixed))+2*length(start.value)
    if(numCores>1) stopCluster(cl)
    return(list(params=res,
                AIC=this.AIC,
                raw=NA))
  }
  # optimize for the best model parameters
  result = do.call("mle2", 
                   c(list(LL, 
                          start = start.value, 
                          data = list(x=x),
                          fixed = fixed,
                          method = method,
                          control = list(trace=TRUE,maxit=2000,reltol=reltol)
                   )))
  # transform the parameters back to the linear scale 
  res = exp(result@fullcoef)
  # stop parallel clusters
  if(numCores>1) stopCluster(cl)
  return(list(params=res,AIC=AIC(result),raw= result))
}

### functions to fit BM and PE1 models along a phylogeny
## phy is a phylo-class object from ape package
## x is a vector of trait values of the tips
## start.value is a list of starting parameters
## numCores is the number of cores when run in parallel
## fixed is a list of fixed parameters
## ignore.zero determines if zero-contrasts are removed from the optimization
## laplace determines if Laplace distribution is used for time-independent variation
## eval.only determines if the function should only evaluate the AIC of starting parameters
## bootstrap is a vector of node indexes that should be used in the optimization 
## reltol is the relative tolerance of the optimizer
# output of these functions is a list of AIC, model parameter and detailed optimization results.
fitPE.phy = function(phy,x,start.value=NULL,numCores=1,fixed=list(),ignore.zero=TRUE,
                     laplace=TRUE,eval.only=FALSE,bootstrap=NULL,reltol=sqrt(.Machine$double.eps),...){
  # initialize node indexes if not supplied
  if(is.null(bootstrap)) bootstrap = 1:phy$Nnode
  # remove zero-contrasts if instructed
  if(ignore.zero) bootstrap = intersect(bootstrap,
                                        which(reconstructAncestralStates(phy = phy,x=x[phy$tip.label],rate = 1,epsilon = 0)$contrast!=0))
  # set probability function to use
  if(laplace){
    dfunc = dPEpoislaplace
  }else{
    dfunc = dPEpoisnorm
  }
  # set likelihood function
  if(numCores>1){
    # Setup parallel backend to use many processors
    cores=detectCores()
    cat(sprintf("There's %d cores available\n",cores))
    #cl <- makeCluster(min(cores[1]-1,numCores),type="PSOCK") # Not to overload
    cl <- makeCluster(min(cores[1]-1,numCores),type = "FORK") # Not to overload
    cat(sprintf("Using %d cores\n",min(cores[1]-1,numCores)))
    registerDoParallel(cl)
    cat(sprintf("%d cluster registered\n",min(cores[1]-1,numCores)))
    LL = function(lambda,size,sigma,epsilon){
      lambda= exp(lambda)
      size = exp(size)
      sigma=exp(sigma)
      epsilon=exp(epsilon)
      rate = lambda*size+sigma
      anc = reconstructAncestralStates(phy = phy,x=x[phy$tip.label],rate = rate,epsilon = epsilon)
      pp = parSapply(cl,bootstrap,function(i){
        dfunc(x = anc$contrast[i],t = anc$l[i],lambda = lambda,size=size,sigma = sigma,epsilon=anc$epsilon[i])
      })
      ll = -sum(log(pp))
      return(ll)
    }
  }else{
    LL = function(lambda,size,sigma,epsilon){
      lambda= exp(lambda)
      size = exp(size)
      sigma=exp(sigma)
      epsilon=exp(epsilon)
      rate = lambda*size+sigma
      anc = reconstructAncestralStates(phy = phy,x=x[phy$tip.label],rate = rate,epsilon = epsilon)
      pp = sapply(bootstrap,function(i){
        dfunc(x = anc$contrast[i],t = anc$l[i],lambda = lambda,size=size,sigma = sigma,epsilon=anc$epsilon[i])
      })
      ll = -sum(log(pp))
      return(ll)
    }
  }
  # if no starting parameters are supplied, set one
  if(is.null(start.value)) start.value = list(lambda=log(1),size=log(1),sigma = log(var(pic(x = x[phy$tip.label],phy = phy))),epsilon=log(1))
  # remove fixed parameters from the starting parameters
  start.value[names(fixed)]=NULL
  # set optimization method
  if(length(start.value)<2){
    method  = "BFGS"
  }else{
    method = "Nelder-Mead"
  }
  # evaluate AIC only if instructed
  if(eval.only){
    res = exp(unlist(c(start.value,fixed)))
    names(res) = c(names(start.value),names(fixed))
    res = res[c("lambda","size","sigma","epsilon")]
    this.AIC = 2*do.call(LL,c(start.value,fixed))+2*length(start.value)
    if(numCores>1) stopCluster(cl)
    return(list(params=res,
                AIC=this.AIC,
                raw=NA))
  }
  # search for the best model parameters
  result = do.call("mle2", 
                   c(list(LL, 
                          start = start.value, 
                          data = list(x=x),
                          fixed = fixed,
                          method = method,
                          control = list(trace=TRUE,maxit=2000,reltol=reltol)
                   )))
  # transform model parameters back to the linear scale
  res = exp(result@fullcoef)
  # stop parallel clusters
  if(numCores>1) stopCluster(cl)
  return(list(params=res,AIC=AIC(result),raw= result))
}

### function to simulate trait change between nodes
## l is a vector of branch lengths
## lambda, size, sigma and epsilon are parameters as defined in PEpois.R
## nj specifies the number of jumps between a specific pair of nodes, if supplied
# the output of this function is a list whose elements are vectors of trait change
# introduced by different processes
simulatePET = function(l,align.length=NA,lambda,size,sigma,epsilon,nj=NULL){
  # get the number of trait changes to simulate
  n = length(l)
  # adjust branch length if necessary (not used)
  if(is.na(align.length)){
    t=l
  }else{
    t = rgamma(n=n,shape=l*align.length+1,scale = 1/align.length)
  }
  # generate random number of jumps
  if(is.null(nj)) nj = rpois(n = n,lambda = lambda*t)
  # generate random trait change by each process
  dxbm = rnorm(n=n)*sqrt(t*sigma)
  dxpe = rnorm(n=n)*sqrt(nj*size)
  ep = runif(n = n)
  ebm = qnorm(p =ep)*sqrt(epsilon)
  elaplace = qlaplace(p = ep)*sqrt(epsilon/4)
  return(list(l = l, t = t,nj=nj,
              dxbm = dxbm, dxpe = dxpe,
              ep=ep,ebm = ebm, elaplace = elaplace))
}

### functions to simulate BM and PE1 trait value over phylogeny
## phy is a phylo-class object from ape package
## lambda, size, sigma and epsilon are parameters as defined in PEpois.R
## start is the trait value of the common ancestor (root)
# the output of this function is a vector of trait values
# BM and PE1 model
simulatePET.phy = function(phy,lambda,size,sigma,epsilon,start=0){
  # initialize trait values
  trait.sim = numeric(length = Ntip(phy)+Nnode(phy))+start
  # set node labels if not already exist
  if(is.null(phy$node.label)){
    names(trait.sim) = c(phy$tip.label,as.character(Ntip(phy)+(1:phy$Nnode)))
  }else{
    names(trait.sim) = c(phy$tip.label,phy$node.label)
  }
  # reorder edges for root-to-tip simulation 
  phy = reorder.phylo(phy,order = "postorder")
  # simulate trait changes between nodes
  trait.delta = simulatePET(l = phy$edge.length,lambda = lambda,size=size,sigma = sigma,epsilon = epsilon)
  # add up trait changes to get trait values from root to tip
  for(i in 1:length(phy$edge.length)){
    l = length(phy$edge.length)+1-i
    trait.sim[phy$edge[l,2]] = trait.sim[phy$edge[l,1]] + trait.delta$dxbm[l] + trait.delta$dxpe[l]
    # add time-independent variation for tips
    if(phy$edge[l,2]<=Ntip(phy)) trait.sim[phy$edge[l,2]] = trait.sim[phy$edge[l,2]] + trait.delta$elaplace[l]
  }
  return(trait.sim)
}

### function to summarize model parameters
## AIC is a vector of AIC values from the optimizer
## scale is the scale used to linearly transform the data
## sample.size is the number of nodes that is evaluated in the likelihood function
# output value of this function is a vector of corrected AIC values
# note that this function should not change the best model for any trait
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

### function to correct for the effect of scaling
## AIC is a vector of AIC values from the optimizer
## scale is the scale used to linearly transform the data
## sample.size is the number of nodes that is evaluated in the likelihood function
# output value of this function is a vector of corrected AIC values
# note that this function should not change the best model for any trait
normalize.AIC = function(AIC,scale=1,sample.size=6667){
  return(AIC-log(scale)*sample.size)
}
