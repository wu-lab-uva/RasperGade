#' @import bbmle
#' @import doParallel
#' @title  Find only tip-pair contrasts in a tree
#' @description  Find only tip-pair contrasts in a tree for a reconstruction-free analysis
#' @export
#' @rdname findTerminalPIC
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

#' @title  Reconstruct ancestral states under BM and time-independent variation
#' @description  Reconstruct ancestral states under BM and time-independent variation for optimization of model parameters
#' @export
#' @rdname reconstructAncestralStates
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

#' @title  Wrapper function for optimization
#' @description  Fit pulsed evolution with reliable results
#' @export
#' @rdname fitModel2Step
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

#' @export
#' @rdname fitModel2Step
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

#' @title  Fit pulsed evolution model
#'
#' @description  Fit pulsed evolution model in a maximum-likelihood framework
#'
#' @param phy a phylo-class object from the `ape` package
#' @param x trait change in `fitPE` and tip trait values in `fitPE.phy`
#' @param l branch length of time/divergence
#' @export
#' @rdname fitPE
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

#' @export
#' @rdname fitPE
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

