#' @export
productNormalDensity = function(mean,var){
  if(length(mean)!=length(var)) stop("Mean and variance for every normal density function should be provided.")
  sig2 = 1/sum(1/var)
  mu = sum(mean/var)*sig2
  S = sqrt(2*pi*sig2)/prod(sqrt(2*pi*var))*exp(-0.5*(sum(mean^2/var)-mu^2/sig2))
  return(c(mean=mu,var=sig2,scale = S))
}
#' @export
compoundNormalMoments = function(probs,means,vars){
  probs = probs/sum(probs)
  compound.mean = sum(probs*means)
  compound.var = sum(probs*vars)+sum(probs*(means^2))-compound.mean^2
  compound.skewness = (sum(probs*(means-compound.mean)^3)+sum(3*probs*(means-compound.mean)*vars))/
    sqrt(compound.var)^3
  compound.kurtosis=(sum(probs*(means-compound.mean)^4)+
                       sum(6*probs*(means-compound.mean)^2*vars)+
                       sum(probs*3*vars^2))/
    compound.var^2-3
  return(c(mean=compound.mean,var=compound.var,skewness = compound.skewness,kurtosis = compound.kurtosis))
}
#' @export
findBestNormalApproximation1 = function(probs,means,vars){
  true.moments = compoundNormalMoments(probs = probs,means = means,vars = vars)
  return(data.frame(mean=unname(true.moments[1]),var=unname(true.moments[2]),probs=1))
}
#' @export
findBestNormalApproximation2 = function(probs,means,vars,numCores=1){
  if(length(probs)<2) return(findBestNormalApproximation1(probs=probs,means=means,vars=vars))
  true.moments = compoundNormalMoments(probs = probs,means = means,vars = vars)
  sample.range = sort(unique(means),decreasing = FALSE)
  sample.range = sort(c(sample.range,(sample.range[-1]+sample.range[-length(sample.range)])/2),decreasing = FALSE)
  sample.likelihood = c(0,sapply(sample.range,function(x){
    sum(dnorm(x = x,mean = means,sd = sqrt(vars))*probs)
  }),0)
  sample.peak = which(sapply(2:(length(sample.range)+1),function(i){
    (sample.likelihood[i]>=sample.likelihood[i+1])&(sample.likelihood[i]>=sample.likelihood[i-1])
  }))
  if(length(sample.peak)<2){
    means.order = sort(abs(means-sample.range[sample.peak]),index.return=TRUE)$ix
    if(length(probs)>1e3){
      sample.break = round(seq(from=1,to = length(means.order)-1,length.out = 1000))
    }else{
      sample.break = 1:(length(means.order)-1)
    }
    moments.diff = mclapply(sample.break,function(i){
      this.idx = means.order[1:i]
      moments.1 = compoundNormalMoments(probs = probs[this.idx]/sum(probs[this.idx]),means = means[this.idx],vars = vars[this.idx])
      moments.2 = compoundNormalMoments(probs = probs[-(this.idx)]/sum(probs[-(this.idx)]),means = means[-(this.idx)],vars = vars[-(this.idx)])
      new.df = data.frame(mean=unname(c(moments.1[1],moments.2[1])),var=unname(c(moments.1[2],moments.2[2])),
                          probs=c(sum(probs[this.idx]),sum(probs[-(this.idx)])))
      new.likelihood = c(0,sapply(sample.range,function(x){
        sum(dnorm(x = x,mean = new.df$mean,sd = sqrt(new.df$var))*new.df$probs)
      }),0)
      rss = sum((sample.likelihood-new.likelihood)^2)
      return(list(diff=rss,df=new.df))
    },mc.cores = numCores)
  }else{
    peak.combn = combn(x = sample.peak,m = 2,simplify = FALSE)
    moments.diff = mclapply(peak.combn,function(i){
      new.df = do.call(rbind,lapply(i,function(j){
        this.group = which(((means-mean(sample.range[i]))*(sample.range[j]-mean(sample.range[i])))>0)
        half.group = which(((means-mean(sample.range[i]))*(sample.range[j]-mean(sample.range[i])))==0)
        this.moment = compoundNormalMoments(probs = c(probs[this.group],probs[half.group]/2),
                              means = means[c(this.group,half.group)],vars = vars[c(this.group,half.group)])
        data.frame(mean=unname(this.moment[1]),var=unname(this.moment[2]),
                   probs=sum(c(probs[this.group],probs[half.group]/2)))
      }))
      new.likelihood = c(0,sapply(sample.range,function(x){
        sum(dnorm(x = x,mean = new.df$mean,sd = sqrt(new.df$var))*new.df$probs)
      }),0)
      rss = sum((sample.likelihood-new.likelihood)^2)
      return(list(diff=rss,df=new.df))
    },mc.cores = numCores)
  }
  best.diff = which.min(sapply(moments.diff,function(x){x$diff}))
  return(moments.diff[[best.diff]]$df)
}

#' @export
marginalReconstruction = function(x,l){
  xs = productNormalDensity(x,l)
  return(c(x=unname(xs["mean"]),var=unname(xs["var"]),scale=unname(xs["scale"])))
}
#' @export
marginalReconstructionWithPE = function(x,l,lambdas,sizes,sigma,epsilons,margin=1e-6,numCores=1,asymptotic=5){
  # create jump profiles
  nJ = lapply(l,function(ll){
    this.nJ = lapply(lambdas,function(lambda){
      if((ll*lambda)>asymptotic) return(list(N=ll*lambda,probs=1))
      rg = qpois(p = c(margin/2,1-margin/2),lambda = ll*lambda)
      nn = rg[1]:rg[2]
      dd = dpois(x = nn,lambda = ll*lambda)
      return(list(N=nn,probs=dd))
    })
    full.index = do.call(expand.grid,lapply(this.nJ,function(y){1:length(y$N)}))
    this.nJ.df = do.call(rbind,lapply(1:dim(full.index)[1],function(i){
      data.frame(n=t(sapply(1:dim(full.index)[2],function(j){
        this.nJ[[j]]$N[full.index[i,j]]
      })),probs=prod(sapply(1:dim(full.index)[2],function(j){
        this.nJ[[j]]$probs[full.index[i,j]]
      })))
    }))
    return(this.nJ.df)
  })
  # get all possible combination of jumps
  idxs = do.call(expand.grid,lapply(nJ,function(nn){
    1:dim(nn)[1]
  }))
  # for each combination, reconstruct ancestral state and get error profile
  # when set parameters, only run in parallel if there are more than 2 jump processes
  profile = do.call(rbind,mclapply(1:dim(idxs)[1],function(k){
    i = unlist(idxs[k,])
    pp = prod(sapply(1:length(i),function(j){
      nJ[[j]]$probs[i[j]]
    }))
    ll = sapply(1:length(i),function(j){
      this.nJ = unlist(nJ[[j]][i[j],])
      sum(this.nJ[-length(this.nJ)]*sizes)+l[j]*sigma+epsilons[j]
    })
    res = marginalReconstruction(x,ll)
    return(data.frame(x=unname(res[1]),var=unname(res[2]),scale=unname(res[3]),prior=pp))
  },mc.cores = numCores))
  # mixing ancestral state distribution together
  sumP = sum(profile$prior)
  if(sumP<0.99) warning("Total probability does not converge to one. Try decreasing margin.",immediate. = TRUE)
  profile$probs = profile$scale*profile$prior
  profile$probs = profile$probs/sum(profile$probs)
  meanX = sum(profile$x*profile$probs)
  varX = sum(profile$x^2*profile$probs) - meanX^2 + sum(profile$var*profile$probs)
  #
  if(any(is.nan(c(meanX,varX)))) stop("NaNs produced!")
  if(varX<=0){
    warning("Negative variance produced, enforcing to zero.",immediate. = TRUE)
    varX = 0
  }
  return(list(x = unname(meanX), var = unname(varX),profile=profile))
}
#' @export
crossValidationWithPE = function(FMR,add.epsilon=TRUE,laplace=FALSE,numApprox=1,margin=1e-6,numCores=1,asymptotic=5){
  cat("Conducting leave-one-out cross validation...\n")
  # the internal function to predict hidden states for each tip
  cv.func = function(tip,FMR){
    this.match = which(FMR$marginal$orientation==tip)
    error.index = do.call(expand.grid,lapply(1:length(this.match),function(i){
      1:dim(FMR$error[[this.match[i]]]$approximate)[1]
    }))
    rec = do.call(rbind,lapply(1:dim(error.index)[1],function(i){
      this.rec = marginalReconstructionWithPE(x = sapply(1:length(this.match),function(j){
        FMR$error[[this.match[j]]]$approximate$x[error.index[i,j]]
      }),
      l = FMR$phy$edge.length[FMR$marginal$edge[this.match]],
      epsilons = sapply(1:length(this.match),function(j){
        FMR$error[[this.match[j]]]$approximate$var[error.index[i,j]]
      }),
      sigma = FMR$params$sigma,lambdas=FMR$params$lambdas,sizes = FMR$params$sizes,margin = margin,asymptotic = asymptotic)
      this.probs = prod(sapply(1:length(this.match),function(j){
        FMR$error[[this.match[j]]]$approximate$probs[error.index[i,j]]
      }))
      this.error = this.rec$profile
      this.error$probs = this.error$probs*this.probs
      return(this.error)
    }))
    rec$probs = rec$probs/sum(rec$probs)
    rec.moments = compoundNormalMoments(probs = rec$probs,means = rec$x,vars = rec$var)
    this.hsp = data.frame(node=tip,label = FMR$phy$tip.label[tip],
                          x = unname(rec.moments[1]),
                          var = unname(rec.moments[2]),
                          stringsAsFactors = FALSE)
    return(list(error=rec,hsp=this.hsp))
  }
  # run in parallel if required
  cv.error = mclapply(1:Ntip(FMR$phy),cv.func,FMR=FMR,mc.cores = numCores)
  # parse results
  cv = do.call(rbind,lapply(cv.error,function(x){x$hsp}))
  cv.error = lapply(cv.error,function(x){x$error})
  # add time-independent variation to error profiles in needed
  if(add.epsilon){
    cv$var = cv$var +FMR$params$epsilon/2
    if(numApprox>1&laplace){
      cv.error = lapply(cv.error,function(x){
        return(data.frame(x=c(x$x,x$x),
                          var=c(x$var+FMR$params$epsilon/2*exp(-1),x$var+FMR$params$epsilon/2*exp(1)),
                          probs=c(x$probs*exp(1)/(1+exp(1)),x$probs/(1+exp(1)))))
      })
    }else{
      for(tip in 1:Ntip(FMR$phy)){
        cv.error[[tip]]$var = cv.error[[tip]]$var+FMR$params$epsilon/2
      }
    }
  }
  #
  return(list(cv=cv,error=cv.error))
}
#' @export
globalReconstructionWithPE = function(FMR,add.epsilon=TRUE,laplace=FALSE,numApprox=1,margin=1e-6,numCores=1,asymptotic=5){
  cat("Reconstructing global ancestral states...\n")
  global.func = function(n,FMR){
    this.match = rev(which(FMR$marginal$orientation==n))
    error.index = do.call(expand.grid,lapply(1:length(this.match),function(i){
      1:dim(FMR$error[[this.match[i]]]$approximate)[1]
    }))
    rec = do.call(rbind,lapply(1:dim(error.index)[1],function(i){
      this.rec = marginalReconstructionWithPE(x = sapply(1:length(this.match),function(j){
        FMR$error[[this.match[j]]]$approximate$x[error.index[i,j]]
      }),
      l = FMR$phy$edge.length[FMR$marginal$edge[this.match]],
      epsilons = sapply(1:length(this.match),function(j){
        FMR$error[[this.match[j]]]$approximate$var[error.index[i,j]]
      }),
      sigma = FMR$params$sigma,lambdas=FMR$params$lambdas,sizes = FMR$params$sizes,margin = margin,asymptotic = asymptotic)
      this.probs = prod(sapply(1:length(this.match),function(j){
        FMR$error[[this.match[j]]]$approximate$probs[error.index[i,j]]
      }))
      this.error = this.rec$profile
      this.error$probs = this.error$probs*this.probs
      return(this.error)
    }))
    rec$probs = rec$probs/sum(rec$probs)
    rec.moments = compoundNormalMoments(probs = rec$probs,means = rec$x,vars = rec$var)
    this.ace = data.frame(node=n,label = FMR$phy$node.label[n-Ntip(FMR$phy)],
                          x = unname(rec.moments[1]),
                          var = unname(rec.moments[2]),
                          stringsAsFactors = FALSE)
    return(list(error=rec,ace=this.ace))
  }
  global.ace = mclapply((1:Nnode(FMR$phy))+Ntip(FMR$phy),global.func,FMR=FMR,mc.cores = numCores)
  global.error = lapply(global.ace,function(x){x$error})
  global.ace = do.call(rbind,lapply(global.ace,function(x){x$ace}))
  if(add.epsilon){
    global.ace$var = global.ace$var +FMR$params$epsilon/2
    if(numApprox>1&laplace){
      global.error = lapply(global.error,function(x){
        return(data.frame(x=c(x$x,x$x),
                          var=c(x$var+FMR$params$epsilon/2*exp(-1),x$var+FMR$params$epsilon/2*exp(1)),
                          probs=c(x$probs*exp(1)/(1+exp(1)),x$probs/(1+exp(1)))))
      })
    }else{
      for(n in 1:Nnode(FMR$phy)){
        global.error[[n]]$var = global.error[[n]]$var+FMR$params$epsilon/2
      }
    }
  }
  return(list(ace=global.ace,error= global.error))
}
#' @export
fullMarginalReconstructionWithPE = function(phy,x,params,laplace=FALSE,approximate=1,margin=1e-6,numCores=1,asymptotic=5){
  # check if the phylogeny is a rooted, bifurcating tree
  if(!is.binary.phylo(phy)) stop("The phylogeny should be bifurcating!")
  if(!is.rooted(phy)) stop("The phylogeny should be rooted!")
  if(is.null(phy$node.label)) phy = makeNodeLabel(phy)
  x = x[phy$tip.label]
  # get the order of edges from tip to root
  edge.index = reorder.phylo(x = phy,order = "postorder",index.only = TRUE)
  node.order = unique(phy$edge[edge.index,1])
  # parse parameters
  numPE = length(params)/2-1
  if(numPE>0){
    lambdas = params[(1:numPE)*2-1]
    sizes = params[(1:numPE)*2]
  }else{
    lambdas = 0
    sizes = 0
  }
  if(numPE>2){
    # for 3 or more jump processes, run marginal reconstruction in parallel
    # to maximize the number of cores used
    numCores.error = min(numCores,approximate^2)
    numCores.marginal = floor(numCores/numCores.error)
  }else{
    if(numPE>1){
      # for 2 jump processes, run reconstruction in parallel by normal approximations
      # the maximum number of cores used will be (number of normal approximations)^2
      numCores.error = min(numCores,approximate^2)
      numCores.marginal = floor(numCores/numCores.error)
    }else{
      # for 1 jump processes, run reconstruction in the main process to avoid performance issues
      numCores.marginal = 1
      numCores.error = 1
    }
  }
  sigma = params[2*numPE+1]
  epsilon = params[2*numPE+2]
  # initialize tip values
  mrtv = data.frame(focal = 1:Ntip(phy),
                    orientation = sapply(1:Ntip(phy),getLatestAncestor,phy=phy),
                    edge = sapply(1:Ntip(phy),findEdges,phy=phy,ancestor=TRUE),
                    index = 1:Ntip(phy),
                    x = x[1:Ntip(phy)],
                    var=epsilon/2,
                    stringsAsFactors = FALSE)
  # initialize error profile
  error.marginal = lapply(1:Ntip(phy),function(i){
    if(laplace&approximate>1){
      approximate.profile =data.frame(x=unname(x[i]),var=epsilon/2*c(exp(-1),exp(1)),probs=c(exp(1)/(1+exp(1)),1/(1+exp(1))))
    }else{
      approximate.profile = data.frame(x=x[i],var=epsilon/2,probs=1)
    }
    return(list(approximate=approximate.profile))
  })
  # reconstruct marginal ancestral states
  cat("Reconstructing marginal ancestral states...\n")
  for(n in node.order){
    this.match = which(mrtv$orientation==n)
    error.index = do.call(expand.grid,lapply(1:length(this.match),function(i){
      1:dim(error.marginal[[this.match[i]]]$approximate)[1]
    }))
    rec = mclapply(1:dim(error.index)[1],function(i){
      this.rec = marginalReconstructionWithPE(x = sapply(1:length(this.match),function(j){
        error.marginal[[this.match[j]]]$approximate$x[error.index[i,j]]
        }),
                                   l = phy$edge.length[mrtv$edge[this.match]],
                                   epsilons = sapply(1:length(this.match),function(j){
                                     error.marginal[[this.match[j]]]$approximate$var[error.index[i,j]]
                                   }),
                                   sigma = sigma,lambdas=lambdas,sizes = sizes,numCores = numCores.marginal,margin = margin,asymptotic=asymptotic)
      this.probs = prod(sapply(1:length(this.match),function(j){
        error.marginal[[this.match[j]]]$approximate$probs[error.index[i,j]]
      }))
      this.error = this.rec$profile
      this.error$probs = this.error$probs*this.probs
      return(this.error)
    },mc.cores = numCores.error)
    rec = do.call(rbind,rec)
    rec$probs = rec$probs/sum(rec$probs)
    rec.moments = compoundNormalMoments(probs = rec$probs,means = rec$x,vars = rec$var)
    this.edge = findEdges(phy=phy,node = n,ancestor=TRUE)
    if(length(this.edge)<1) this.edge = NA
    this.mrtv = data.frame(focal = n,
                           orientation = getLatestAncestor(phy = phy,x = n),
                           edge = this.edge,
                           index = dim(mrtv)[1]+1,
                           x = unname(rec.moments[1]),
                           var = unname(rec.moments[2]),
                           stringsAsFactors = FALSE)
    mrtv = rbind(mrtv,this.mrtv)
    if(approximate>1){
      this.approximate = findBestNormalApproximation2(probs = rec$probs,means = rec$x,vars = rec$var,numCores = numCores.marginal)
      this.approx.df = data.frame(x=unname(this.approximate$mean),
                                  var=unname(this.approximate$var),
                                  probs=unname(this.approximate$probs))
    }else{
      this.approximate = findBestNormalApproximation1(probs = rec$probs,means = rec$x,vars = rec$var)
      this.approx.df = data.frame(x=unname(this.approximate$mean),
                                  var=unname(this.approximate$var),
                                  probs=unname(this.approximate$probs))
    }
    error.marginal[[this.mrtv$index]] = list(approximate=this.approx.df)
  }
  #
  cat("Reconstructing rerooting marginal ancestral states...\n")
  for(n in rev(node.order)){
    this.edges = findEdges(phy=phy,node=n,ancestor = FALSE)
    for(ed in this.edges){
      this.match = rev(which((mrtv$orientation==n)&(mrtv$edge!=ed)))
      error.index = do.call(expand.grid,lapply(1:length(this.match),function(i){
        1:dim(error.marginal[[this.match[i]]]$approximate)[1]
      }))
      rec = mclapply(1:dim(error.index)[1],function(i){
        this.rec = marginalReconstructionWithPE(x = sapply(1:length(this.match),function(j){
          error.marginal[[this.match[j]]]$approximate$x[error.index[i,j]]
        }),
                                                l = phy$edge.length[mrtv$edge[this.match]],
                                                epsilons = sapply(1:length(this.match),function(j){
                                                  error.marginal[[this.match[j]]]$approximate$var[error.index[i,j]]
                                                }),
                                                sigma = sigma,lambdas=lambdas,sizes = sizes,numCores=numCores.marginal,margin = margin,asymptotic=asymptotic)
        this.probs = prod(sapply(1:length(this.match),function(j){
          error.marginal[[this.match[j]]]$approximate$probs[error.index[i,j]]
        }))
        this.error = this.rec$profile
        this.error$probs = this.error$probs*this.probs
        return(this.error)
      },mc.cores = numCores.error)
      rec = do.call(rbind,rec)
      rec$probs = rec$probs/sum(rec$probs)
      rec.moments = compoundNormalMoments(probs = rec$probs,means = rec$x,vars = rec$var)
      this.mrtv = data.frame(focal = n,
                             orientation = phy$edge[ed,2],
                             edge = ed,
                             index = dim(mrtv)[1]+1,
                             x = unname(rec.moments[1]),
                             var = unname(rec.moments[2]),
                             stringsAsFactors = FALSE)
      mrtv = rbind(mrtv,this.mrtv)
      if(approximate>1){
        this.approximate = findBestNormalApproximation2(probs = rec$probs,means = rec$x,vars = rec$var,numCores = numCores.marginal)
        this.approx.df = data.frame(x=unname(this.approximate$mean),
                                    var=unname(this.approximate$var),
                                    probs=unname(this.approximate$probs))
      }else{
        this.approximate = findBestNormalApproximation1(probs = rec$probs,means = rec$x,vars = rec$var)
        this.approx.df = data.frame(x=unname(this.approximate$mean),
                                    var=unname(this.approximate$var),
                                    probs=unname(this.approximate$probs))
      }
      error.marginal[[this.mrtv$index]] = list(approximate=this.approx.df)
    }
  }
  return(list(phy=phy,marginal = mrtv,error=error.marginal,params=list(lambdas=lambdas,sizes=sizes,sigma=sigma,epsilon=epsilon)))
}
#' @export
discretizeResult = function(res,error=NULL,epsilon=0,laplace=FALSE){
  new.res = lapply(1:dim(res)[1],function(i){
    df = data.frame(node=res$node[i],label=res$label[i],x=round(res$x[i]))
    if(laplace){
      if(is.null(error)){
        df$probs = sapply(res$var,function(x){pconv.norm.laplace(q = 0.5,sigma = x,epsilon = epsilon/2)})-
          sapply(error,function(x){pconv.norm.laplace(q = -0.5,sigma = x,epsilon = epsilon/2)})
      }else{
        df$probs = pMixNormalLaplace(q = round(res$x[i])+0.5,mean = error[[i]]$x,
                              sd = sqrt(error[[i]]$var),laplace.var = epsilon/2,probs = error[[i]]$probs)-
          pMixNormalLaplace(q = round(res$x[i])-0.5,mean = error[[i]]$x,
                     sd = sqrt(error[[i]]$var),laplace.var = epsilon/2,probs = error[[i]]$probs)
      }
    }else{
      if(is.null(error)){
        df$probs = pnorm(q = round(res$x[i])+0.5,mean = res$x[i],sd = sqrt(res$var[i]+epsilon/2))-
          pnorm(q = round(res$x[i])-0.5,mean = res$x[i],sd = sqrt(res$var[i]+epsilon/2))
      }else{
        df$probs = pMixNormal(q = round(res$x[i])+0.5,mean = error[[i]]$x,
                              sd = sqrt(error[[i]]$var+epsilon/2),probs = error[[i]]$probs)-
          pMixNormal(q = round(res$x[i])-0.5,mean = error[[i]]$x,
                     sd = sqrt(error[[i]]$var+epsilon/2),probs = error[[i]]$probs)
      }
    }
    return(df)
  })
  new.res = do.call(rbind,new.res)
  return(new.res)
}
#' @export
convert.pulsR.parameters = function(params){
  new.params = c(params[2],params[c(3,1,10)]^2)*c(1,1,1,2)
  names(new.params) = c("lambda","size","sigma","epsilon")
  return(new.params)
}
