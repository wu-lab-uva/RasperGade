#' @title  Reconstruct marginal ancestral state from a node's immediate descendants
#'
#' @description Reconstruct marginal ancestral state using immediate descendants under BM or pulsed evolution model
#'
#' @param x vector of means of mixed normal distributions
#' @param l vector of standard deviations of mixed normal distributions
#' @param lambdas vector of jump frequencies, has the same length as the number of compound Poisson processes
#' @param sizes vector of variances of trait change per jump, has the same length as the number of compound Poisson processes
#' @param sigma rate of BM, has a length of 1
#' @param epsilons vector of uncertainty not explained by the branches, has the same length as the number of immediate descendants
#' @param margin total probability mass of the number of jump not taken into account per compound Poisson process
#' @param numCores number of cores to use
#' @param asymptotic threshold of expected number of jumps beyond which the trait change is assumed to be normal
#' @return `marginalReconstruction` returns a vector of the mean variance and the scaling factor
#' @return `marginalReconstructionWithPE` returns a list of the mean, variance and a data frame of the normal mixture
#' @export
#' @rdname marginalReconstruction
marginalReconstruction = function(x,l){
  xs = productNormalDensity(x,l)
  return(c(x=unname(xs["mean"]),var=unname(xs["var"]),scale=unname(xs["scale"])))
}

#' @export
#' @rdname marginalReconstruction
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

#' @title  Reconstruct all marginal ancestral state in a phylogeny
#'
#' @description Reconstruct all marginal ancestral states, including those require rerooting, in a phylogeny
#'
#' @param x named vector of tip trait values; names should match the tip labels; NA should not be included
#' @param phy phylo-class object from `ape` package
#' @param params a vector of pulsed evolution parameters
#' @param laplace logical, when true, time-independent variation follows Laplace distribution
#' @param approximate number of normal distributions to approximate the uncertainty distribution for internal nodes
#' @param margin total probability mass of the number of jump not taken into account per compound Poisson process
#' @param numCores number of cores to use
#' @param asymptotic threshold of expected number of jumps beyond which the trait change is assumed to be normal
#' @details
#' Parameter vector should be in the format of lambda,size,...,sigma,epsilon
#' Currently `laplace` should always set to be `FALSE`.
#' @return returns a list that contains all the information to calculate the global ancestral states and hidden states
#' @export
#' @rdname fullMarginalReconstructionWithPE
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

