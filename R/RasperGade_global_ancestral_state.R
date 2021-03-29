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
continuousAceByPIC = function(x,phy,rate=1){
  # this function calculates the global ancestral state of all internal nodes by rerooting the tree
  # it is the BM-ML solution to all ancestral states of internal nodes
  #
  # initialize variables and node labels
  res = numeric(phy$Nnode)
  if(is.null(phy$node.label)){
    phy$node.label = as.character(Ntip(phy)+1:phy$Nnode)
    names(res) = phy$node.label
  }else{
    names(res) = phy$node.label
  }
  this.var = res
  # the root ancestral state is the same as the marginal ancestral state from ace function
  res[1] = ace(phy = phy,x = x[phy$tip.label],type = "continuous",method = "pic")$ace[1]
  pic.phy = pic(x = x[phy$tip.label],phy = phy,rescaled.tree = TRUE)[[2]]
  this.var[1] = 1/sum(1/pic.phy$edge.length[pic.phy$edge[,1]==(Ntip(pic.phy)+1)])
  # for other internal nodes, reroot the tree to get the global ancestral states
  for(i in phy$node.label[-1]){
    new.phy = root(phy,node=which(phy$node.label==i)+Ntip(phy),resolve.root=TRUE)
    res[i] = ace(phy = new.phy,x = x[new.phy$tip.label],type = "continuous",method = "pic")$ace[1]
    pic.phy = pic(x = x[new.phy$tip.label],phy = new.phy,rescaled.tree = TRUE)[[2]]
    this.var[i] = 1/sum(1/pic.phy$edge.length[pic.phy$edge[,1]==(Ntip(pic.phy)+1)])
  }
  # adjust for rate of BM
  ace.df = data.frame(node = Ntip(phy)+1:Nnode(phy),label=phy$node.label,x=res,var=this.var*rate,row.names = phy$node.label)
  return(ace.df)
}
