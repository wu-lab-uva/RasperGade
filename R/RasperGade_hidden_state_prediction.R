#' @importFrom castor asr_independent_contrasts
#' @title  Predict hidden state under the BM model
#' @description Predict hidden state under the BM model given a phylogeny, some known states, and the rate of BM
#' @param trait named vector of tip trait values; names should match the tip labels; missing values (NA or not present) are hidden states
#' @param phy phylo-class object from `ape` package
#' @return a list of two elements:
#' @return `$summary` is a data frame listing the means and variances of the hidden states
#' @return `$error` is a list of error distributions where each element is a data frame
#' @export
#' @rdname predict_hidden_state_by_pic
predict_hidden_state_by_pic = function(phy, trait){
  # get names of unknown species
  query.tip = setdiff(phy$tip.label,names(trait[!is.na(trait)]))
  # predict hidden states via rerooting
  res = lapply(1:length(query.tip),function(j){
    this.tree = drop.tip(phy = phy,tip = query.tip[-j])
    this.tree = reroot_tree_at_tip(this.tree,query.tip[j])
    this.pred = asr_independent_contrasts(tree = this.tree,tip_states = trait[this.tree$tip.label],
                                         weighted = TRUE,include_CI = TRUE)
    return(data.frame(node=-1,label=query.tip[j],x=this.asr$ancestral_states[1],
                      var=this.asr$standard_errors[1]^2,stringsAsFactors = FALSE))
  })
  res.error = lapply(res,function(x){
    data.frame(x=x$x,var=x$var,scale=1,prior=1,probs=1)
  })
  return(list(summary=do.call(rbind,res),error = res.error))
}

#' @title  Predict hidden state under the PE model
#' @description Predict hidden state under the PE model given a phylogeny, some known states, and the PE model
#' @param FMR the returned data structure from function `fullMarginalReconstructionWithPE`
#' @param query.keys the tracking key of query tips
#' @param laplace logical, if true, Laplace distribution is used for time-independent variation
#' @param numApprox the number of normal distributions to approximate the Laplace distribution
#' @param margin the total probability mass that the number of jumps omitted in a compound Poisson process
#' @param numCores the number of cores to run in parallel
#' @param asymptotic the threshold of expected number of jumps on a branch beyond which normal distribution is assumed
#' @return a list of two elements:
#' @return `$summary` is a data frame listing the means and variances of the hidden states
#' @return `$error` is a list of error distributions#' @details to be added
#' @export
#' @rdname predictHiddenStateWithPE
predictHiddenStateWithPE = function(FMR,query.keys,laplace=FALSE,numApprox=1,
                                    margin=1e-6,asymptotic=5,numCores=1){
  cv.error = mclapply(query.keys,function(this.key){
    key.check = sapply(this.key$hash,function(x){is.null(FMR$key[[x]])})
    if(all(key.check)){
      warning("Two-side reference mismatch observed, returning NULL.",immediate. = TRUE)
      return(NULL)
    }
    des.nodes = lapply(this.key$hash,function(x){FMR$key[[x]]})
    if(any(key.check)){
      warning("One-side reference mismatch observed, restoring from the other match.",immediate. = TRUE)
      des.nodes[[which(key.check)]] = rev(des.nodes[[which(!key.check)]])
    }
    this.match = sapply(des.nodes,function(des.node){
      mm = which((FMR$marginal$orientation == des.node[2])&
                   (FMR$marginal$focal == des.node[1]))
      if(length(mm)<1){
        mm = which((FMR$marginal$orientation == getRoot(FMR$phy))&
                     (FMR$marginal$focal == des.node[1]))
      }
      return(mm)
    })
    parent.node = intersect(des.nodes[[1]],
                            sapply(des.nodes[[2]],getAncestor,phy=FMR$phy))
    if(length(parent.node)<1) parent.node = getAncestor(phy = FMR$phy,des.nodes[[1]][1])
    this.scale = FMR$scale[parent.node]
    this.epsilon = FMR$epsilon.branch[parent.node]
    this.key$l = this.key$l/this.scale+c(0,0,this.epsilon)
    this.l = this.key$l[1:2]
    error.index = do.call(expand.grid, lapply(1:length(this.match), 
                                              function(i) {
                                                1:dim(FMR$error[[this.match[i]]]$approximate)[1]
                                              }))
    rec = do.call(rbind, lapply(1:dim(error.index)[1], function(i) {
      this.rec = marginalReconstructionWithPE(x = sapply(1:length(this.match), 
                                                         function(j) {
                                                           FMR$error[[this.match[j]]]$approximate$x[error.index[i,j]]
                                                         }), l = this.l, 
                                              epsilons = sapply(1:length(this.match), function(j) {
                                                FMR$error[[this.match[j]]]$approximate$var[error.index[i,j]]
                                              }), sigma = FMR$params$sigma, lambdas = FMR$params$lambdas, 
                                              sizes = FMR$params$sizes, margin = margin, asymptotic = asymptotic)
      this.probs = prod(sapply(1:length(this.match), function(j) {
        FMR$error[[this.match[j]]]$approximate$probs[error.index[i,j]]
      }))
      this.error = this.rec$profile
      this.error$probs = this.error$probs * this.probs
      return(this.error)
    }))
    rec$probs = rec$probs/sum(rec$probs)
    rec.moments = compoundNormalMoments(probs = rec$probs, 
                                        means = rec$x, vars = rec$var)
    this.hsp = data.frame(node = -1, label = this.key$label, 
                          x = unname(rec.moments[1]), var = unname(rec.moments[2]), 
                          stringsAsFactors = FALSE)
    rec = do.call(rbind, lapply(1, function(i) {
      this.rec = marginalReconstructionWithPE(x = this.hsp$x,
                                              l = this.key$l[3], 
                                              epsilons = this.hsp$var,
                                              sigma = FMR$params$sigma, lambdas = FMR$params$lambdas, 
                                              sizes = FMR$params$sizes, margin = margin, asymptotic = asymptotic)
      this.error = this.rec$profile
      this.error$probs = this.error$probs
      return(this.error)
    }))
    rec$probs = rec$probs/sum(rec$probs)
    rec.moments = compoundNormalMoments(probs = rec$probs, 
                                        means = rec$x, vars = rec$var)
    this.hsp = data.frame(node = -1, label = this.key$label, 
                          x = unname(rec.moments[1]), var = unname(rec.moments[2]), 
                          stringsAsFactors = FALSE)
    return(list(error = rec, hsp = this.hsp))
  },mc.cores = numCores)
  cv = do.call(rbind, lapply(cv.error, function(x) {
    x$hsp
  }))
  cv.error = lapply(cv.error, function(x) {
    x$error
  })
  cv$var = cv$var + FMR$params$epsilon/2
  for (ii in 1:length(cv.error)) {
    cv.error[[ii]]$var = cv.error[[ii]]$var + FMR$params$epsilon/2
  }
  return(list(hsp = cv, error = cv.error))
}

