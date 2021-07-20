#' @title  Predict hidden state under the BM model
#'
#' @description Predict hidden state under the BM model given a phylogeny, some known states, and the rate of BM
#'
#' @param trait named vector of tip trait values; names should match the tip labels; missing values are hidden states
#' @param phy phylo-class object from `ape` package
#' @param rate the rate of BM
#' @return a data frame listing the means and variances of the hidden states
#' @details this function is modified from the function in `picante` package
#' @export
#' @rdname predTraitByPIC
predTraitByPIC = function(phy, trait,rate=1){
  # get names of known and unknown species
  refSpecies = names(trait)
  tarSpecies = phy$tip.label[is.na(match(phy$tip.label,refSpecies))]
  #
  pred = as.data.frame(matrix(nrow=length(tarSpecies), ncol=4, dimnames=list(tarSpecies, c("node","label","x","var"))))
  #
  for (i in tarSpecies) {
    # keep only the focal target species and the references
    this.tree = drop.tip(phy = phy,tip = tarSpecies[tarSpecies!=i])
    # reroot the tree at the focal target species
    this.tree = root(this.tree, i, resolve.root=FALSE)
    # get branch length leading to the focal species in the rerooted tree for standard errors
    focal.bl = this.tree$edge.length[which(this.tree$edge[,2]==which(this.tree$tip.label==i))]
    # trim the focal species
    this.tree = drop.tip(this.tree, i)
    # get the adjusted branch lengths for the root estimate
    pic.tree = pic(x = trait,phy = this.tree,rescaled.tree = TRUE)$rescaled.tree
    root.bl = pic.tree$edge.length[pic.tree$edge[,1]==getRoot(pic.tree)]
    # use PIC (BM-ML analytic solution) to estimate trait value and error
    est = ace(x = trait[this.tree$tip.label], phy = this.tree, method="pic",type = "continuous")$ace[1]
    this.var = 1/sum(1/root.bl) + focal.bl
    pred[i,] <- data.frame(node=which(phy$tip.label==i),label=i,x=unname(est), var=unname(this.var))
  }
  # adjust for rate of BM
  pred$var = pred$var*rate
  return(pred)
}

#' @export
#' @rdname predictHiddenStateWithPE
predictHiddenStateWithPE = function(FMR,query.keys,laplace=FALSE,numApprox=1,
                                    margin=1e-6,asymptotic=5,numCores=1){
  cv.error = lapply(query.keys,function(this.key){
    des.node = FMR$key[[this.key$hash[1]]]
    if(is.null(des.node)){
      this.match = which(is.na(FMR$marginal$orientation)&(FMR$marginal$focal == getRoot(FMR$phy)))
      this.l = this.key$l[2]
    }else{
      parent.node = getLatestAncestor(FMR$phy,des.node)
      this.match = c(which((FMR$marginal$orientation == des.node)&(FMR$marginal$focal == parent.node)),
                     which((FMR$marginal$focal == des.node)&(FMR$marginal$orientation == parent.node)))
      this.scale = FMR$phy$edge.length[FMR$phy$edge[,2]==des.node]/
        sum(this.key$l[1:2])
      this.key$l = this.key$l*this.scale
      this.l = this.key$l[1:2]
    }
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
  })
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
