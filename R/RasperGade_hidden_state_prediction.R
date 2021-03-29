#' @title  Predict hidden state under the BM model
#'
#' @description Reconstruct marginal ancestral state using immediate descendants under BM or pulsed evolution model
#'
#' @param trait named vector of tip trait values; names should match the tip labels; hidden states to be predicted have NA values
#' @param phy phylo-class object from `ape` package
#' @param rate the rate of BM
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