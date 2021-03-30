#' @title  Simulate trait under pulsed evolution
#'
#' @description  Simulate trait under pulsed evolution
#'
#' @param phy a phylo-class object from the `ape` package
#' @param l branch length of time/divergence
#' @param lambda jump frequency
#' @param size variance of trait change per jump
#' @param sigma rate of BM
#' @param epsilon variance of time-independent variation
#' @export
#' @rdname simulatePET
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
  ebm =  rnorm(n=n)*sqrt(epsilon)
  return(list(l = l, t = t,nj=nj,
              dxbm = dxbm, dxpe = dxpe,
              ebm = ebm))
}

#' @export
#' @rdname simulatePET
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
    if(phy$edge[l,2]<=Ntip(phy)) trait.sim[phy$edge[l,2]] = trait.sim[phy$edge[l,2]] + trait.delta$ebm[l]
  }
  return(trait.sim)
}
