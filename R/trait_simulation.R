#
simulateBMTrait = function(phy,start=0){
  trait.sim = numeric(length = Ntip(phy)+Nnode(phy))+start
  if(is.null(phy$node.label)){
    names(trait.sim) = c(phy$tip.label,as.character(Ntip(phy)+(1:phy$Nnode)))
  }else{
    names(trait.sim) = c(phy$tip.label,phy$node.label)
  }
  phy = reorder.phylo(phy,order = "postorder")
  trait.delta = rnorm(length(phy$edge.length))*sqrt(phy$edge.length)
  for(i in 1:length(phy$edge.length)){
    l = length(phy$edge.length)+1-i
    trait.sim[phy$edge[l,2]] = trait.sim[phy$edge[l,1]] + trait.delta[l]
  }
  return(trait.sim)
}
#
