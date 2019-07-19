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
simulateRateInPhylogeny = function(phy,model,rate,probs){
  res.phy = phy
  if(is.null(model)){
    multipliers = sample(x = rate,size = length(phy$edge.length),replace = TRUE,prob = probs)
  }else{
    func = match.fun(paste0("r",model$distr))
    multipliers = do.call(func,c(list(n=length(phy$edge.length)),model$fullparams[model$params]))
  }
  if(model$homogenous.sisters){
    for(n in Ntip(phy)+1:Nnode(phy)){
      res.phy$edge.length[res.phy$edge[,1]==n] = phy$edge.length[phy$edge[,1]==n]*multipliers[n-Ntip(phy)]*
        do.call(model$scalefunc,c(list(l=sum(phy$edge.length[phy$edge[,1]==n])),model$fullparams))
    }
  }else{
    res.phy$edge.length = phy$edge.length*multipliers*do.call(model$scalefunc,c(list(l=phy$edge.length),model$fullparams))
  }
  return(res.phy)
}
#
simulateTrait = function(phy,start=0,BM.tree = FALSE,...){
  res.phy = simulateRateInPhylogeny(phy,...)
  trait = simulateBMTrait(res.phy,start)
  if(BM.tree){
    return(list(trait=trait,BM.phy=BM.phy))
  }else{
    return(trait)
  }
}
