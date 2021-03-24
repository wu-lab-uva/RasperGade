# This script contains some useful function extensions for tree manipulation and analyses from the ape package.
# By Yingnan Gao, 20210117;  Email: yg5ap@virginia.edu

# load required library
library(ape)
#
getLatestAncestor = function (phy, x){
  # x should be the index of the focal node, not the node/tip label
  if (x == Ntip(phy) + 1) return(NA)
  i <- which(phy$edge[, 2] == x)
  return(phy$edge[i, 1])
}
#
getAncestor = function (phy, x,level=1,trace=FALSE){
  # returning the ancestor of a specified degree, compatible with other functions with the same name
  # e.g. Kemble's 2012 paper
  if(trace) all.anc = numeric(level)
  anc = x
  for(i in 1:level){
    anc = getLatestAncestor(phy=phy,x=anc)
    if(trace) all.anc[i]=anc
    if(is.na(anc)) break
  }
  if(trace){
    return(all.anc[-which(all.anc==0)])
  }else{
    return(anc) 
    }
}
#
getRoot = function(phy){
  return(Ntip(phy) + 1)
}
#
getNextDescendants = function (phy, x){
  # x should be the index
  if (x <= Ntip(phy)) 
    return(NA)
  i <- which(phy$edge[, 1] == x)
  return(phy$edge[i, 2])
}
#
findSisters <- function(phy,x){
  # x should be the index
  parent = getLatestAncestor(phy,x)
  if(is.na(parent)) return(NA)
  twin = getNextDescendants(phy,parent)
  sis = twin[twin!=x]
  return(sis)
}
#
findEdges = function(phy,node,ancestor = FALSE){
  # node should be the index
  if(ancestor){
    edge.index = which(phy$edge[,2]==node)
  }else{
    edge.index = which(phy$edge[,1]==node)
  }
  return(edge.index)
}
#
getConnected = function(phy,x){
  conn = c(phy$edge[phy$edge[,2]==x,1],phy$edge[phy$edge[,1]==x,2])
  return(conn)
}
#
predTraitByPIC = function(phy, trait,rate=1){
  # predicts trait value from known extant species using PIC (BM-ML)
  # this function is modified from Kembel's function in picante package
  #
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
#
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
#