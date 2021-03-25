# This script contains some useful function extensions for tree manipulation and analyses from the ape package.
# By Yingnan Gao, 20210325;  Email: yg5ap@virginia.edu

# load required library
library(ape)

# function to get the immediate ancestor of a node/tip
## phy is a phylo-class object
## x is the index of the tip/node whose immediate ancestor's index will be returned
## output value of this function is the index of the immediate ancestor
## if the root is supplied as x, NA will be returned
### note that this function is NOT vectorized
getLatestAncestor = function (phy, x){
  if (x == Ntip(phy) + 1) return(NA)
  i <- which(phy$edge[, 2] == x)
  return(phy$edge[i, 1])
}

# function to get the ancestor of a node/tip with a specified depth
## phy is a phylo-class object
## x is the index of the tip/node whose immediate ancestor's index will be returned
## level is the depth between the ancestor and the focal node
## trace determines if all indices along the path should be returned
## output value of this function is a vector of the index/indices of the ancestor(s)
### this function is NOT vectorized
### this function should be compatible with other functions with the same name (e.g., Kemble et al 2012)
getAncestor = function (phy, x,level=1,trace=FALSE){
  if(trace) all.anc = numeric(level)
  anc = x
  for(i in 1:level){
    anc = getLatestAncestor(phy=phy,x=anc)
    if(trace) all.anc[i]=anc
    if(is.na(anc)) break
  }
  if(trace){
    return(all.anc[which((all.anc!=0)&(!is.na(all.anc)))])
  }else{
    return(anc) 
  }
}

# function to get the root of a rooted phylogeny
## phy is a phylo-class object
## output value of this function is the index of the root
### for unrooted tree, this function will also return an index, but it is not the true root of the tree
getRoot = function(phy){
  return(Ntip(phy) + 1)
}

# function to get the immediate descendants of a node/tip
## phy is a phylo-class object
## x is the index of the tip/node whose immediate descendants' indices will be returned
## output value of this function is a vector of the indices of the immediate descendants
## if a tip is supplied as x, NA will be returned
### note that this function is NOT vectorized
getNextDescendants = function (phy, x){
  if (x <= Ntip(phy))
    return(NA)
  i <- which(phy$edge[, 1] == x)
  return(phy$edge[i, 2])
}

# function to get the sister node/tip of a node/tip
## phy is a phylo-class object
## x is the index of the tip/node whose sisters' indices will be returned
## output value of this function is a vector of the indices of the sisters
## if the root is supplied as x, NA will be returned
### note that this function is NOT vectorized
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
