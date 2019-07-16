#
getRoot = function(phy){
  return(Ntip(phy) + 1)
}
#
getLatestAncestor = function (phy, x){
  # x should be the index of the focal node, not the node/tip label
  if (x == Ntip(phy) + 1) return(NA)
  i = which(phy$edge[, 2] == x)
  return(phy$edge[i, 1])
}
#
getAncestor = function (phy, x,level=1,trace=FALSE){
  # returning the ancestor of a specified degree, compatiable with other functions with the same name
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
getNextDescendants = function (phy, x){
  # x should be the index
  if (x <= Ntip(phy)) 
    return(NA)
  i <- which(phy$edge[, 1] == x)
  return(phy$edge[i, 2])
}
#
getSisters <- function(phy,x){
  # x should be the index
  parent = getLatestAncestor(phy,x)
  if(is.na(parent)) return(NA)
  twin = getNextDescendants(phy,parent)
  sis = twin[twin!=x]
  return(sis)
}
#
getEdges = function(phy,node,ancestor = FALSE){
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
