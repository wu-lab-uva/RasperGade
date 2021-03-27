#' @import ape
#' @title  Find the ancestor of a node or tip
#'
#' @description  Utility functions to find the ancestors of a node or tip
#'
#' @param phy a phylo-class object from the `ape` package
#' @param x the index of the focal node/tip
#' @param level the depth between the ancestor and the focal node/tip
#' @param trace logical, if FALSE (the default), only the earliest ancestor is returned, otherwise, all ancestors along the path are returned
#' @return A vector of ancestor indices is returned. If root is supplied, `NA` will be returned
#' @seealso [getRoot][getNextDescendants][findSisters][getConnected]
#' @export
#' @rdname getAncestor
getLatestAncestor = function (phy, x){
  if (x == Ntip(phy) + 1) return(NA)
  i <- which(phy$edge[, 2] == x)
  return(phy$edge[i, 1])
}

#' @export
#' @rdname getAncestor
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

#' @title  Find the root of a phylogeny
#'
#' @description  An utility function to find the root of a tree
#'
#' @param phy a phylo-class object from the `ape` package
#' @return The index of the root
#' @seealso [getAncestor][getNextDescendants][findSisters][getConnected]
#' @export
#' @rdname getRoot
getRoot = function(phy){
  return(Ntip(phy) + 1)
}

#' @title  Find the immediate descendants of a node
#'
#' @description  An utility function to find the immediate descendants of a node
#'
#' @param phy a phylo-class object from the `ape` package
#' @param x the index of the focal node
#' @return The indices of the immediate descendants. If a tip is supplied, `NA`` will be returned.
#' @seealso [getRoot][getAncestor][findSisters][getConnected]
#' @export
#' @rdname getNextDescendants
getNextDescendants = function (phy, x){
  if (x <= Ntip(phy))
    return(NA)
  i <- which(phy$edge[, 1] == x)
  return(phy$edge[i, 2])
}

#' @title  Find the sister nodes/tips of a node/tip
#'
#' @description  An utility function to find the sister nodes/tips of a node/tip
#'
#' @param phy a phylo-class object from the `ape` package
#' @param x the index of the focal node
#' @return The indices of the sisters. If the root is supplied, `NA`` will be returned.
#' @return If a non-root node has no sister, an empty numeric vector will be returned.
#' @seealso [getRoot][getAncestor][getNextDescendants][getConnected]
#' @export
#' @rdname findSisters
findSisters <- function(phy,x){
  # x should be the index
  parent = getLatestAncestor(phy,x)
  if(is.na(parent)) return(NA)
  twin = getNextDescendants(phy,parent)
  sis = twin[twin!=x]
  return(sis)
}

#' @title  Find the edges/nodes connected to a node
#'
#' @description  An utility function to find the sister nodes/tips of a node/tip
#'
#' @param phy a phylo-class object from the `ape` package
#' @param x the index of the focal node
#' @param node the index of the focal node
#' @param ancestor logical, if FALSE (the default), the edges originating from the node are returned, otherwise, the edge leading to the node is returned.
#' @return `findEdges` returns the indices of the edges. Note that these indices are different from those of nodes/tips
#' @return `getConnected` returns the indices of the connected nodes.
#' @seealso [getRoot][getAncestor][getNextDescendants][findSisters]
#' @export
#' @rdname getConnected
findEdges = function(phy,node,ancestor = FALSE){
  # node should be the index
  if(ancestor){
    edge.index = which(phy$edge[,2]==node)
  }else{
    edge.index = which(phy$edge[,1]==node)
  }
  return(edge.index)
}

#' @export
#' @rdname getConnected
getConnected = function(phy,x){
  conn = c(phy$edge[phy$edge[,2]==x,1],phy$edge[phy$edge[,1]==x,2])
  return(conn)
}
