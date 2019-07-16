#
getPIC = function(phy,x,remove.zero=FALSE,tip.only=TRUE){
  if(is.null(phy$node.label)) phy = makeNodeLabel(phy)
  x = x[phy$tip.label]
  tip.node = setdiff(sapply(1:Ntip(phy),getLatestAncestor,phy=phy),sapply(Ntip(phy)+1:phy$Nnode,getLatestAncestor,phy=phy))
  term.l = lapply(tip.node,function(n){phy$edge.length[phy$edge[,1]==n]})
  term.x = lapply(tip.node,function(n){x[phy$edge[phy$edge[,1]==n,2]]})
  phy.pic = pic(x = x,phy = phy)
  term.pic = phy.pic[tip.node-Ntip(phy)]
  if(remove.zero){
    p0 = sum(term.pic==0)/length(term.pic)
    idx = term.pic!=0
    return(list(pic=term.pic[idx],x=term.x[idx],l=term.l[idx],node=tip.node[idx],p0=p0))
  }else{
    return(list(pic=term.pic,x=term.x,l=term.l,node=tip.node))
  }
}
#
