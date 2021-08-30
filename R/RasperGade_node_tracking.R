#' @title  Calculating the MD5 sum of an unordered character vector
#'
#' @description  Calculate the MD5 sum of an unordered character vector
#'
#' @export
#' @rdname sort_and_digest
sort_and_digest = function(str){
  digest::digest(object = paste0(sort(str),collapse=""))
}

#' @title  Calculating traceable node labels
#'
#' @description  Calculate a traceable node label that can be maintained among trees
#'
#' @export
#' @rdname make_traceable_node_labels
make_traceable_edge_label = function(phy){
  edge.hash = new.env(hash = TRUE)
  node.tips = get_descendant_tips_for_each_node(phy)
  edge.key = lapply(1:length(phy$edge.length),function(i){
    des.idx = phy$edge[i,2]
    tip.label.set1 = node.tips[[des.idx]]
    tip.label.set2 = setdiff(phy$tip.label,tip.label.set1)
    node.key1 = sort_and_digest(tip.label.set1)
    node.key2 = sort_and_digest(tip.label.set2)
    return(c(node.key1,node.key2))
  })
  for(i in 1:length(phy$edge.length)){
    if(is.rooted.phylo(phy)&(phy$edge[i,1]==getRoot(phy))){
      edge.hash[[edge.key[[i]][1]]] = c(phy$edge[i,2],findSisters(phy,phy$edge[i,2]))
      edge.hash[[edge.key[[i]][2]]] = rev(edge.hash[[edge.key[[i]][1]]])
    }else{
      edge.hash[[edge.key[[i]][1]]] = c(phy$edge[i,2],phy$edge[i,1])
      edge.hash[[edge.key[[i]][2]]] = c(phy$edge[i,1],phy$edge[i,2])
    }
  }
  return(edge.hash)
}

#' @title  Partition test set from reference data
#'
#' @description  Prepare the reference data for cross-validation at different degree of uncertainty
#'
#' @export
#' @rdname partition_test_set_with_minimum_distance
partition_test_set_with_minimum_distance = function(phy,trait,test.set,distance=0,location=TRUE){
  if(is.character(test.set)) test.set = match(test.set,phy$tip.label)
  train.set = setdiff(1:Ntip(phy),test.set)
  set.dist = castor::find_nearest_tips(tree = phy,only_descending_tips = FALSE,
                               target_tips = test.set,as_edge_counts = FALSE)
  close.ref = which(set.dist$nearest_distance_per_tip<=distance)
  new.phy = drop.tip(phy,union(test.set,close.ref))
  new.insert.phy = drop.tip(phy,setdiff(close.ref,test.set))
  new.query.hash = list()
  if(location) new.query.hash = locate_query_sequences(new.insert.phy,phy$tip.label[test.set])
  return(list(ref=list(phy=new.phy,dat=trait[new.phy$tip.label]),
              test=phy$tip.label[test.set],phy=new.insert.phy,
              hash=new.query.hash$hash))
}

#' @title  Locate query sequences
#'
#' @description  Find the position of insertion in a phylogeny with query tips inserted
#'
#' @export
#' @rdname locate_query_sequences
locate_query_sequences = function(query.phy,query.tips){
  if(is.character(query.tips)) query.tips = match(query.tips,query.phy$tip.label)
  query.location = lapply(1:length(query.tips),function(i){
    new.query.phy = unroot(drop.tip(phy = query.phy,tip = query.tips[-i]))
    new.query.phy = root.phylo(phy = new.query.phy,
                               outgroup = which(new.query.phy$tip.label==query.phy$tip.label[query.tips[i]]),
                               resolve.root = FALSE)
    tip.idx = which(new.query.phy$tip.label==query.phy$tip.label[query.tips[i]])
    insert.idx = getLatestAncestor(new.query.phy,tip.idx)
    this.des = setdiff(getNextDescendants(new.query.phy,insert.idx),tip.idx)
    tips.per.nodes = get_descendant_tips_for_each_node(new.query.phy)
    des.hash = sapply(this.des,function(x){sort_and_digest(tips.per.nodes[[x]])})
    des.l = sapply(this.des,function(x){
      new.query.phy$edge.length[which(new.query.phy$edge[,2]==x)]
    })
    leading.l = new.query.phy$edge.length[which(new.query.phy$edge[,2]==tip.idx)]
    return(list(hash=des.hash,
                l=c(des.l,leading.l),
                label=query.phy$tip.label[query.tips[i]]))
  })
  return(list(query=query.phy$tip.label[query.tips],hash=query.location))
}