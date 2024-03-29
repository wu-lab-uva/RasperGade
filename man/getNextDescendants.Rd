% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ape_utils.R
\name{getNextDescendants}
\alias{getNextDescendants}
\alias{get_descendant_tips_for_each_node}
\alias{get_descendant_nodes_for_each_node}
\title{Find the descendants of a node}
\usage{
getNextDescendants(phy, x)

get_descendant_tips_for_each_node(phy)

get_descendant_nodes_for_each_node(phy)
}
\arguments{
\item{phy}{a phylo-class object from the \code{ape} package}

\item{x}{the index of the focal node}
}
\value{
\code{getNextDescendants} returns the indices of the immediate descendants. If a tip is supplied, `NA`` will be returned.

\code{get_descendant_tips_for_each_node} and \code{get_descendant_nodes_for_each_node} returns a list of descendants.

The list is ordered so that the index matches the tip/node index in the phylo-class object.
}
\description{
Utility functions to find the descendants of a node
}
\seealso{
[getRoot][getAncestor][findSisters][getConnected]
}
