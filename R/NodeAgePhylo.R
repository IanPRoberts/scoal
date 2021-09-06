#' Phylo object from node ages
#'
#' Produces a phylo object under a structured coalescent model using the age of
#' each node
#'
#' @param edge Edge matrix for a structured coalescent tree
#' @param edge.deme Label for the deme in which each edge lies
#' @param node.ages Time of each node in the coalescent tree
#' @param tip.label (optional) Label for each tip; if not entered, labeled 1:n.leaf
#' @param plot.phylo logical; if FALSE (default) plot is not produced
#'
#' @return An object of class \code{phylo} (from package \code{ape}) storing the given structured coalescent process
#'
#' @export

node.age.phylo <- function(edge, edge.deme, node.ages, tip.label = NA, plot.phylo = FALSE){
  all.nodes <- sort(unique(as.vector(edge)))
  node.freq <- table(match(as.vector(edge),all.nodes))  #Frequency of each node featured in edge matrix

  leaf.nodes <- all.nodes[node.freq == 1]  #Leaf nodes have frequency 1
  n.leaf <- length(leaf.nodes)

  Nnode <- dim(edge)[1] - n.leaf + 1  #Number of interior nodes

  edge.length <- abs(node.ages[edge[,1]] - node.ages[edge[,2]])  #Difference in node age between ends of each edge

  if(is.na(tip.label) == TRUE){
    tip.label <- 1:n.leaf
  }

  phylo <- list()
  class(phylo) <- 'phylo'
  phylo$edge <- edge
  phylo$edge.length <- edge.length
  phylo$Nnode <- Nnode
  phylo$tip.label <- tip.label

  if (plot.phylo == TRUE){
    n.deme <- max(edge.deme)
    color.palette <- rainbow(n.deme)
    edge.color <- color.palette[edge.deme]
    plot(phylo, edge.color = edge.color)
    axisPhylo(1, root.time = node.ages[n.leaf+1], backward = FALSE)
  }
  return(phylo)
}
