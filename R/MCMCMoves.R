#' Migration Pair Birth MCMC Move
#'
#' Performs a Migration pair birth move (Ewing et al. 2004). Adds two migration
#' nodes on an edge selected uniformly at random from a structured coalescent
#' process, allocating a deme for the added edge such that a migration event
#' does not target its origin deme.
#'
#' @param edge Edge matrix for a structured coalescent tree
#' @param edge.deme Label for the deme in which each edge lies
#' @param node.ages Time of each node in the coalescent tree
#' @param n.deme Number of possible demes in the process
#'
#' @return List with 3 elements, first element the edge matrix, second element
#' the node ages and third element giving the edge demes

mig.pair.birth <- function(edge, edge.deme, node.ages, n.deme){
  all.nodes <- (unique(as.vector(edge)))

  #Update edge matrix
  selected.edge <- sample(seq_along(edge[,1]),1)  #Samples 1 row of edge matrix uniformly at random
  new.nodes <- max(all.nodes) + c(1,2)  #New migration node IDs to add
  new.node.ages <- sort(runif(2, min = min(node.ages[edge[selected.edge,]]), max = max(node.ages[edge[selected.edge,]]))) #Node ages distributed uniformly along selected edge
  edge <- rbind(edge, as.vector(new.nodes), c(edge[selected.edge,1],new.nodes[1]))  #Add new edges to edge matrix
  edge[selected.edge,1] <- new.nodes[2]  #Update remaining edge endpoints


  #Update edge deme labels
  selected.deme <- sample.vector((1:n.deme)[-edge.deme[selected.edge]],1)
  edge.deme <- c(edge.deme, selected.deme, edge.deme[selected.edge])

  #Update node ages
  node.ages <- c(node.ages, new.node.ages)

  output <- list()
  output$edge <- edge
  output$node.ages <- node.ages
  output$edge.deme <- edge.deme

  return(output)
}

#' Migration Pair Death MCMC Move
#'
#' Performs a Migration pair death move (Ewing et al. 2004). Deletes two
#' migration nodes if they lie between two edges in the same deme.
#'
#' @param edge Edge matrix for a structured coalescent tree
#' @param edge.deme Label for the deme in which each edge lies
#' @param node.ages Time of each node in the coalescent tree
#' @param n.deme Number of possible demes in the process
#'
#' @return List with 3 elements, first element the edge matrix, second element
#' the node ages and third element giving the edge demes

mig.pair.death <- function(edge, edge.deme, node.ages, n.deme){
  # Identify node types based on frequencies in edge matrix
    # Leaves have freq 1
    # Root and migration nodes have freq 2 (Root is n.leaf + 1)
    # Coal nodes have freq 3

  all.nodes <- sort(unique(as.vector(edge)))
  node.freq <- table(match(as.vector(edge),all.nodes))  #Frequency of each node featured in edge matrix

  leaf.nodes <- all.nodes[node.freq == 1]  #Leaf nodes have frequency 1
  n.leaf <- length(leaf.nodes)
  root.node <- n.leaf + 1  #Root node is n.leaf + 1 for ape phylo-style object
  coalescence.nodes <- all.nodes[node.freq == 3]  #Coalescence nodes have frequency 3
  migration.nodes <- all.nodes[! all.nodes %in% c(root.node, coalescence.nodes, leaf.nodes)]

  selected.edge <- sample(seq_along(edge[,1]),1)  #Samples 1 row of edge matrix uniformly at random

  if ((edge[selected.edge,1] %in% migration.nodes == TRUE) && (edge[selected.edge,1] %in% migration.nodes == TRUE )){
    #DO EVENT
  } else{
    "REJECT"
    #return("REJECT")
  }

}
