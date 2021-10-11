#' Vector sampling
#'
#' Samples a vector in the same way as base::sample except if the vector has
#' length 1, the value entered is always returned
#'
#' @inheritParams base::sample
#'
#' @return A vector of length \code{size} with elements drawn from \code{x}
#'
#' @export

sample.vector <- function(x, size, replace = FALSE, prob = NULL){
  #Samples a vector. If the vector has length 1, always returns the single value
  if (length(x) == 1){
    x
  }
  else{
    sample(x, size, replace, prob)
  }
}



#' Parent node
#'
#' Given a node from a phylo object, identifies the node ID of the parent node
#'
#' @param node Node to find parent node of
#' @param edge Edge matrix from a \code{phylo} object
#' @param node.ages Time of each node in the \code{phylo} object
#'
#' @return Node ID of the parent node. If the input node is the root, returns Inf
#'
#' @export

parent.node <- function(node, edge, node.ages){
  all.nodes <- sort(unique(as.vector(edge)))
  node.freq <- table(match(as.vector(edge),all.nodes))  #Frequency of each node featured in edge matrix

  leaf.nodes <- all.nodes[node.freq == 1]  #Leaf nodes have frequency 1
  n.leaf <- length(leaf.nodes)
  root.node <- n.leaf + 1  #Root node is n.leaf + 1 for ape phylo-style object

  if (node == root.node){
    return(Inf)
  }
  else{
    possible.parent <- c(edge[which(edge[,1] == node),2], edge[which(edge[,2] == node),1])
    parent <- possible.parent[node.ages[possible.parent] < node.ages[node]]
    return(parent)
  }
}



#' Child nodes
#'
#' Given a node from a \code{phylo} object, identifies the node IDs of any direct
#' child nodes
#'
#' @inheritParams parent.node
#'
#' @return Vector of node IDs of any child nodes. If the input node is a leaf, returns NA
#'
#' @export

child.nodes <- function(node, edge, node.ages){
  all.nodes <- sort(unique(as.vector(edge)))
  node.freq <- table(match(as.vector(edge),all.nodes))  #Frequency of each node featured in edge matrix

  leaf.nodes <- all.nodes[node.freq == 1]  #Leaf nodes have frequency 1

  if (node %in% leaf.nodes){
    return(NA)
  }
  else{
    possible.children <- c(edge[which(edge[,1] == node),2], edge[which(edge[,2] == node),1])
    children <- possible.children[node.ages[possible.children] > node.ages[node]]
    return(children)
  }
}

#' Get edge ID
#'
#' Returns edge ID in an edge matrix from a \code{phylo} object given two nodes
#'
#' @param start First endpoint of the edge
#' @param end Second endpoint of the edge
#' @param edge Edge matrix from a \code{phylo} object
#'
#' @return Edge ID for the edge connecting \code{start} and \code{end}
#'
#' @export

get.edge.id <- function(start, end, edge){
  edge.id <- which(((edge[,1] == start) & (edge[,2] == end)) | ((edge[,1] == end) & (edge[,2] == start)))
  return(edge.id)
}

#' Plot Structured Coalescent Tree
#'
#' Plots a structured coalescent tree input as a phylo object augmented with node demes
#'
#' @param phylo phylo object augmented with node demes
#'
#' @export

structured.plot <- function(phylo, n.demes = NA){
  edge <- phylo$edge

  if (is.na(n.demes)){
    n.demes <- max(unique(phylo$node.deme))
  }
  color.palette <- rainbow(n.demes)
  edge.color <- rep(NA,dim(edge)[1])
  for (i in 1 : dim(edge)[1]){
    edge.color[i] <- color.palette[phylo$node.deme[edge[i,2]]]
  }

  plot(phylo, edge.color = edge.color, no.margin = TRUE, edge.width = 2, show.tip.label = FALSE)
}
