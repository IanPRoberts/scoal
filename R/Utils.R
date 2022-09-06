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

structured.plot <- function(x, n_deme = NA){
  if (! ('phylo' %in% class(x))){
    phylo <- ed.to.phylo(x)
  } else{
    phylo <- x
  }

  edge <- phylo$edge

  if (is.na(n_deme)){
    n_deme <- max(unique(phylo$node.deme))
  }

  if (any(phylo$node.deme == 0)){
    color.palette <- c(rainbow(n_deme), "black")
    phylo$node.deme[phylo$node.deme == 0] <- n_deme + 1
  } else{
    color.palette <- rainbow(n_deme)
  }

  edge.color <- color.palette[phylo$node.deme[edge[,2]]]

  plot(phylo, edge.color = edge.color, no.margin = TRUE, edge.width = 2, show.tip.label = FALSE)
  par(mar = 0.1 + c(5,4,4,2))
}

#' Normal Distribution with Reflecting Boundary
#'
#' Generates normally distributed variates subject to reflecting boundary
#' conditions
#'
#' @param n Number of variates to generate
#' @param mean Mean of the normal distribution
#' @param sd Standard deviation of the normal distribution
#' @param lower Lower reflecting boundary (possibly infinite)
#' @param upper Upper reflecting boundary (possibly infinite)
#'
#' @export

rnorm.reflect <- function(n, mean = 0, sd = 1, lower = -Inf, upper = Inf){
  sample <- rnorm(n, mean, sd)

  if (upper < lower){
    stop('upper < lower')
  }

  for (i in 1:n){
    while ((sample[i] < lower) | (sample[i] > upper)){
      if (sample[i] < lower){
        sample[i] <- lower + (lower - sample[i])
      } else{
        sample[i] <- upper - (sample[i] - upper)
      }
    }
  }

  return(sample)
}

