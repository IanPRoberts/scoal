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

#' The Truncated Exponential Distribution
#'
#' Density, distribution function and random generation for the exponential
#' distribution with rate \code{rate} truncated above by \code{c}
#'
#' @param x,q vector of quantiles
#' @param n number of observations
#' @param rate vector of rates
#' @param c truncation cutoff
#'
#' @export

dexp.trunc <- function(x, c = Inf, rate = 1){
  #Density function
  dexp <- rate * exp(-rate * x) / (1 - exp(-rate * c))
  dexp
}

#' @rdname dexp.trunc

pexp.trunc <- function(q, c = Inf, rate = 1){
  #Distribution function
  pexp <- (1 - exp(-rate * q)) / (1 - exp(-rate*c))
  pexp
}

#' @rdname dexp.trunc

rexp.trunc <- function(n, c = Inf, rate = 1){
  #i.i.d. sample of size n
  rexp <- runif(n)
  rexp <- - (1/rate) * log(1 - (1 - exp(-rate * c)) * rexp)
  rexp
}
