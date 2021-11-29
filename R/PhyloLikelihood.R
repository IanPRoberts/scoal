#' Calculates Likelihood of a Coalescent Tree
#'
#' Returns the likelihood and log-likelihood of a (possibly heterochronous)
#' coalescent tree
#'
#' @param phylo object of class \code{phylo}
#' @param effective.pop effective population size from which the sample is taken
#' @param gen.length mean generation length of the sampled individuals
#'
#' @return A vector (\code{log-likelihood}, \code{likelihood}) giving the log-likelihood and likelihood of \code{phylo}
#'
#' @export

phylo.likelihood <- function(phylo, effective.pop, gen.length){
  lambda <- effective.pop * gen.length
  n <- length(phylo$tip.label)

  node.heights <- node.depth.edgelength(phylo)  #Horizontal distance from root
  event.times <- sort(unique(node.heights))  #Times of events, root at time=0
  time.increments <- diff(event.times)

  #Augments phylo$edge to add edge "start" and "end" times
  aug.edge <- matrix(0,nrow = (2 * (n-1)), ncol = 4)
  aug.edge[,c(1,2)] <- phylo$edge

  for (i in 1 : (2 * (n-1))){
    aug.edge[i,c(3,4)] <- c(node.heights[aug.edge[i,c(1,2)]])
  }

  #Calculate number of edges in tree midway between event times
  check.time <- event.times[-length(event.times)] + 0.5 * time.increments
  k <- rep(0,length(time.increments))  #Number of edges in tree between events
  for (i in 1 : length(time.increments)){
    k[i] <- sum( (aug.edge[,3] < check.time[i]) & (aug.edge[,4] > check.time[i]) )
  }

  likelihood <- 0  #Initialise log likelihood at 0
  likelihood <- sum(- (k * ( k-1) / (2 * lambda)) * time.increments) - (n-1) * log(lambda)
  return(c(likelihood, exp(likelihood)))  #Output (log-likelihood, likelihood)
}



#' Calculates the likelihood of a structured coalescent tree
#'
#' Returns the likelihood and log-likelihood of a tree generated under the
#' structured coalescent process
#'
#' @param phylo object of class \code{phylo} augmented with deme labels for each node
#' @param effective.pop effective population sizes from which the sample is taken
#' @param gen.length mean generation length of the sampled individuals
#' @param migration.matrix matrix of migration rates between demes
#'
#' @return A vector (\code{log-likelihood}, \code{likelihood}) giving the log-likelihood and likelihood of \code{phylo}
#'
#' @export


structured.likelihood <- function(phylo, effective.pop, gen.length, migration.matrix){
  lambda <- effective.pop * gen.length
  n <- length(phylo$tip.label) #Number of tips
  n.deme <- dim(migration.matrix)[1] #Number of demes
  n.migrations <- phylo$Nnode - n + 1 #Number of migration events
  diag(migration.matrix) <- 0  #Prevent self-migrations

  if (length(lambda) == 1){
    lambda <- rep(lambda,n.deme)
  }

  #Identify node type (migration/coalescence/leaf)
  leaf.nodes <- 1:n
  migration.nodes <- phylo$edge[!phylo$edge[,1] %in% phylo$edge[duplicated(phylo$edge[,1]),1],1] #Unique elements in col 1 of edge matrix
  coalescence.nodes <- ((n+1):(n.migrations+2*n-1))[!((n+1):(n.migrations+2*n-1) %in% migration.nodes)] #All remaining nodes

  node.heights <- node.depth.edgelength(phylo)  #Horizontal distance from root
  event.times <- sort(unique(node.heights))  #Times of events, root at time=0
  time.increments <- diff(event.times)

  #Augments phylo$edge to add edge "start" and "end" times
  n.edges <- dim(phylo$edge)[1]
  aug.edge <- matrix(0,nrow = n.edges, ncol = 4)
  aug.edge[,c(1,2)] <- phylo$edge
  aug.edge[,c(3,4)] <- node.heights[aug.edge[,c(1,2)]]

  #Number of lineages in each deme between each event time
  check.time <- event.times[-length(event.times)] + 0.5 * time.increments
  k <- matrix(0, nrow = length(check.time), ncol = n.deme)
  for (i in 1 : length(time.increments)){
    current.edges <- which((aug.edge[,3] < check.time[i]) & (aug.edge[,4] > check.time[i]) )
    for (j in 1 : n.deme){
      k[i,j] <- sum(phylo$node.deme[aug.edge[current.edges,2]] == j)
    }
  }

  #Number of coalescence events in each deme
  c <- rep(0,n.deme)
  for (j in 1 : n.deme){
    c[j] <- sum(phylo$node.deme[coalescence.nodes] == j)
  }

  #Number of migration events between each pair of demes
  m <- matrix(0,nrow = n.deme, ncol = n.deme)
  for (nodes in migration.nodes){
    i <- phylo$node.deme[aug.edge[which(aug.edge[,1] == nodes),2]]
    j <- phylo$node.deme[nodes]
    m[i,j] <- m[i,j] + 1
  }

  #Likelihood computation
  likelihood <- 0
  likelihood <- - sum(rowSums(t(t(k * (k-1)) / (2 * lambda)) + t(t(k) * rowSums(migration.matrix))) * time.increments) -
    sum(c * log(lambda)) + sum(log(migration.matrix ^ m))
  return(c(likelihood, exp(likelihood)))  #Output (log-likelihood, likelihood)
}

#' Calculates the likelihood of a structured coalescent tree
#'
#' Returns the likelihood and log-likelihood of a tree generated under the
#' structured coalescent process
#'
#' @param ED Extended data representation of a structured coalescent process
#'  @param effective.pop effective population sizes from which the samples are taken
#' @param gen.length mean generation length of the sampled individuals
#' @param migration.matrix matrix of migration rates between demes
#'
#' @return A vector (\code{log-likelihood}, \code{likelihood}) giving the log-likelihood and likelihood of \code{phylo}
#'
#' @export

ed.likelihood <- function(ED, effective.pop, gen.length, migration.matrix){
  n.deme <- dim(migration.matrix)[1]
  diag(migration.matrix) <- 0  #Prevent self-migrations

  lambda <- effective.pop * gen.length
  if (length(lambda) == 1){
    lambda <- rep(lambda,n.deme)
  }

  all.nodes <- ED[,1]
  leaf.nodes <- all.nodes[is.na(ED[,3])]
  root.node <- all.nodes[is.na(ED[,2])]
  coalescence.nodes <- all.nodes[!is.na(ED[,4])]
  coalescence.nodes <- coalescence.nodes[coalescence.nodes != root.node]
  migration.nodes <- all.nodes[! all.nodes %in% c(root.node, leaf.nodes, coalescence.nodes)]

  node.heights <- ED[,6]
  event.times <- sort(unique(node.heights))  #Times of events, root at time=0
  time.increments <- diff(event.times)

  #Number of lineages in each deme between each event time
  k <- matrix(0, nrow = length(event.times) - 1, ncol = n.deme)
  k[1,ED[root.node,5]] <- 2
  for (i in 2 : (length(event.times) - 1)){
    current.rows <- which(ED[,6] == event.times[i])
    k[i,] <- k[i-1,]
    if (length(current.rows) > 1){ #Multiple leaves added simultaneously
      for (j in current.rows){
        k[i, ED[j, 5]] <- k[i, ED[j, 5]] + 1
      }
    } else{
      if (current.rows %in% migration.nodes){ #Migration event
        k[i, ED[current.rows, 5]] <- k[i, ED[current.rows, 5]] - 1
        current.child <- ED[current.rows, 3]
        k[i, ED[current.child, 5]] <- k[i, ED[current.child, 5]] + 1
      } else if (current.rows %in% coalescence.nodes){ #Coalescence event
        k[i, ED[current.rows, 5]] <- k[i, ED[current.rows, 5]] + 1
      } else{ #Single leaf added
        k[i, ED[current.rows, 5]] <- k[i, ED[current.rows, 5]] - 1
      }
    }
  }

  nc <- ed.node.count(ED, n.deme)
  c <- nc$c
  m <- nc$m

  #Likelihood computation
  likelihood <- 0
  likelihood <- - sum(rowSums(t(t(k * (k-1)) / (2 * lambda)) + t(t(k) * rowSums(migration.matrix))) * time.increments) -
    sum(c * log(lambda)) + sum(log(migration.matrix ^ m))
  return(c(likelihood, exp(likelihood)))  #Output (log-likelihood, likelihood)
}
