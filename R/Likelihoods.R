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

ed.likelihood <- function(ED, effective.pop, gen.length, migration.matrix, node.indices){
  n.deme <- length(effective.pop)
  lambda <- effective.pop * gen.length
  deme_decomp <- DemeDecompC(ED, n.deme, node.indices)
  node_count <- NodeCountC(ED, n.deme, node.indices)

  k <- deme_decomp$k
  time_increments <- deme_decomp$time.increments
  c <- node_count$c
  m <- node_count$m

  deme_lengths <- colSums(k * time_increments)
  coal_constants <- colSums(k * (k-1) * time_increments) / ( 2 * lambda)
  mm_row_sum <- rowSums(migration.matrix)
  log_mig_mat <- log(migration.matrix)
  diag(log_mig_mat) <- 0

  like <- sum(m * log_mig_mat) - sum(c * log(effective.pop)) - sum(coal_constants) - sum(deme_lengths * mm_row_sum)

  return(list(log.likelihood = like, likelihood = exp(like)))
}


#' Calculates the DTA likelihood for a structured coalescent tree
#'
#' Returns the DTA likelihood and log-DTA-likelihood of a structured coalescent tree
#'
#' @param ED Extended data representation of a structured coalescent process
#'  @param effective_population effective population sizes from which the samples are taken
#' @param gen_length mean generation length of the sampled individuals
#' @param migration_matrix matrix of migration rates between demes
#' @param node_indices Vector giving row indices for each node label
#'
#' @return A vector (\code{log-likelihood}, \code{likelihood}) giving the log-likelihood and likelihood of \code{phylo}
#'
#' @export

dta.likelihood <- function(ED, effective_population, gen_length, migration_matrix, node_indices){
  root_row <- which(is.na(ED[,2]))

  f_mm <- forward.migration.matrix(migration_matrix, effective_population)
  f_mm_rowsums <- rowSums(f_mm)
  log_likelihood <- 0

  for (i in (1 : dim(ED)[1])[-root_row]){
    node_parent_row <- node_indices[ED[i,2]]
    time_increment <- ED[i, 6] - ED[node_parent_row, 6]

    parent_deme <- ED[i, 5]
    if (is.na(ED[i,3])){ #Leaf node
      node_deme <- parent_deme
    } else{
      node_deme <- ED[node_indices[ED[i,3]], 5]
    }

    log_likelihood <- log_likelihood - f_mm_rowsums[parent_deme] * time_increment

    if (parent_deme != node_deme){
      log_likelihood <- log_likelihood + log(f_mm[parent_deme, node_deme])
    }
  }
  return(list(log.likelihood = log_likelihood, likelihood = exp(log_likelihood)))
}

#' Calculates synthetic likelihood for a coalescent tree with migration history
#'
#' Computes a synthetic likelihood for a coalescent tree with migration history
#' with a Poisson marginal distribution on the number of migrations
#'
#' @param ED Extended data representation of a structured coalescent process
#' @param rate Poisson rate for synthetic likelihood
#' @param effective.pop Effective population sizes
#' @param migration.matrix Migration matrix
#' @param node.indices
#'
#' @return A vector (\code{log-likelihood}, \code{likelihood}) giving the log-likelihood and likelihood of \code{phylo}
#'
#' @export

synth.likelihood <- function(ED, effective.pop, rate, migration.matrix, node.indices){
  root.row <- which(is.na(ED[,2]))
  n.deme <- length(effective.pop)
  tree.length <- 0

  for (i in (1 : dim(ED)[1])[-root.row]){
    node.parent <- ED[i,2]
    parent.row <- node.indices[node.parent]
    tree.length <- tree.length + ED[i,6] - ED[parent.row,6]
  }

  M <- dim(ED)[1] - root.row
  like <- M * (log(rate) - log(tree.length) - log(n.deme - 1)) - log(n.deme) - rate

  return(list(log.likelihood = like, likelihood = exp(like)))
}
