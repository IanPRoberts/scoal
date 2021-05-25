#' Calculates Likelihood of a Coalescent Tree
#'
#' Returns the likelihood and log-likelihood of a (possibly heterochronous)
#' coalescent tree
#'
#' @param phylo object of class \code{phylo}
#' @param effective.pop effective population size from which the sample is taken
#' @param gen.length generation length of the sampled individuals
#'
#' @return A vector (\code{log-likelihood}, \code{likelihood}) giving the log-likelihood and likelihood of \code{phylo}
phylo.likelihood <- function(phylo, effective.pop, gen.length){
  lambda <- effective.pop * gen.length
  n <- length(phylo$tip.label)

  node.heights <- node.depth.edgelength(phylo)  #Horizontal distance from root
  event.times <- sort(unique(node.heights))  #Times of events, root at time=0
  time.increments <- event.times[2 : length(event.times)] - event.times[1 : (length(event.times) - 1)]

  #Augments phylo$edge to add edge "start" and "end" times
  aug.edge <- matrix(0,nrow = (2 * (n-1)), ncol = 4)
  aug.edge[,c(1,2)] <- phylo$edge

  for (i in 1 : (2 * (n-1))){  ##### REPLACE LOOP WITH APPLY? #####
    aug.edge[i,c(3,4)] <- c(node.heights[aug.edge[i,c(1,2)]])
  }

  #Check how many edges in tree midway between event times
  check.time <- event.times[-length(event.times)] + 0.5 * time.increments

  k <- rep(0,length(time.increments))  #Number of edges in tree between events
  for (i in 1 : length(time.increments)){
    k[i] <- sum( (aug.edge[,3] < check.time[i]) & (aug.edge[,4] > check.time[i]) )
  }

  likelihood <- 0  #Initialise log likelihood at 0
  likelihood <- sum(- (k * ( k-1) / (2 * lambda)) * time.increments) - (n-1) * log(lambda)
  c(likelihood, exp(likelihood))  #Output (log-likelihood, likelihood)
}
