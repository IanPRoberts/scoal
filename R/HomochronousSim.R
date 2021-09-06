#' Simulation of Homochronous Coalescent Trees
#'
#' Simulates a coalescent tree by estimating coalescence times backwards in time
#' for a homochronous sample taken all at one time.
#'
#' @param Data vector of labels to be assigned to the tips
#' @param EffectivePop effective population size from which the sample is taken
#' @param GenLength generation length of the sampled individuals
#'
#' @return An object of class \code{phylo} (from package \code{ape})
#'
#' @export

Homochronous.Sim <- function(Data,EffectivePop,GenLength){
  lambda <- EffectivePop * GenLength
  n <- length(Data)

  Nnode <- n-1  #Number of internal nodes
  edge <- matrix(NA,n+Nnode -1,2)  #Edge matrix
  edge.length <- numeric(n+Nnode - 1)  #Edge length vector

  node.age <- numeric(n+Nnode)  #Age of node (total distance from time 0)
  time <- 0

  active <- 1:n  #Active nodes (at time 0)
  new.node <- 2*n-1  #New node to connect to

  for (i in 1:(n-1)){
    node.sample <- sample(seq_along(active),2)  #Indices of nodes to merge
    coal.time <- rexp(1,choose(n+1-i,2)/lambda)  #Time until next coalescence
    time <- time + coal.time
    node.age[new.node] <- time

    rows <- c(2*i - 1, 2*i)  #Rows of edge matrix to change

    #Update edge matrix and edge.length
    edge[rows,1] <- new.node
    edge[rows,2] <- active[node.sample]
    edge.length[rows] <- time - node.age[active[node.sample]]

    #Update active nodes
    active <- c(active[-node.sample],new.node)

    new.node <- new.node - 1 #new.node for next coalescence
  }

  Phylo.sim <- list()
  class(Phylo.sim) <- 'phylo'
  Phylo.sim$tip.label <- Data
  Phylo.sim$Nnode <- Nnode
  Phylo.sim$edge <- edge
  Phylo.sim$edge.length <- edge.length

  Phylo.sim
}
