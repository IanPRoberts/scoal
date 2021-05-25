#' Simulation of Heterochronous Coalescent Trees
#'
#' Simulates a coalescent tree by estimating coalescence times backwards in time
#' for a heterochronous sample taken at different time.
#'
#' @param Data nx2 matrix with first column giving the tip labels and second column the time at which sample was taken
#' @param EffectivePop effective population size from which the sample is taken
#' @param GenLength generation length of the sampled individuals
#' @param phylo.plot logical; if FALSE (default) plot is not produced
#'
#' @return An object of class \code{phylo} (from package \code{ape}) augmented with the likelihood and log-likelihood of the simulated tree

Heterochronous.Sim <- function(Data, effective.pop, gen.length, phylo.plot=FALSE){
  lambda <- effective.pop * gen.length

  n <- dim(Data)[1]

  Nnode <- n-1  #Number of internal nodes
  edge <- matrix(NA,n+Nnode -1,2)  #Edge matrix
  edge.length <- numeric(n+Nnode - 1)  #Edge length vector
  likelihood <- 0  #Likelihood of tree simulation, calculated as log likelihood and exponentiated at end

  tip.age <- max(Data[,2]) - Data[,2]  #Ages of tips from most recent tip

  node.height <- rep(Inf,n+Nnode)  #Height of node from furthest tips
  node.height[1:n] <- tip.age
  time <- 0

  active <- which(tip.age <= time)  #Active sample nodes (at time 0)
  new.node <- 2*n-1  #New node to connect to
  height.next.tip <- suppressWarnings(min(tip.age[which(tip.age > time)]))

  count <- 1

  #Tree until all observed lineages have been introduced
  while (height.next.tip < Inf){
    k <- length(active)  #Current Number Lineages
    exp.rate <- k * (k-1) / (2 * lambda)
    prob.coalesce <- pexp(height.next.tip - time, rate = exp.rate)  #Probability coalescence before next node added

    if (runif(1) > prob.coalesce){  #New nodes added before coalescence
      time <- height.next.tip
      active <- c(active, which(tip.age == height.next.tip))
      height.next.tip <- suppressWarnings(min(tip.age[which(tip.age > time)]))  #Update next node time
      likelihood <- likelihood + log(1 - prob.coalesce)
    }
    else {  #Coalescence before new nodes added
      node.sample <- sample(seq_along(active),2)  #Indices of nodes to merge

      #Time until next coalescence is (upper) truncated exponential with cutoff
      #height.next.tip - time; see LaTeX Workings for details

      coal.time <- rexp.trunc(1, height.next.tip - time, exp.rate)
      likelihood <- likelihood + log(dexp(coal.time,exp.rate)) - log(k * (k-1) / 2)

      time <- time + coal.time
      node.height[new.node] <- time

      #Update edge matrix and edge.length
      rows <- c(2*count - 1, 2*count)  #Rows of edge matrix to change
      edge[rows,1] <- new.node
      edge[rows,2] <- active[node.sample]
      edge.length[rows] <- time - node.height[active[node.sample]]

      active <- c(active[-node.sample],new.node)  #Update active nodes
      new.node <- new.node - 1 #new.node for next coalescence
      count <- count + 1
    }
  }

  #Tree after all lineages have been introduced
  k <- length(active)
  for (i in 1 : (k-1)){
    node.sample <- sample(seq_along(active),2)  #Indices of nodes to merge
    exp.rate <- (k+1-i) * (k-i) /(2 * lambda)
    coal.time <- rexp(1,exp.rate)  #Time until next coalescence
    likelihood <- likelihood + log(dexp(coal.time, exp.rate)) - log((k+1-i) * (k-i) /2)

    time <- time + coal.time
    node.height[new.node] <- time

    #Update edge matrix and edge.length
    rows <- c(2*(i + count) - 3, 2*(i+count)-2)  #Rows of edge matrix to change
    edge[rows,1] <- new.node
    edge[rows,2] <- active[node.sample]
    edge.length[rows] <- time - node.height[active[node.sample]]

    #Update active nodes
    active <- c(active[-node.sample],new.node)

    new.node <- new.node - 1 #new.node for next coalescence
  }

  Phylo.sim <- list()
  class(Phylo.sim) <- 'phylo'
  Phylo.sim$tip.label <- Data[,1]
  Phylo.sim$Nnode <- Nnode
  Phylo.sim$edge <- edge
  Phylo.sim$edge.length <- edge.length
  Phylo.sim$log.likelihood <- likelihood
  Phylo.sim$likelihood <- exp(likelihood)

  if (phylo.plot == TRUE){
    plot(Phylo.sim, show.tip.label = FALSE)
    axisPhylo(1, root.time = max(Data[,2]) - time, backward = FALSE)
  }
  Phylo.sim
}
