#' Simulation of Heterochronous Structured Coalescent
#'
#' Simulates a tree under the structured coalescent model
#'
#' @param data  nx3 matrix with first column giving the tip labels, second column the time at which the sample was taken and third column the initial deme of each sample point
#' @param effective.pop effective population size from which the sample is taken
#' @param gen.length generation length of the sampled individuals
#' @param n.deme total number of demes
#' @param migration.matrix matrix of migration rates between demes
#' @param plot.phylo logical; if FALSE (default) plot is not produced
#'
#' @return An object of class \code{phylo} (from package \code{ape}) augmented with the likelihood and log-likelihood of the simulated tree, and the deme of each node in the tree. Additional nodes are added to account for migration events between coalescences
#'
#' @export

Structured.sim <- function(data, effective.pop, gen.length, n.deme, migration.matrix, plot.phylo = FALSE){
  lambda <- effective.pop * gen.length

  if (length(lambda) == 1){
    lambda <- rep(lambda,n.deme)
  }

  n <- dim(data)[1]

  edge.list <- list()  #Listing tree edges
  edge.length <- numeric(0)  #Vector of edge lengths
  node.deme <- c(data[,3], NA)  #Initial deme of each leaf and NA for deme of root

  likelihood <- 0  #Log-likelihood of tree

  tip.age <- max(data[,2]) - data[,2]  #Ages of leaves from newest leaf

  node.height <- c(tip.age, NA)
  time <- 0

  new.node <- n+2

  active <- which(tip.age <= time)  #Active nodes (at time 0)
  height.next.tip <- suppressWarnings(min(tip.age[which(tip.age > time)]))
  count <- 1

  K <- numeric(n.deme)  #Number of lineages in each deme
  for (i in 1 : n.deme){
    K[i] <- sum(node.deme[active] == i)
  }

  diag(migration.matrix) <- 0  #Ensure no self-migration is possible (migrations i -> i prevented)

  while (height.next.tip < Inf){
    migration.rate <- sum(t(migration.matrix) %*% K)
    coalescence.rate <- sum(choose(K,2)/lambda)
    event.rate <- migration.rate + coalescence.rate  #Combined event rate of coalescence and migration

    event.prob <- pexp(height.next.tip - time, rate = event.rate)  #Probability event before next node added

    if (runif(1) > event.prob){  #New nodes added before next event
      time <- height.next.tip
      active <- c(active, which(tip.age == height.next.tip))
      height.next.tip <- suppressWarnings(min(tip.age[which(tip.age > time)]))  #Update next node time

      for (i in 1 : n.deme){
        K[i] <- sum(node.deme[active] == i)
      }
      likelihood <- likelihood + log(1-event.prob)
    }
    else{  #Event occurs before new nodes added
      event.time <- rexp.trunc(1, height.next.tip - time, event.rate)
      likelihood <- likelihood - event.rate * event.time

      if (runif(1) <= coalescence.rate/event.rate){  #Coalescence event
        coalescence.deme <- sample.int(n.deme, 1, prob = choose(K,2)/lambda)  #Select coalescence deme
        node.sample <- sample.vector(intersect(active,which(node.deme == coalescence.deme)),2)  #Select coalescing lineages
        node.deme <- c(node.deme, coalescence.deme)  #Assigning new node's deme

        time <- time + event.time
        node.height <- c(node.height, time)

        likelihood <- likelihood - log(lambda[coalescence.deme])

        #Update edge.list and edge.length
        edge.list[[count]] <- c(new.node,node.sample[1])
        edge.list[[count + 1]] <- c(new.node, node.sample[2])
        edge.length <- c(edge.length, time - node.height[node.sample])
        count <- count + 2

        K[coalescence.deme] <- K[coalescence.deme] - 1  #Number of lineages in coalescence deme decreases by 1
        active <- c(active[!(active %in% node.sample)], new.node)
        new.node <- new.node + 1
      }
      else{  #Migration event
        origin.deme <- sample.int(n.deme, 1, prob = K * rowSums(migration.matrix))
        target.deme <- sample.int(n.deme,1, prob = migration.matrix[origin.deme,])  #Select target deme
        node.sample <- sample.vector(intersect(active,which(node.deme == origin.deme)),1)  #Select lineage to migrate

        node.deme <- c(node.deme, target.deme)  #Assigning new node's deme

        time <- time + event.time
        node.height <- c(node.height, time)

        edge.list[[count]] <- c(new.node, node.sample)
        edge.length <- c(edge.length, time - node.height[node.sample])
        count <- count + 1

        likelihood <- likelihood + log(migration.matrix[origin.deme,target.deme])

        K[origin.deme] <- K[origin.deme] -1; K[target.deme] <- K[target.deme] + 1  #Number of lineages in target increases by 1, number in origin decreases by 1
        active <- c(active[!(active %in% node.sample)], new.node)
        new.node <- new.node + 1
      }
    }
  }

  while(length(active) > 1){
    migration.rate <- sum(t(migration.matrix) %*% K)
    coalescence.rate <- sum(choose(K,2)/lambda)
    event.rate <- migration.rate + coalescence.rate  #Combined event rate of coalescence and migration

    event.time <- rexp.trunc(1, height.next.tip - time, event.rate)
    likelihood <- likelihood - event.rate * event.time

    if (runif(1) <= coalescence.rate/event.rate){  #Coalescence event
      coalescence.deme <- sample.int(n.deme, 1, prob = choose(K,2)/lambda)  #Select coalescence deme
      node.sample <- sample.vector(intersect(active,which(node.deme == coalescence.deme)),2)  #Select coalescing lineages
      node.deme <- c(node.deme, coalescence.deme)  #Assigning new node's deme

      time <- time + event.time
      node.height <- c(node.height, time)

      likelihood <- likelihood - log(lambda[coalescence.deme])

      #Update edge.list and edge.length
      edge.list[[count]] <- c(new.node,node.sample[1])
      edge.list[[count + 1]] <- c(new.node, node.sample[2])
      edge.length <- c(edge.length, time - node.height[node.sample])
      count <- count + 2

      K[coalescence.deme] <- K[coalescence.deme] - 1  #Number of lineages in coalescence deme decreases by 1
      active <- c(active[!(active %in% node.sample)], new.node)
      new.node <- new.node + 1
    }
    else{  #Migration event
      origin.deme <- sample.int(n.deme, 1, prob = K * rowSums(migration.matrix))
      target.deme <- sample.int(n.deme,1, prob = migration.matrix[origin.deme,])  #Select target deme
      node.sample <- sample.vector(intersect(active,which(node.deme == origin.deme)),1)  #Select lineage to migrate

      node.deme <- c(node.deme, target.deme)  #Assigning new node's deme

      time <- time + event.time
      node.height <- c(node.height, time)

      edge.list[[count]] <- c(new.node, node.sample)
      edge.length <- c(edge.length, time - node.height[node.sample])
      count <- count + 1

      likelihood <- likelihood + log(migration.matrix[origin.deme,target.deme])

      K[origin.deme] <- K[origin.deme] -1; K[target.deme] <- K[target.deme] + 1  #Number of lineages in target increases by 1, number in origin decreases by 1
      active <- c(active[!(active %in% node.sample)], new.node)
      new.node <- new.node + 1
    }
  }

  edge <- do.call(rbind,edge.list)  #Construct edge matrix from ledge list
  edge[which(edge == max(edge))] <- n + 1  #Set root to node n+1

  Nnode <- max(edge) - n

  node.deme[n+1] <- node.deme[new.node - 1]  #Assigning deme of root at location n+1
  node.deme <- node.deme[1:(n + Nnode)]  #Removing final entry in node.deme

  Phylo.sim <- list()
  class(Phylo.sim) <- 'phylo'
  Phylo.sim$tip.label <- data[,1]
  Phylo.sim$Nnode <- Nnode
  Phylo.sim$edge <- edge
  Phylo.sim$edge.length <- edge.length
  Phylo.sim$node.deme <- node.deme
  Phylo.sim$log.likelihood <- likelihood
  Phylo.sim$likelihood <- exp(likelihood)

  if (plot.phylo == TRUE){
    color.palette <- rainbow(n.deme)
    edge.color <- rep(NA,dim(edge)[1])
    for (i in 1 : dim(edge)[1]){
      edge.color[i] <- color.palette[node.deme[edge[i,2]]]
    }

    plot(Phylo.sim, edge.color = edge.color, edge.width = 2)
    axisPhylo(1, root.time = max(data[,2]) - time, backward = FALSE)
  }
  return(Phylo.sim)
}
