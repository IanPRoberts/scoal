#' Simulation of Heterochronous Structured Coalescent
#'
#' Simulates a tree under the structured coalescent model
#'
#' @param data $n \times 3$ matrix with first column giving the tip labels, second column the time at which the sample was taken and third column the initial deme of each sample point
#' @param effective.pop effective population size from which the sample is taken
#' @param gen.length generation length of the sampled individuals
#' @param n.deme total number of demes
#' @param migration.matrix matrix of migration rates between demes
#' @param phylo.plot logical; if FALSE (default) plot is not produced
#'
#' @return An object of class \code{phylo} (from package \code{ape}) augmented with the likelihood and log-likelihood of the simulated tree, and the deme of each node in the tree. Additional nodes are added to account for migration events between coalescences

Structured.sim <- function(data, effective.pop, gen.length, n.deme, migration.matrix, phylo.plot = FALSE){
  lambda <- effective.pop * gen.length

  n <- dim(data)[1]

  edge.list <- list()  #Listing edges in the tree
  edge.length <- numeric(0)  #Edge length vector
  node.deme <- c(data[,3], NA)  #Initial deme of each node and NA for deme of root

  likelihood <- 0  #Likelihood of tree simulation, calculated as log likelihood and exponentiated at end

  tip.age <- max(data[,2]) - data[,2]  #Ages of tips from most recently sampled tip

  node.height <- c(tip.age, NA)
  time <- 0

  new.node <- n+2

  active <- which(tip.age <= time)  #Active sample nodes (at time 0)
  height.next.tip <- suppressWarnings(min(tip.age[which(tip.age > time)]))
  count <- 1

  K <- numeric(n.deme)  #Number of lineages in each deme
  for (i in 1 : n.deme){
    K[i] <- sum(node.deme[active] == i)
  }

  while (height.next.tip < Inf){
    migration.rate <- sum(t(migration.matrix) %*% K)
    coalescence.rate <- sum(choose(K,2)/effective.pop)/gen.length
    event.rate <- migration.rate + coalescence.rate  #Combined event rate of coalescences and migrations

    event.prob <- pexp(height.next.tip - time, rate = event.rate)  #Probability coalescence/migration event occurs before next node added

    if (runif(1) > event.prob){  #New nodes added before next coalescence/migration
      time <- height.next.tip
      active <- c(active, which(tip.age == height.next.tip))
      height.next.tip <- suppressWarnings(min(tip.age[which(tip.age > time)]))  #Update next node time

      for (i in 1 : n.deme){
        K[i] <- sum(node.deme[active] == i)
      }
    }
    else{  #Coalescence/migration event occurs before new nodes added
      event.time <- rexp.trunc(1, height.next.tip - time, event.rate)  #Time of next event

      if (runif(1) <= coalescence.rate/event.rate){  #Next event is a coalescence
        possible.coal.demes <- which(K >= 2)  #Possible demes in which a coalescence may occur
        coalescence.deme <- sample.vector(possible.coal.demes,1, prob = K[possible.coal.demes])  #Selected deme for coalescence
        node.sample <- sample.vector(intersect(active,which(node.deme == coalescence.deme)),2)  #Nodes to coalesce

        node.deme <- c(node.deme, coalescence.deme)  #Deme of new.node

        time <- time + event.time
        node.height <- c(node.height, time)

        #Update edge.list and edge.length
        edge.list[[count]] <- c(new.node,node.sample[1])
        edge.list[[count + 1]] <- c(new.node, node.sample[2])
        edge.length <- c(edge.length, time - node.height[node.sample])
        count <- count + 2

        K[coalescence.deme] <- K[coalescence.deme] - 1  #Number of lineages in coalescence deme decreases by 1
        active <- c(active[!(active %in% node.sample)], new.node)
        new.node <- new.node + 1
      }
      else{  #Next event is a migration event
        possible.origin.demes <- which(K >= 1)  #Possible demes where a migration may originate
        origin.deme <- sample.vector(possible.origin.demes,1,prob = K[possible.origin.demes])  #Deme for migration origin

        target.deme <- sample.int(n.deme,1, prob = migration.matrix[origin.deme,])  #Deme for migration target
        node.sample <- sample.vector(intersect(active,which(node.deme == origin.deme)),1)  #Node to migrate from origin to target

        node.deme <- c(node.deme, target.deme)  #Deme of new.node

        time <- time + event.time
        node.height <- c(node.height, time)

        edge.list[[count]] <- c(new.node, node.sample)
        edge.length <- c(edge.length, time - node.height[node.sample])
        count <- count + 1

        K[origin.deme] <- K[origin.deme] -1; K[target.deme] <- K[target.deme] + 1  #Number of lineages in target increases by 1, number in origin decreases by 1
        active <- c(active[!(active %in% node.sample)], new.node)
        new.node <- new.node + 1
      }
    }
  }

  while(length(active) > 1){
    migration.rate <- sum(t(migration.matrix) %*% K)
    coalescence.rate <- sum(choose(K,2)/effective.pop)/gen.length
    event.rate <- migration.rate + coalescence.rate  #Combined event rate of coalescences and migrations

    event.time <- rexp.trunc(1, height.next.tip - time, event.rate)  #Time of next event

    if (runif(1) <= coalescence.rate/event.rate){  #Next event is a coalescence
      possible.coal.demes <- which(K >= 2)  #Possible demes in which a coalescence may occur
      coalescence.deme <- sample.vector(possible.coal.demes,1, prob = K[possible.coal.demes])  #Selected deme for coalescence
      node.sample <- sample.vector(intersect(active,which(node.deme == coalescence.deme)),2)  #Nodes to coalesce
      node.deme <- c(node.deme, coalescence.deme)  #Deme of new.node is same as coalescence.deme

      time <- time + event.time
      node.height <- c(node.height, time)

      #Update edge.list and edge.length
      edge.list[[count]] <- c(new.node,node.sample[1])
      edge.list[[count + 1]] <- c(new.node, node.sample[2])
      edge.length <- c(edge.length, time - node.height[node.sample])
      count <- count + 2

      K[coalescence.deme] <- K[coalescence.deme] - 1  #Number of lineages in coalescence deme decreases by 1
      active <- c(active[!(active %in% node.sample)], new.node)
      new.node <- new.node + 1
    }
    else{  #Next event is a migration event
      possible.origin.demes <- which(K >= 1)  #Possible demes where a migration may originate
      origin.deme <- sample.vector(possible.origin.demes,1,prob = K[possible.origin.demes])  #Deme for migration origin
      target.deme <- sample.int(n.deme,1, prob = migration.matrix[origin.deme,])  #Deme for migration target

      node.sample <- sample.vector(intersect(active,which(node.deme == origin.deme)),1)  #Node to migrate from origin to target
      node.deme <- c(node.deme, target.deme)  #Deme of new.node is same as target.deme

      time <- time + event.time
      node.height <- c(node.height, time)

      edge.list[[count]] <- c(new.node, node.sample)
      edge.length <- c(edge.length, time - node.height[node.sample])
      count <- count + 1

      K[origin.deme] <- K[origin.deme] -1; K[target.deme] <- K[target.deme] + 1  #Number of lineages in target increases by 1, number in origin decreases by 1
      active <- c(active[!(active %in% node.sample)], new.node)
      new.node <- new.node + 1
    }
  }

  edge <- do.call(rbind,edge.list)  #Construct edge matrix from list of edges
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
  #Phylo.sim$log.likelihood <- likelihood
  #Phylo.sim$likelihood <- exp(likelihood)

  color.palette <- rainbow(n.deme)
  edge.color <- rep(NA,dim(edge)[1])
  for (i in 1 : dim(edge)[1]){
    edge.color[i] <- color.palette[node.deme[edge[i,2]]]
  }

  plot(Phylo.sim, edge.color = edge.color)
  axisPhylo(1, root.time = max(data[,2]) - time, backward = FALSE)

  Phylo.sim
}
