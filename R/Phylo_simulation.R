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
#'
#' @export

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
  class(Phylo.sim) <- c('str_phylo', 'phylo')
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

#' Simulation of Heterochronous Structured Coalescent
#'
#' Simulates a tree under the structured coalescent model
#'
#' @param leaf_data  nx3 matrix with first column giving the tip labels, second column the time at which the sample was taken and third column the initial deme of each sample point
#' @param coal_rate effective population size from which the sample is taken
#' @param time_scale generation length of the sampled individuals
#' @param n_deme total number of demes
#' @param mig_mat matrix of migration rates between demes
#' @param plot_phylo logical; if FALSE (default) plot is not produced
#'
#' @return An object of class \code{phylo} (from package \code{ape}) augmented with the likelihood and log-likelihood of the simulated tree, and the deme of each node in the tree. Additional nodes are added to account for migration events between coalescences
#'
#' @export

scaled_sim <- function(leaf_data, coal_rate, time_scale, n_deme, mig_mat, plot_phylo = FALSE){
  n <- dim(leaf_data)[1]

  edge_list <- list()  #Listing tree edges
  edge_length <- numeric(0)  #Vector of edge lengths
  node_deme <- c(leaf_data[,3], NA)  #Initial deme of each leaf and NA for deme of root

  likelihood <- 0  #Log-likelihood of tree

  tip_age <- max(leaf_data[,2]) - leaf_data[,2]  #Ages of leaves from newest leaf

  node_height <- c(tip_age, NA)
  time <- 0

  new_node <- n+2

  active <- which(tip_age <= time)  #Active nodes (at time 0)
  height_next_tip <- suppressWarnings(min(tip_age[which(tip_age > time)]))
  count <- 1

  k <- numeric(n_deme)  #Number of lineages in each deme
  for (i in 1 : n_deme){
    k[i] <- sum(node_deme[active] == i)
  }

  diag(mig_mat) <- 0  #Ensure no self-migration is possible (migrations i -> i prevented)

  while((length(active) > 1) || (height_next_tip < Inf)){
    migration_rate <- sum(t(mig_mat) %*% k) * time_scale
    coalescence_rate <- sum(choose(k,2) * coal_rate * time_scale)
    event_rate <- migration_rate + coalescence_rate  #Combined event rate of coalescence and migration

    event_prob <- pexp(height_next_tip - time, rate = event_rate)  #Probability event before next node added

    if (runif(1) > event_prob){  #New nodes added before next event
      likelihood <- likelihood - event_rate * (height_next_tip - time)
      time <- height_next_tip
      active <- c(active, which(tip_age == height_next_tip))
      height_next_tip <- suppressWarnings(min(tip_age[which(tip_age > time)]))  #Update next node time

      for (i in 1 : n_deme){
        k[i] <- sum(node_deme[active] == i)
      }
    } else{  #Event occurs before new nodes added
      event_time <- rexp.trunc(1, height_next_tip - time, event_rate)
      likelihood <- likelihood - event_rate * event_time

      if (runif(1) <= coalescence_rate/event_rate){  #Coalescence event
        coalescence_deme <- sample.int(n_deme, 1, prob = choose(k,2)*coal_rate)
        node_sample <- sample.vector(intersect(active,which(node_deme == coalescence_deme)),2)  #Select coalescing lineages
        node_deme <- c(node_deme, coalescence_deme)  #Assigning new node's deme

        time <- time + event_time
        node_height <- c(node_height, time)

        likelihood <- likelihood + log(coal_rate[coalescence_deme] * time_scale)

        #Update edge_list and edge_length
        edge_list[[count]] <- c(new_node,node_sample[1])
        edge_list[[count + 1]] <- c(new_node, node_sample[2])
        edge_length <- c(edge_length, time - node_height[node_sample])
        count <- count + 2

        k[coalescence_deme] <- k[coalescence_deme] - 1  #Number of lineages in coalescence deme decreases by 1
        active <- c(active[!(active %in% node_sample)], new_node)
        new_node <- new_node + 1
      } else{  #Migration event
        origin_deme <- sample.int(n_deme, 1, prob = k * rowSums(mig_mat))
        target_deme <- sample.int(n_deme,1, prob = mig_mat[origin_deme,])
        node_sample <- sample.vector(intersect(active,which(node_deme == origin_deme)),1)  #Select migrating lineage

        node_deme <- c(node_deme, target_deme)  #Assigning new node's deme

        time <- time + event_time
        node_height <- c(node_height, time)

        edge_list[[count]] <- c(new_node, node_sample)
        edge_length <- c(edge_length, time - node_height[node_sample])
        count <- count + 1

        likelihood <- likelihood + log(mig_mat[origin_deme,target_deme] * time_scale)

        k[origin_deme] <- k[origin_deme] -1
        k[target_deme] <- k[target_deme] + 1  #Number of lineages in target increases by 1, number in origin decreases by 1
        active <- c(active[!(active %in% node_sample)], new_node)
        new_node <- new_node + 1
      }
    }
  }

  edge <- do.call(rbind,edge_list)  #Construct edge matrix from ledge list
  edge[which(edge == max(edge))] <- n + 1  #Set root to node n+1

  Nnode <- max(edge) - n

  node_deme[n+1] <- node_deme[new_node - 1]  #Assigning deme of root at location n+1
  node_deme <- node_deme[1:(n + Nnode)]  #Removing final entry in node_deme

  Phylo_sim <- list()
  class(Phylo_sim) <- c('str_phylo', 'phylo')
  Phylo_sim$tip.label <- leaf_data[,1]
  Phylo_sim$Nnode <- Nnode
  Phylo_sim$edge <- edge
  Phylo_sim$edge.length <- edge_length
  Phylo_sim$node.deme <- node_deme
  Phylo_sim$log.likelihood <- likelihood
  Phylo_sim$likelihood <- exp(likelihood)

  if (plot_phylo == TRUE){
    color_palette <- rainbow(n_deme)
    edge_color <- rep(NA,dim(edge)[1])
    for (i in 1 : dim(edge)[1]){
      edge_color[i] <- color_palette[node_deme[edge[i,2]]]
    }

    plot(Phylo_sim, edge.color = edge_color, edge.width = 2)
    axisPhylo(1, root.time = max(leaf_data[,2]) - time, backward = FALSE)
  }
  return(Phylo_sim)
}
