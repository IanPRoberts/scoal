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
  class(Phylo_sim) <- 'phylo'
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
