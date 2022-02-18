#' Migration Pair Birth MCMC Move
#'
#' Performs a Migration pair birth move (Ewing et al. 2004). Adds two migration
#' nodes on an edge selected uniformly at random from a structured coalescent
#' process, allocating a deme for the added edge such that a migration event
#' does not target its origin deme.
#'
#' @param ED Extended data object; matrix with columns Node ID, parent, child 1, child 2, deme, node age
#' @param n.deme Number of distinct demes in the population
#'
#' @return Updated extended data object with the proposal from the migration pair birth move
#'
#' @export

ed.mig.pair.birth <- function(ED, n.deme, node.indices){
  root.node <- which(is.na(ED[,2]))

  #Sample non-root node to obtain edge <sampled.node, node.parent>
  selected.node <- sample(ED[-root.node, 1], 1)
  selected.row <- node.indices[selected.node]
  parent.node <- ED[selected.row, 2]
  parent.row <- node.indices[parent.node]

  new.nodes <- max(ED[,1]) + c(1,2)  #New node IDs
  new.node.ages <- sort(runif(2, min = ED[parent.row, 6], max = ED[selected.row, 6]), decreasing = TRUE) #New node ages

  new.deme <- sample.vector((1:n.deme)[-ED[selected.row,5]], 1)

  ED[selected.row, 2] <- new.nodes[1]
  ED[parent.row, 2 + which(ED[parent.row,3:4] == selected.node)] <- new.nodes[2]
  ED <- rbind(ED,
              c(new.nodes[1], new.nodes[2], selected.node, NA, new.deme, new.node.ages[1]),
              c(new.nodes[2], parent.node, new.nodes[1], NA, ED[selected.row, 5], new.node.ages[2])
              )

  n.nodes <- dim(ED)[1]
  prop.ratio <- (n.deme - 1) * (n.nodes - 1) * (ED[selected.row, 6] - ED[parent.row, 6])^2 / (2 * (n.nodes + 1))
  log.prop.ratio <- log(n.deme - 1) + log(n.nodes - 1) + 2 * log(abs(ED[selected.row, 6] - ED[parent.row, 6])) - log(2) - log(n.nodes + 1)

  if (new.nodes[2] > length(node.indices)){
    new.node.indices <- numeric(new.nodes[2])
    new.node.indices[1:length(node.indices)] <- node.indices
  } else{
    new.node.indices <- node.indices
  }
  new.node.indices[new.nodes] <- c(dim(ED)[1] - 1, dim(ED)[1])
  return(list(ED = ED, prop.ratio = prop.ratio, log.prop.ratio = log.prop.ratio, node.indices = new.node.indices))
}

#' Migration Pair Death MCMC Move
#'
#' Performs a Migration pair death move (Ewing et al. 2004). Deletes two
#' migration nodes if they lie between two edges in the same deme.
#'
#' @param ED Extended data object; matrix with columns Node ID, parent, child 1, child 2, deme, node age
#' @param n.deme Number of distinct demes in the population
#'
#' @return Updated extended data object with the proposal from the migration pair death move
#'
#' @export

ed.mig.pair.death <- function(ED, n.deme, node.indices){
  leaf.nodes <- ED[is.na(ED[,3]), 1]
  root.node <- ED[is.na(ED[,2]), 1]
  coalescence.nodes <- ED[!is.na(ED[,4]), 1]
  coalescence.nodes <- coalescence.nodes[coalescence.nodes != root.node]
  migration.nodes <- ED[ is.na(ED[,4]) & (!is.na(ED[,3])) ,1]

  #Sample non-root node to obtain edge <sampled.node, node.parent>
  selected.node <- sample(ED[-root.node, 1], 1)
  selected.row <- node.indices[selected.node]
  parent.node <- ED[selected.row, 2]
  parent.row <- node.indices[parent.node]

  if ((ED[selected.row, 1] %in% migration.nodes) &&
      (ED[parent.row, 1] %in% migration.nodes) &&
      (ED[node.indices[ED[selected.row, 3]], 5] == ED[parent.row, 5])){
    #Both ends of the edge are migration nodes, and the demes are consistent to remove the pair of nodes
    parent2.node <- ED[parent.row, 2]
    parent2.row <- node.indices[parent2.node]
    child.node <- ED[selected.row, 3]
    child.row <- node.indices[child.node]

    prop.ratio <- 2 * (n.nodes - 1) / ((n.nodes - 3) * (n.deme - 1) * (ED[child.row, 6] - ED[parent2.row, 6])^2)
    log.prop.ratio <- log(2) + log(n.nodes - 1) - log(n.nodes - 3) - log(n.deme - 1) - 2 * log(abs(ED[child.row, 6] - ED[parent2.row, 6]))

    ED[parent2.row, 2 + which(ED[parent2.row, 3:4] == parent.node)] <- child.node
    ED[child.row, 2] <- parent2.node
    ED <- ED[-c(selected.row, parent.row),]

    node.indices[c(selected.node, parent.node)] <- 0
    index.changes.1 <- (node.indices > selected.row)
    index.changes.2 <- (node.indices > parent.row)
    node.indices[index.changes.1] <- node.indices[index.changes.1] - 1
    node.indices[index.changes.2] <- node.indices[index.changes.2] - 1

    return(list(ED = ED, prop.ratio = prop.ratio, log.prop.ratio = log.prop.ratio, node.indices = node.indices))
  } else {
    return(list(ED = ED, prop.ratio = 0, log.prop.ratio = -Inf, node.indices = node.indices))
  }
}
