#' Migration Birth MCMC Move
#'
#' Performs a Migration birth move (Ewing et al. 2004). Adds a migration
#' node between a randomly selected ancestral node and its parent, allocating a
#' deme for the new edge consistent with the surrounding edges.
#'
#'
#' @param ED Extended data object; matrix with columns Node ID, parent, child 1, child 2, deme, node age
#' @param n.deme Number of distinct demes in the population
#'
#' @return Updated extended data object with the proposal from the migration birth move
#'
#' @export

ed.mig.birth <- function(ED, n.deme, fix.leaf.deme = TRUE, node.indices){
  root.node <- ED[is.na(ED[,2]), 1]
  coalescence.nodes <- ED[!is.na(ED[,4]),1]
  coalescence.nodes <- coalescence.nodes[coalescence.nodes != root.node]
  leaf.nodes <- ED[(is.na(ED[,3])) & (is.na(ED[,4])), 1]
  migration.nodes <- ED[ is.na(ED[,4]) & (!is.na(ED[,3])) ,1]

  parent.rows <- node.indices[ED[,2]]
  parent.times <- ED[parent.rows, 6]
  parent.times[root.node] <- 0
  edge.length <- ED[,6] - parent.times

  tree.length <- sum(edge.length)  #Total tree length
  new.location <- runif(1, 0, tree.length)

  child.row <- min(which(cumsum(edge.length) >= new.location))
  child.node <- ED[child.row, 1]

  parent.node <- ED[child.row, 2]
  parent.row <- node.indices[parent.node]

  #Calculating subtree containing edge <new.node, child.node> terminating in migration or leaf nodes
  subtree.nodes <- child.node
  active.nodes <- child.node
  while (any(active.nodes %in% coalescence.nodes)){
    active.nodes <- active.nodes[active.nodes %in% coalescence.nodes]

    active.rows <- node.indices[active.nodes]
    subtree.nodes <- c(subtree.nodes, ED[active.rows, 3:4]) #Add children of active.nodes to subtree
    active.nodes <- ED[active.rows, 3:4]
  }

  if ((fix.leaf.deme == TRUE) && (any(subtree.nodes %in% leaf.nodes))){
    # REJECT
    return(list(ED = ED, prop.ratio = 0, log.prop.ratio = -Inf, node.indices = node.indices))
  } else{
    #Continue proposal
    subtree.leaves <- subtree.nodes[subtree.nodes %in% migration.nodes]
    old.deme <- ED[child.row, 5]

    proposal.deme <- sample.vector((1:n.deme)[- old.deme], 1)  #Propose deme update without accounting for surrounding demes

    subtree.leaf.rows <- node.indices[subtree.leaves]
    subtree.leaf.children <- ED[subtree.leaf.rows, 3]
    subtree.leaf.children.rows <- node.indices[subtree.leaf.children]

    if (any(ED[subtree.leaf.children.rows, 5] == proposal.deme)){
      #REJECT
      return(list(ED = ED, prop.ratio = 0, log.prop.ratio = -Inf, node.indices = node.indices))
    }

    #Update deme across subtree
    subtree.rows <- node.indices[subtree.nodes]
    ED[subtree.rows, 5] <- proposal.deme

    new.node <- max(ED[,1]) + 1

    #Update child.row with new parent
    ED[child.row, 2] <- new.node

    #Update parent row with new child
    which.child <- which(ED[parent.row, 3:4] == child.node)
    ED[parent.row, 2 + which.child] <- new.node  #which.child either 1 or 2 so add 2 to get to col no.

    #Add row for new.node
    cumulative.edge.length <- c(0, cumsum(edge.length)[-length(edge.length)])  #Offset cumsum(edge.length) by 1 in case new.location is on first edge in edge length
    new.age <- ED[child.row, 6] - (new.location - max(cumulative.edge.length[cumulative.edge.length < new.location])) #Age of new.node (currently have distance along tree from new.location)

    ED <- rbind(ED, c(new.node, parent.node, child.node, NA, old.deme, new.age))
    prop.ratio <- (n.deme - 1) * tree.length / (length(migration.nodes) + 1)
    log.prop.ratio <- log(n.deme - 1) + log(tree.length) - log(length(migration.nodes) + 1)

    if (new.node > length(node.indices)){
      new.node.indices <- numeric(new.node)
      new.node.indices[1:length(node.indices)] <- node.indices
    } else{
      new.node.indices <- node.indices
    }
    new.node.indices[new.node] <- dim(ED)[1]

    return(list(ED= ED, prop.ratio = prop.ratio, log.prop.ratio = log.prop.ratio, node.indices = new.node.indices))
  }
}


#' Migration Death MCMC Move
#'
#' Performs a Migration death move (Ewing et al. 2004). Selects an ancestral node
#' uniformly at random and removes the parent of the selected node, allocating
#' demes as necessary for modified edges.
#'
#'
#' @param ED Extended data object; matrix with columns Node ID, parent, child 1, child 2, deme, node age
#' @param n.deme Number of distinct demes in the population
#'
#' @return Updated extended data object with the proposal from the migration death move
#'
#' @export

ed.mig.death <- function(ED, n.deme, fix.leaf.deme = TRUE, node.indices){
  root.node <- ED[is.na(ED[,2]), 1]
  coalescence.nodes <- ED[!is.na(ED[,4]),1]
  coalescence.nodes <- coalescence.nodes[coalescence.nodes != root.node]
  leaf.nodes <- ED[(is.na(ED[,3])) & (is.na(ED[,4])), 1]
  migration.nodes <- ED[ is.na(ED[,4]) & (!is.na(ED[,3])) ,1]

  parent.rows <- node.indices[ED[,2]]
  parent.times <- ED[parent.rows, 6]
  parent.times[root.node] <- 0
  edge.length <- ED[,6] - parent.times

  tree.length <- sum(edge.length)  #Total tree length

  if (length(migration.nodes) == 0){
    #REJECT
    return(list(ED = ED, prop.ratio = 0, log.prop.ratio = -Inf, node.indices = node.indices))
  }
  selected.node <- sample.vector(migration.nodes, 1)
  selected.row <- node.indices[selected.node]

  child.node <- ED[selected.row, 3]
  child.row <- node.indices[child.node]

  #Calculating subtree containing edge <new.node, child.node> terminating in migration or leaf nodes
  subtree.nodes <- child.node
  active.nodes <- child.node
  while (any(active.nodes %in% coalescence.nodes)){
    active.nodes <- active.nodes[active.nodes %in% coalescence.nodes]
    active.rows <- node.indices[active.nodes]
    subtree.nodes <- c(subtree.nodes, ED[active.rows, 3:4]) #Add children of active.nodes to subtree
    active.nodes <- ED[active.rows, 3:4]
  }

  if ((fix.leaf.deme == TRUE) && (any(subtree.nodes %in% leaf.nodes))){
    # REJECT
    return(list(ED = ED, prop.ratio = 0, log.prop.ratio = -Inf, node.indices = node.indices))
  } else{
    #Continue proposal
    subtree.leaves <- subtree.nodes[subtree.nodes %in% migration.nodes]

    subtree.leaf.rows <- node.indices[subtree.leaves]
    subtree.leaf.children <- ED[subtree.leaf.rows, 3]
    subtree.leaf.children.rows <- node.indices[subtree.leaf.children]
    forbidden.demes <- unique(ED[subtree.leaf.children.rows, 5])

    proposal.deme <- ED[selected.row, 5]  #Propose deme immediately above subtree
    if (proposal.deme %in% forbidden.demes){
      #REJECT - cannot change deme
      return(list(ED = ED, prop.ratio = 0, log.prop.ratio = -Inf, node.indices = node.indices))
    }

    parent.node <- ED[selected.row, 2]
    parent.row <- node.indices[parent.node]

    #Update deme across subtree
    subtree.rows <- node.indices[subtree.nodes]
    ED[subtree.rows, 5] <- proposal.deme

    ED[child.row, 2] <- parent.node  #Update child.node parent

    #Update parent.node child
    which.child <- which(ED[parent.row, 3:4] == selected.node)
    ED[parent.row, 2 + which.child] <- child.node

    node.indices[selected.node] <- 0
    index.changes <- (node.indices > selected.row)
    node.indices[index.changes] <- node.indices[index.changes] - 1

    ED <- ED[-selected.row,]
    prop.ratio <- length(migration.nodes)/ ((n.deme - 1) * tree.length)
    log.prop.ratio <- log(length(migration.nodes)) - log(n.deme - 1) - log(tree.length)
    return(list(ED = ED, prop.ratio = prop.ratio, log.prop.ratio = log.prop.ratio, node.indices = node.indices))
  }
}
