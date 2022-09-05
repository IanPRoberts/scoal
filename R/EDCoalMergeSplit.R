#' Coalescent Node Merge Proposal
#'
#' Performs a coalescent node merge move (Ewing et al. 2004). Merges two migration
#' nodes immediately below a coalescent node and place above coalescent node
#'
#' @param ED Extended data object; matrix with columns Node ID, parent, child 1, child 2, deme, node age
#' @param n.deme Number of distinct demes in the population
#' @param node.indices Vector giving row indices for node labels
#'
#' @return Updated extended data object with the proposal from the migration pair birth move
#'
#' @export

ed.coal.merge <- function(ED, n.deme, node.indices){
  root.node <- ED[is.na(ED[,2]), 1]
  coalescence.nodes <- ED[!is.na(ED[,4]),1]

  selected.node <- sample.vector(coalescence.nodes, 1)
  selected.row <- node.indices[selected.node]

  child.nodes <- ED[selected.row, 3:4]
  child.rows <- node.indices[child.nodes]
  if (!all((is.na(ED[child.rows,4])) & (!is.na(ED[child.rows,3])))){
    #REJECT as cannot delete coalescent node
    return(list(ED = ED, prop.ratio = 0, log.prop.ratio = -Inf, node.indices = node.indices))
  }

  child.child <- ED[child.rows, 3]
  child.child.rows <- node.indices[child.child]
  exterior.demes <- ED[child.child.rows, 5]

  if (exterior.demes[1] != exterior.demes[2]){
    #REJECT as exterior demes are not the same
    return(list(ED = ED, prop.ratio = 0, log.prop.ratio = -Inf, node.indices = node.indices))
  }

  if (selected.node == root.node){
    prop.ratio <- 1 / ((n.deme - 1)* (ED[child.child.rows[1], 6] - ED[selected.row, 6]) * (ED[child.child.rows[2], 6] - ED[selected.row, 6]) )
    log.prop.ratio <- -log(n.deme - 1) - log(abs(ED[child.child.rows[1], 6] - ED[selected.row, 6])) - log(abs(ED[child.child.rows[2], 6] - ED[selected.row, 6]))

    ED[selected.row, c(3:5)] <- c(child.child, exterior.demes[1])
    ED[child.child.rows, 2] <- selected.node
    ED <- ED[!ED[,1] %in% child.nodes, ] #Remove child migration nodes

    node.indices[child.nodes] <- 0
    index.changes.1 <- (node.indices > child.rows[1])
    index.changes.2 <- (node.indices > child.rows[2])
    node.indices[index.changes.1] <- node.indices[index.changes.1] - 1
    node.indices[index.changes.2] <- node.indices[index.changes.2] - 1

    return(list(ED = ED, prop.ratio = prop.ratio, log.prop.ratio = log.prop.ratio, node.indices = node.indices))
  } else{
    parent.node <- ED[selected.row, 2]
    parent.row <- node.indices[parent.node]

    new.node <- max(ED[,1]) + 1
    new.node.time <- runif(1, ED[parent.row, 6], ED[selected.row, 6])
    old.deme <- ED[selected.row, 5]

    ED[selected.row, c(2:5)] <- c(new.node, child.child, exterior.demes[1])
    ED[c(child.child.rows[1], child.child.rows[2]), 2] <- selected.node

    which.child <- which(ED[parent.row, 3:4] == selected.node)
    ED[parent.row, 2 + which.child] <- new.node

    prop.ratio <- (ED[selected.row, 6] - ED[parent.row ,6]) / ((ED[child.child.rows[1], 6] - ED[selected.row, 6]) * (ED[child.child.rows[2], 6] - ED[selected.row, 6]))
    log.prop.ratio <- log(abs(ED[selected.row, 6] - ED[parent.row ,6])) - log(abs(ED[child.child.rows[1], 6] - ED[selected.row, 6])) - log(abs(ED[child.child.rows[2], 6] - ED[selected.row, 6]))


    ED <- ED[!ED[,1] %in% child.nodes, ] #Remove child migration nodes
    ED <- rbind(ED, c(new.node, parent.node, selected.node, NA, old.deme, new.node.time)) #Add new parent migration node

    node.indices[child.nodes] <- 0
    index.changes.1 <- (node.indices > child.rows[1])
    index.changes.2 <- (node.indices > child.rows[2])
    node.indices[index.changes.1] <- node.indices[index.changes.1] - 1
    node.indices[index.changes.2] <- node.indices[index.changes.2] - 1

    if (new.node > length(node.indices)){
      new.node.indices <- numeric(new.node)
      new.node.indices[1:length(node.indices)] <- node.indices
    } else{
      new.node.indices <- node.indices
    }
    new.node.indices[new.node] <- dim(ED)[1]

    return(list(ED = ED, prop.ratio = prop.ratio, log.prop.ratio = log.prop.ratio, node.indices = new.node.indices))
  }
}


#' Coalescent Node Split Proposal
#'
#' Performs a coalescent node split move (Ewing et al. 2004). Splits a migration
#' node immediately above a coalescent node into two migration nodes placed immediately
#' below the coalescent node
#'
#' @param ED Extended data object; matrix with columns Node ID, parent, child 1, child 2, deme, node age
#' @param n.deme Number of distinct demes in the population
#' @param node.indices Vector giving row indices for node labels
#'
#' @return Updated extended data object with the proposal from the migration pair birth move
#'
#' @export

ed.coal.split <- function(ED, n.deme, node.indices){
  root.node <- ED[is.na(ED[,2]), 1]
  coalescence.nodes <- ED[!is.na(ED[,4]),1]

  selected.node <- sample.vector(coalescence.nodes, 1)
  selected.row <- node.indices[selected.node]

  if (selected.node == root.node){
    child.nodes <- ED[selected.row, 3:4]
    child.rows <- node.indices[child.nodes]

    new.nodes <- max(ED[,1]) + c(1,2)
    new.node.times <- sapply(ED[child.rows, 6], runif, n = 1, min = ED[selected.row, 6])

    proposal.deme <- sample.vector((1:n.deme)[- ED[selected.row, 5]], 1)  #Sample new deme not equal to current deme of the root

    prop.ratio <- (n.deme - 1)* (ED[child.rows[1], 6] - ED[selected.row, 6]) * (ED[child.rows[2], 6] - ED[selected.row, 6])
    log.prop.ratio <- log(n.deme - 1) + log(abs(ED[child.rows[1], 6] - ED[selected.row, 6])) + log(abs(ED[child.rows[2], 6] - ED[selected.row, 6]))
    ED[selected.row, 3:5] <- c(new.nodes, proposal.deme)
    ED[child.rows, 2] <- new.nodes
    ED <- rbind(ED, matrix(c(new.nodes, rep(selected.node,2), child.nodes, rep(NA,2), rep(proposal.deme, 2), new.node.times), 2))

    if (new.nodes[2] > length(node.indices)){
      new.node.indices <- numeric(new.nodes[2])
      new.node.indices[1:length(node.indices)] <- node.indices
    } else{
      new.node.indices <- node.indices
    }
    new.node.indices[new.nodes] <- c(dim(ED)[1] - 1, dim(ED)[1])

    return(list(ED = ED, prop.ratio = prop.ratio, log.prop.ratio = log.prop.ratio, node.indices = new.node.indices))
  } else{
    parent.node <- ED[selected.row, 2]
    parent.row <- node.indices[parent.node]

    if (!((is.na(ED[parent.row,4])) & (!is.na(ED[parent.row,3])))){
      #REJECT
      return(list(ED = ED, prop.ratio = 0, log.prop.ratio = -Inf, node.indices = node.indices))
    }

    child.nodes <- ED[selected.row, 3:4]
    child.rows <- node.indices[child.nodes]

    parent.parent.node <- ED[parent.row, 2]
    parent.parent.row <- node.indices[parent.parent.node]
    new.nodes <- max(ED[,1]) + c(1,2)
    new.node.times <- sapply(ED[child.rows, 6], runif, n = 1, min = ED[selected.row, 6])

    prop.ratio <- ((ED[child.rows[1], 6] - ED[selected.row, 6]) * (ED[child.rows[2], 6] - ED[selected.row, 6])) / (ED[selected.row, 6] - ED[parent.parent.row ,6])
    log.prop.ratio <- log(abs(ED[child.rows[1], 6] - ED[selected.row, 6])) + log(abs(ED[child.rows[2], 6] - ED[selected.row, 6])) - log(abs(ED[selected.row, 6] - ED[parent.parent.row ,6]))

    which.child <- which(ED[parent.parent.row, 3:4] == parent.node)
    ED[parent.parent.row, 2 + which.child] <- selected.node
    ED[child.rows, 2] <- new.nodes
    ED[selected.row, 2:5] <- c(parent.parent.node, new.nodes, ED[parent.row, 5])
    ED <- ED[-parent.row, ]
    ED <- rbind(ED, matrix(c(new.nodes, rep(selected.node, 2), child.nodes, rep(NA, 2), rep(ED[selected.row, 5],2), new.node.times), 2))

    node.indices[parent.node] <- 0
    index.changes <- (node.indices > parent.row)
    node.indices[index.changes] <- node.indices[index.changes] - 1

    if (new.nodes[2] > length(node.indices)){
      new.node.indices <- numeric(new.nodes[2])
      new.node.indices[1:length(node.indices)] <- node.indices
    } else{
      new.node.indices <- node.indices
    }
    new.node.indices[new.nodes] <- c(dim(ED)[1] - 1, dim(ED)[1])

    return(list(ED = ED, prop.ratio = prop.ratio, log.prop.ratio = log.prop.ratio, node.indices = new.node.indices))
  }
}
