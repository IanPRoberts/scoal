#' Block Recolouring MCMC Proposal
#'
#' Performs a block recolouring MCMC proposal on a structured coalescent
#' migration history. A block of the migration history currently in a single
#' deme is migrated into a new deme selected uniformly from all other demes. If
#' an inconsistent migration history is obtained, the move will be rejected.
#'
#' @param ED Extended data object; matrix with columns Node ID, parent, child 1, child 2, deme, node age
#' @param n.deme Number of distinct demes in the population
#' @param fix.leaf.deme Logical; if TRUE, the deme of each leaf in the tree is fixed. Any proposal which attempts to change a leaf deme is automatically rejected
#'
#' @return Updated extended data object with the proposal from the block recolouring proposal
#'
#' @export

ed.block.recolour <- function(ED, n.deme, fix.leaf.deme = TRUE, node.indices){
  root.node <- ED[is.na(ED[,2]), 1]
  leaf.nodes <- ED[(is.na(ED[,3])) & (is.na(ED[,4])), 1]
  migration.nodes <- ED[ is.na(ED[,4]) & (!is.na(ED[,3])) ,1]

  if (length(migration.nodes) == 0){
    if (fix.leaf.deme == TRUE){
      #REJECT - no migration nodes and fixed leaf demes -> cannot change deme
      return(list(ED = ED, prop.ratio = 0, log.prop.ratio = -Inf, node.indices = node.indices))
    } else{ #No migration nodes -> all same deme
      proposal.deme <- sample.vector((1:n.deme)[-ED[1, 5]], 1)
      ED[,5] <- proposal.deme
      return(list(ED = ED, prop.ratio = 1, log.prop.ratio = 0, node.indices = node.indices))
    }
  }

  parent.rows <- node.indices[ED[,2]]
  parent.times <- ED[parent.rows, 6]
  parent.times[root.node] <- 0
  edge.length <- ED[,6] - parent.times

  tree.length <- sum(edge.length)  #Total tree length
  new.location <- runif(1, 0, tree.length)  #Uniform location along tree.length

  child.row <- min(which(cumsum(edge.length) >= new.location))
  child.node <- ED[child.row, 1]

  #Sweep upwards until hitting the root or a migration node
  subtree.root <- ED[child.row, 2]
  while (! any(subtree.root == c(migration.nodes, root.node))){
    subtree.root.row <- node.indices[subtree.root]
    subtree.root <- ED[subtree.root.row, 2]
  }

  #Sweep downwards until hitting migration nodes
  subtree.root.row <- node.indices[subtree.root]
  if (subtree.root == root.node){
    active.nodes <- ED[subtree.root.row, 3:4]
    subtree.nodes <- c(subtree.root, active.nodes)
  } else{
    active.nodes <- ED[subtree.root.row, 3]
    subtree.nodes <- c(subtree.root, active.nodes)
  }

  # active.rows <- node.indices[active.nodes]
  # while (length(active.rows) > 0){
  #   subtree.nodes <- c(subtree.nodes, ED[active.rows, 3:4])
  #   active.nodes <- ED[active.rows, 3:4]
  #   active.rows <- node.indices[active.nodes]
  #   active.rows <- active.rows[!is.na(ED[active.rows, 4])]
  # }

  while (! all(active.nodes %in% migration.nodes)){
    active.nodes <- active.nodes[! (active.nodes %in% c(migration.nodes, leaf.nodes))]
    active.rows <- node.indices[active.nodes]
    subtree.nodes <- c(subtree.nodes, ED[active.rows, 3:4])
    active.nodes <- ED[active.rows, 3:4]
  }

  if ((fix.leaf.deme == TRUE) && (any(subtree.nodes %in% leaf.nodes))){
    # REJECT
    return(list(ED = ED, prop.ratio = 0, log.prop.ratio = -Inf, node.indices = node.indices))
  } else{
    #Continue proposal
    subtree.leaves <- subtree.nodes[subtree.nodes %in% migration.nodes]
    if (subtree.root %in% subtree.leaves){
      subtree.leaves <- subtree.leaves[subtree.leaves != subtree.root]
    }
    old.deme <- ED[child.row, 5]

    proposal.deme <- sample.vector((1:n.deme)[- old.deme], 1)  #Propose deme update without accounting for surrounding demes

    subtree.leaf.rows <- node.indices[subtree.leaves]
    subtree.leaf.children <- ED[subtree.leaf.rows, 3]
    subtree.leaf.children.rows <- node.indices[subtree.leaf.children]

    if (any(ED[subtree.leaf.children.rows, 5] == proposal.deme)){
      #REJECT - induces self-migration
      return(list(ED = ED, prop.ratio = 0, log.prop.ratio = -Inf, node.indices = node.indices))
    }

    if (subtree.root != root.node){
      if (ED[subtree.root.row, 5] == proposal.deme){
        #REJECT
        return(list(ED = ED, prop.ratio = 0, log.prop.ratio = -Inf, node.indices = node.indices))
      }
    }

    #Update deme across subtree
    if (subtree.root == root.node){
      subtree.rows <- node.indices[subtree.nodes]
    } else{
      subtree.rows <- node.indices[subtree.nodes[subtree.nodes != subtree.root]]
    }
    ED[subtree.rows, 5] <- proposal.deme

    return(list(ED = ED, prop.ratio = 1, log.prop.ratio = 0, node.indices = node.indices))
  }
}
