#' Migration Birth MCMC Move
#'
#' Selections a location uniformly over a structured coalescent genealogy and
#' attempts to add a migration node, reallocating demes to obtain a consistent
#' migration history
#'
#'
#' @param ED Extended data object; matrix with columns Node ID, parent, child 1, child 2, deme, node age
#' @param n.deme Number of distinct demes in the population
#'
#' @return Updated extended data object with the proposal from the migration birth move
#'
#' @export

ed.mig.birth.new <- function(ED, n.deme){
  root.node <- ED[is.na(ED[,2]),1]
  coalescence.nodes <- ED[!is.na(ED[,4]),1]
  coalescence.nodes <- coalescence.nodes[coalescence.nodes != root.node]
  migration.nodes <- ED[is.na(ED[,4]) & (! is.na(ED[,3])), 1]
  leaf.nodes <- ED[is.na(ED[,3]),1]

  edge.length <- numeric(dim(ED)[1] - 1)

  for (i in (1:dim(ED)[1])[-root.node]){
    edge.length[i] <- ED[i,6] - ED[ED[i,2],6]
  }

  tree.length <- sum(edge.length)
  new.position <- runif(1, 0, tree.length)

  edge.child <- min(which(cumsum(edge.length) >= new.position))
  edge.parent <- which(ED[,1] == ED[edge.child, 2])

  subtree.nodes <- ED[edge.child,1]
  active.nodes <- ED[edge.child,1]
  while(TRUE %in% (active.nodes %in% coalescence.nodes)){
    active.nodes <- active.nodes[active.nodes %in% coalescence.nodes]

    for (i in active.nodes){
      i.row <- which(ED[,1] == i)
      if (is.na(ED[i.row,4]) == TRUE){
        if (is.na(ED[i.row,3]) == TRUE){  #Leaf nodes
          active.nodes <- active.nodes[active.nodes != i]
        } else{  #Migration nodes
          active.nodes <- c(active.nodes[active.nodes != i], ED[i.row, 3])
          subtree.nodes <- c(subtree.nodes, ED[i.row, 3])
        }
      } else{  #Coalescence nodes
        active.nodes <- c(active.nodes[active.nodes != i], ED[i.row, 3:4])
        subtree.nodes <- c(subtree.nodes, ED[i.row, 3:4])
      }
    }
  }

  if (TRUE %in% (subtree.nodes %in% leaf.nodes)){
    return("REJECT")
  } else{
    subtree.leaf <- subtree.nodes[subtree.nodes %in% migration.nodes]
    forbidden.deme <- ED[edge.child,5]  #Deme of parent edge

    for (i in subtree.leaf){
      leaf.row <- which(ED[,1] == i)
      forbidden.deme <- c(forbidden.deme, ED[which(ED[,1] == ED[leaf.row, 3]), 5])
    }

    forbidden.deme <- unique(forbidden.deme)

    if (length(forbidden.deme) == n.deme){
      return("REJECT")
    } else{
      new.node <- max(ED[,1]) + 1
      new.row <- dim(ED)[1] + 1
      ED <- rbind(ED, c(new.node, rep(NA,5)))
      ED[new.row, 6] <- ED[edge.child, 6] - (new.position - sum(edge.length[1:(edge.child-1)])) #Update new.node age
      ED[new.row, 5] <- ED[edge.child, 5]
      ED[new.row, 3] <- ED[edge.child, 1]
      ED[new.row, 2] <- ED[edge.parent, 1]

      ED[edge.child, 2] <- new.node

      new.deme <- sample.vector((1:n.deme)[-forbidden.deme],1)
      for (i in subtree.nodes){
        ED[which(ED[,1] == i),5] <- new.deme
      }

      ED[edge.parent, which(ED[edge.parent, 1:4] == edge.child)] <- new.node #Update child of parent.node
      return(ED)
    }
  }
}

#' Migration Death MCMC Move
#'
#' Performs a Migration death move. Selects a migration event
#' uniformly at random and removes the parent of the selected node, pulling the
#' deme from the parent edge downwards if a consistent deme labelling can be
#' obtained.
#'
#'
#' @param ED Extended data object; matrix with columns Node ID, parent, child 1, child 2, deme, node age
#' @param n.deme Number of distinct demes in the population
#'
#' @return Updated extended data object with the proposal from the migration death move
#'
#' @export

ed.mig.death.new <- function(ED, n.deme){
  root.node <- ED[is.na(ED[,2]),1]
  coalescence.nodes <- ED[!is.na(ED[,4]),1]
  coalescence.nodes <- coalescence.nodes[coalescence.nodes != root.node]
  migration.nodes <- ED[is.na(ED[,4]) & (! is.na(ED[,3])), 1]
  leaf.nodes <- ED[is.na(ED[,3]),1]

  selected.node <- sample.vector(migration.nodes, 1)  #Sample migration node
  selected.row <- which(ED[,1] == selected.node)  #Row in ED in case of non-sequential node labels
  #node.parent <- ED[selected.row, 2]
  #parent.row <- which(ED[,1] == node.parent)

  node.child <- ED[selected.node, 3] #Child of selected.node
  subtree.nodes <- c(selected.node, node.child)
  active.nodes <- c(selected.node, node.child)

  while (TRUE %in% (active.nodes %in% coalescence.nodes)){
    active.nodes <- active.nodes[active.nodes %in% coalescence.nodes]

    for (i in active.nodes){
      i.row <- which(ED[,1] == i)
      if (is.na(ED[i.row,4]) == TRUE){
        if (is.na(ED[i.row,3]) == TRUE){  #Leaf nodes
          active.nodes <- active.nodes[active.nodes != i]
        } else{  #Migration nodes
          active.nodes <- c(active.nodes[active.nodes != i], ED[i.row, 3])
          subtree.nodes <- c(subtree.nodes, ED[i.row, 3])
        }
      } else{  #Coalescence nodes
        active.nodes <- c(active.nodes[active.nodes != i], ED[i.row, 3:4])
        subtree.nodes <- c(subtree.nodes, ED[i.row, 3:4])
      }
    }
  }

    if (TRUE %in% (subtree.nodes %in% leaf.nodes)){
      return("REJECT")
    } else{
      #Check whether proposal.deme provides a consistent deme labelling
      proposal.deme <- ED[selected.row, 5]  #Deme of edge above selected.node
      exterior.deme <- numeric(0)
      subtree.leaf <- subtree.nodes[subtree.nodes %in% migration.nodes[migration.nodes != selected.node]]
      for (i in subtree.leaf){
        leaf.row <- which(ED[,1] == i)
        exterior.deme <- c(exterior.deme, ED[which(ED[,1] == ED[leaf.row, 3]), 5])
      }

      if (proposal.deme %in% exterior.deme){
        #Proposal would introduce a migration from a deme into itself -> reject proposal
        return("REJECT")
      } else{
        node.parent <- ED[selected.row, 2]  #Parent of selected.node
        parent.row <- which(ED[,1] == node.parent)

        child.row <- which(ED[,1] == node.child)

        ED[parent.row, which(ED[parent.row, 1:4] == selected.node)] <- node.child  #Replace selected.node in parent.row with child of selected.node
        ED[child.row, 2] <- node.parent  #Replace parent of node.child with parent of selected.node

        for (i in subtree.nodes){  #Update deme labels within subtree
          ED[which(ED[,1] == i),5] <- proposal.deme
        }

        ED <- ED[-selected.row,]
        return(ED)
      }

    }
  }
}
