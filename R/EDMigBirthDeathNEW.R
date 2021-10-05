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
    edge.length[i] <- ED[i,6] - ED[which(ED[,1] == ED[i,2]),6]
  }

  tree.length <- sum(edge.length)
  new.position <- runif(1, 0, tree.length)

  child.row <- min(which(cumsum(edge.length) >= new.position))
  parent.row <- which(ED[,1] == ED[child.row, 2])

  subtree.nodes <- ED[child.row,1]
  active.nodes <- ED[child.row,1]
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
    return(list(ED = ED, prop.ratio = 0))
  } else{
    subtree.leaf <- subtree.nodes[subtree.nodes %in% migration.nodes]
    forbidden.deme <- ED[child.row,5]  #Deme of subtree interior

    for (i in subtree.leaf){ #Demes adjacent to tips of subtree
      leaf.row <- which(ED[,1] == i)
      forbidden.deme <- c(forbidden.deme, ED[which(ED[,1] == ED[leaf.row, 3]), 5])
    }

    forbidden.deme <- unique(forbidden.deme)

    if (length(forbidden.deme) == n.deme){
      return(list(ED = ED, prop.ratio = 0))
    } else{
      new.node <- max(ED[,1]) + 1
      new.row <- dim(ED)[1] + 1
      ED <- rbind(ED, c(new.node, rep(NA,5)))
      ED[new.row, 6] <- ED[child.row, 6] - (new.position - sum(edge.length[1:(child.row-1)])) #Update new.node age
      ED[new.row, 5] <- ED[child.row, 5]
      ED[new.row, 3] <- ED[child.row, 1]
      ED[new.row, 2] <- ED[parent.row, 1]

      ED[child.row, 2] <- new.node

      new.deme <- sample.vector((1:n.deme)[-forbidden.deme],1)
      for (i in subtree.nodes){
        ED[which(ED[,1] == i),5] <- new.deme
      }

      ED[parent.row, which(ED[parent.row, 1:4] == child.row)] <- new.node #Update child of parent.node

      prop.ratio <- tree.length * (n.deme - length(forbidden.deme)) / (1 + length(migration.nodes))  #Proposal ratio
      return(list(ED = ED, prop.ratio = prop.ratio))
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

  node.child <- ED[selected.row, 3] #Child of selected.node
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
      return(list(ED = ED, prop.ratio = 0))
    } else{
      #Check whether proposal.deme provides a consistent deme labelling
      proposal.deme <- ED[selected.row, 5]  #Deme of edge above selected.node
      exterior.deme <- numeric(0)
      subtree.leaf <- subtree.nodes[subtree.nodes %in% migration.nodes[migration.nodes != selected.node]]
      for (i in subtree.leaf){
        leaf.row <- which(ED[,1] == i)
        exterior.deme <- c(exterior.deme, ED[which(ED[,1] == ED[leaf.row, 3]), 5])
      }
      exterior.deme <- unique(exterior.deme)

      if (proposal.deme %in% exterior.deme){
        #Proposal would introduce a migration from a deme into itself -> reject proposal
        return(list(ED = ED, prop.ratio = 0))
      } else{

        #Tree length needed for proposal ratio
        edge.length <- numeric(dim(ED)[1] - 1)
        for (i in (1:dim(ED)[1])[-root.node]){
          edge.length[i] <- ED[i,6] - ED[which(ED[,1] == ED[i,2]),6]
        }
        tree.length <- sum(edge.length)

        node.parent <- ED[selected.row, 2]  #Parent of selected.node
        parent.row <- which(ED[,1] == node.parent)

        child.row <- which(ED[,1] == node.child)

        ED[parent.row, which(ED[parent.row, 1:4] == selected.node)] <- node.child  #Replace selected.node in parent.row with child of selected.node
        ED[child.row, 2] <- node.parent  #Replace parent of node.child with parent of selected.node

        for (i in subtree.nodes){  #Update deme labels within subtree
          ED[which(ED[,1] == i),5] <- proposal.deme
        }

        ED <- ED[-selected.row,]

        prop.ratio <- length(migration.nodes) / (tree.length * (n.deme - length(exterior.deme) - 1))

        return(list(ED = ED, prop.ratio = prop.ratio))
      }

    }
  }
