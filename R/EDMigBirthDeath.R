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

ed.mig.birth <- function(ED, n.deme){
  all.nodes <- ED[,1]
  leaf.nodes <- all.nodes[is.na(ED[,3])]
  root.node <- all.nodes[is.na(ED[,2])]
  coalescence.nodes <- all.nodes[!is.na(ED[,4])]
  coalescence.nodes <- coalescence.nodes[coalescence.nodes != root.node]
  migration.nodes <- all.nodes[! all.nodes %in% c(root.node, leaf.nodes, coalescence.nodes)]

  selected.node <- sample(all.nodes[-leaf.nodes],1)  #Ancestral node selected uniformly at random
  selected.row <- which(ED[,1] == selected.node)
  node.parent <- ED[selected.row, 2]
  parent.row <- which(ED[,1] == node.parent)

  if (selected.node == root.node){
    return("REJECT")
  } else{
    new.node <- max(all.nodes) + 1
    subtree.nodes <- selected.node
    active.nodes <- selected.node

    while(TRUE %in% (active.nodes %in% coalescence.nodes)){
      active.nodes <- active.nodes[active.nodes %in% coalescence.nodes]

      #Add children of coalescence nodes to subtree
      #Root of subtree for a birth move is always the new node being added, so don't need to sweep upwards through tree
      for (i in active.nodes){
        i.row <- which(ED[,1] == i)
        if (is.na(ED[i.row,4]) == TRUE){
          if (is.na(ED[i.row,3]) == TRUE){
            active.nodes <- active.nodes[active.nodes != i]
          } else{
            active.nodes <- c(active.nodes[active.nodes != i], ED[i.row, 3])
            subtree.nodes <- c(subtree.nodes, ED[i.row, 3])
          }

        } else{
          active.nodes <- c(active.nodes[active.nodes != i], ED[i.row, 3:4])
          subtree.nodes <- c(subtree.nodes, ED[i.row, 3:4])
        }
      }
    }

    #Reject if a leaf node lies within the subtree
    if (TRUE %in% (subtree.nodes %in% leaf.nodes)){
      return("REJECT")
    } else{
      subtree.leaf <- subtree.nodes[subtree.nodes %in% migration.nodes]
      exterior.deme <- ED[selected.row,5]  #Deme of parent edge

      for (i in subtree.leaf){
        leaf.row <- which(ED[,1] == i)
        exterior.deme <- c(exterior.deme, ED[which(ED[,1] == ED[leaf.row, 3]), 5])
      }

      exterior.deme <- unique(exterior.deme)

      if (length(exterior.deme) == n.deme){
        return("REJECT")
      } else{
        interior.deme <- ED[selected.row,5]
        new.deme <- sample.vector((1:n.deme)[-exterior.deme],1)

        #Update ED for output
        new.age <- runif(1, ED[parent.row, 6], ED[selected.row, 6])
        ED <- rbind(ED, c(new.node, node.parent, selected.node, NA, interior.deme, new.age))

        for (i in subtree.nodes){
          ED[which(ED[,1] == i),5] <- new.deme
        }
        ED[selected.row,2] <- new.node #Update parent of selected.row
        ED[parent.row, which(ED[parent.row, ] == selected.node)] <- new.node #Update child of parent.node
        return(ED)
      }
    }
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


############### INCOMPLETE  ##################
ed.mig.death <- function(ED, n.deme){
  all.nodes <- ED[,1]
  leaf.nodes <- all.nodes[is.na(ED[,3])]
  root.node <- all.nodes[is.na(ED[,2])]
  coalescence.nodes <- all.nodes[!is.na(ED[,4])]
  coalescence.nodes <- coalescence.nodes[coalescence.nodes != root.node]
  migration.nodes <- all.nodes[! all.nodes %in% c(root.node, leaf.nodes, coalescence.nodes)]

  selected.node <- sample(all.nodes[-leaf.nodes],1)  #Ancestral node selected uniformly at random
  selected.row <- which(ED[,1] == selected.node)
  node.parent <- ED[selected.row, 2]
  parent.row <- which(ED[,1] == node.parent)
  node.parent2 <- ED[parent.row, 2]

  if ((is.na(node.parent) == TRUE) ||
      (is.na(node.parent2) == TRUE) ||
      !(node.parent %in% migration.nodes)){
    # || necessary here when node.parent is root. || exits logical test without checking second condition if first condition is TRUE. | returns logical(0) if node.parent is root as node.parent2 = numeric()
    # Require node.parent to be migration node for edge <selected.node, node.parent2> to exist in the tree
    return("REJECT")
  } else{
    subtree.root <- node.parent2

    #Identify subtree root by sweeping upwards until either the root or a migration node
    while (! subtree.root %in% c(migration.nodes, root.node)){
      subtree.root <- ED[which(ED[,1] == subtree.root),2]
    }

    subtree.nodes <- subtree.root
    active.nodes <- subtree.root

    #Sweep downwards from subtree root until hitting migration or leaf nodes (ignoring node.parent as a migration node)
    while (FALSE %in% (active.nodes %in% c(migration.nodes[! migration.nodes %in% c(subtree.root, node.parent)], leaf.nodes))){
      for (i in active.nodes){
        i.row <- which(ED[,1] == i)
        if (is.na(ED[i.row,4]) == TRUE){
          if (is.na(ED[i.row,3]) == TRUE){
            active.nodes <- active.nodes[active.nodes != i]
          } else{
            active.nodes <- c(active.nodes[active.nodes != i], ED[i.row, 3])
            subtree.nodes <- c(subtree.nodes, ED[i.row, 3])
          }

        } else{
          active.nodes <- c(active.nodes[active.nodes != i], ED[i.row, 3:4])
          subtree.nodes <- c(subtree.nodes, ED[i.row, 3:4])
        }
      }
    }

    if (TRUE %in% (subtree.nodes %in% leaf.nodes)){
      return("REJECT")
    } else{
      #Update subtree deme labels
      exterior.deme <- ED[which(ED[,1] == subtree.root), 5]

      subtree.leaf <- subtree.nodes[subtree.nodes %in% migration.nodes[! migration.nodes %in% c(node.parent, subtree.root)]]
      for (i in subtree.leaf){
        leaf.row <- which(ED[,1] == i)
        exterior.deme <- c(exterior.deme, ED[which(ED[,1] == ED[leaf.row, 3]), 5])
      }

      exterior.deme <- unique(exterior.deme)

      if (length(exterior.deme) == n.deme){
        return("REJECT")
      } else{
        interior.deme <- edge.deme[parent.edge]
        new.deme <- sample.vector((1:n.deme)[-exterior.deme],1)

        edge.deme[subtree$edges] <- new.deme
        parent.edge <- get.edge.id(selected.node, node.parent, edge)
        edge <- edge[-parent.edge,]
        edge.deme <- edge.deme[-parent.edge]
        node.ages <- node.ages[-node.parent]

        output <- list()
        output$edge <- edge
        output$node.ages <- node.ages
        output$edge.deme <- edge.deme

        return(output)
      }
    }
  }
}
