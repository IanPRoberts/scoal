#' Migration Pair Birth MCMC Move
#'
#' Performs a Migration pair birth move (Ewing et al. 2004). Adds two migration
#' nodes on an edge selected uniformly at random from a structured coalescent
#' process, allocating a deme for the added edge such that a migration event
#' does not target its origin deme.
#'
#' @param edge Edge matrix for a structured coalescent tree
#' @param edge.deme Label for the deme in which each edge lies
#' @param node.ages Time of each node in the coalescent tree
#' @param n.deme Number of possible demes in the process
#'
#' @return List containing the updated edge matrix,
#' node ages and edge demes

mig.pair.birth <- function(edge, edge.deme, node.ages, n.deme){
  all.nodes <- (unique(as.vector(edge)))

  #Update edge matrix
  selected.edge <- sample(seq_along(edge[,1]),1)  #Samples 1 row of edge matrix uniformly at random
  new.nodes <- max(all.nodes) + c(1,2)  #New migration node IDs to add
  new.node.ages <- sort(runif(2, min = min(node.ages[edge[selected.edge,]]), max = max(node.ages[edge[selected.edge,]]))) #Node ages distributed uniformly along selected edge
  edge <- rbind(edge, as.vector(new.nodes), c(edge[selected.edge,1],new.nodes[1]))  #Add new edges to edge matrix
  edge[selected.edge,1] <- new.nodes[2]  #Update remaining edge endpoints


  #Update edge deme labels
  selected.deme <- sample.vector((1:n.deme)[-edge.deme[selected.edge]],1)
  edge.deme <- c(edge.deme, selected.deme, edge.deme[selected.edge])

  #Update node ages
  node.ages <- c(node.ages, new.node.ages)

  output <- list()
  output$edge <- edge
  output$node.ages <- node.ages
  output$edge.deme <- edge.deme

  return(output)
}

#' Migration Pair Death MCMC Move
#'
#' Performs a Migration pair death move (Ewing et al. 2004). Deletes two
#' migration nodes if they lie between two edges in the same deme.
#'
#' @param edge Edge matrix for a structured coalescent tree
#' @param edge.deme Label for the deme in which each edge lies
#' @param node.ages Time of each node in the coalescent tree
#' @param n.deme Number of possible demes in the process
#'
#' @return List containing the updated edge matrix,
#' node ages and edge demes

mig.pair.death <- function(edge, edge.deme, node.ages, n.deme){
  # Identify node types based on frequencies in edge matrix
    # Leaves have freq 1
    # Root and migration nodes have freq 2 (Root is n.leaf + 1)
    # Coal nodes have freq 3

  all.nodes <- sort(unique(as.vector(edge)))
  node.freq <- table(match(as.vector(edge),all.nodes))  #Frequency of each node featured in edge matrix

  leaf.nodes <- all.nodes[node.freq == 1]  #Leaf nodes have frequency 1
  n.leaf <- length(leaf.nodes)
  root.node <- n.leaf + 1  #Root node is n.leaf + 1 for ape phylo-style object
  coalescence.nodes <- all.nodes[node.freq == 3]  #Coalescence nodes have frequency 3
  migration.nodes <- all.nodes[! all.nodes %in% c(root.node, coalescence.nodes, leaf.nodes)]

  selected.edge <- sample(seq_along(edge[,1]),1)  #Samples 1 row of edge matrix uniformly at random

  #Assuming selected.edge is a migration node, there will be one child and one parent edge from selected.edge
  parent.edge <- which(edge[,2] == edge[selected.edge,1])
  child.edge <- which(edge[,1] == edge[selected.edge,2])

  if ((edge[selected.edge,1] %in% migration.nodes == TRUE) &&
      (edge[selected.edge,2] %in% migration.nodes == TRUE) &&
      (edge.deme[child.edge] == edge.deme[parent.edge])){
    edge[child.edge,1] <- edge[parent.edge,1]  #Update remaining edge in edge matrix
    edge.deme <- edge.deme[-c(selected.edge, parent.edge)]
    node.ages <- node.ages[-edge[selected.edge,]]

    edge <- edge[-c(selected.edge, parent.edge),]  #Remove selected.edge and parent.edge from edge matrix

    output <- list()
    output$edge <- edge
    output$node.ages <- node.ages
    output$edge.deme <- edge.deme

    return(output)
  } else{
    return("REJECT")
  }
}

#' Migration Birth MCMC Move
#'
#' Performs a Migration birth move (Ewing et al. 2004). Adds a migration
#' node between a randomly selected ancestral node and its parent, allocating a
#' deme for the new edge consistent with the surrounding edges.
#'
#'
#' @param edge Edge matrix for a structured coalescent tree
#' @param edge.deme Label for the deme in which each edge lies
#' @param node.ages Time of each node in the coalescent tree
#' @param n.deme Number of possible demes in the process
#'
#' @return List containing the updated edge matrix,
#' node ages and edge demes

mig.birth <- function(edge, edge.deme, node.ages, n.deme){
  all.nodes <- sort(unique(as.vector(edge)))
  node.freq <- table(match(as.vector(edge),all.nodes))  #Frequency of each node featured in edge matrix

  leaf.nodes <- all.nodes[node.freq == 1]  #Leaf nodes have frequency 1
  n.leaf <- length(leaf.nodes)
  root.node <- n.leaf + 1  #Root node is n.leaf + 1 for ape phylo-style object
  coalescence.nodes <- all.nodes[node.freq == 3]  #Coalescence nodes have frequency 3
  migration.nodes <- all.nodes[! all.nodes %in% c(root.node, coalescence.nodes, leaf.nodes)]

  selected.node <- sample(all.nodes[-leaf.nodes],1)  #Ancestral node selected uniformly at random
  node.parent <- parent.node(selected.node, edge, node.ages)

  if (selected.node == root.node){
    return("REJECT")
  } else{
    new.node <- max(all.nodes) + 1
    parent.edge <- get.edge.id(selected.node,node.parent,edge)
    subtree <- list()  #List containing edges and nodes in maximal connected subtree containing <selected.node, new.node>
    subtree$nodes <- c(selected.node, new.node)
    subtree$edges <- parent.edge  #Edge <selected.node, node.parent> will be modified to <selected.node, new.node> later without changing edge ID if proposal is accepted
    active.nodes <- selected.node

    #Identify maximal connected subtree with tips at migration or leaf nodes
    while (TRUE %in% (active.nodes %in% coalescence.nodes)){
      active.nodes <- active.nodes[active.nodes %in% coalescence.nodes]
      #Add children of coalescence nodes to subtree
      for (i in active.nodes){
        i.children <- child.nodes(i,edge,node.ages)
        subtree$nodes <- c(subtree$nodes, i.children)
        for (j in i.children){
          subtree$edges <- c(subtree$edges, get.edge.id(i,j,edge))
        }
        active.nodes <- c(active.nodes[-1], i.children)  #[-1] removes first element
      }
    }
    #Reject if a leaf node lies within the subtree
    if (TRUE %in% (subtree$nodes %in% leaf.nodes)){
      return("REJECT")
    } else{
      subtree.leaf <- subtree$nodes[subtree$nodes %in% migration.nodes]
      exterior.deme <- edge.deme[parent.edge]
      for (i in subtree.leaf){
        terminal.child <- child.nodes(i,edge, node.ages)  #Identify child of each subtree leaf outside of subtree
        external.edge <- get.edge.id(i,terminal.child, edge)  #Edge ID of
        exterior.deme <- c(exterior.deme, edge.deme[external.edge])
      }

      exterior.deme <- unique(exterior.deme)

      if (length(exterior.deme) == n.deme){
        return("REJECT")
      } else{
        interior.deme <- edge.deme[parent.edge]
        new.deme <- sample.vector((1:n.deme)[-exterior.deme],1)

        # Update edge matrix
        edge[parent.edge,which(edge[parent.edge,] == node.parent)] <- new.node
        edge <- rbind(edge, c(node.parent, new.node))

        # Update edge deme
        edge.deme <- c(edge.deme, interior.deme)
        edge.deme[subtree$edges] <- new.deme

        # Add node age for new node
        node.ages <- c(node.ages, runif(1, min = node.ages[node.parent], max = node.ages[selected.node]))

        output <- list()
        output$edge <- edge
        output$node.ages <- node.ages
        output$edge.deme <- edge.deme

        return(output)
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
#' @param edge Edge matrix for a structured coalescent tree
#' @param edge.deme Label for the deme in which each edge lies
#' @param node.ages Time of each node in the coalescent tree
#' @param n.deme Number of possible demes in the process
#'
#' @return List containing the updated edge matrix,
#' node ages and edge demes

mig.death <- function(edge, edge.deme, node.ages, n.deme){
  all.nodes <- sort(unique(as.vector(edge)))
  node.freq <- table(match(as.vector(edge),all.nodes))  #Frequency of each node featured in edge matrix

  leaf.nodes <- all.nodes[node.freq == 1]  #Leaf nodes have frequency 1
  n.leaf <- length(leaf.nodes)
  root.node <- n.leaf + 1  #Root node is n.leaf + 1 for ape phylo-style object
  coalescence.nodes <- all.nodes[node.freq == 3]  #Coalescence nodes have frequency 3
  migration.nodes <- all.nodes[! all.nodes %in% c(root.node, coalescence.nodes, leaf.nodes)]

  selected.node <- sample(all.nodes[-leaf.nodes],1)  #Ancestral node selected uniformly at random
  node.parent <- parent.node(selected.node, edge, node.ages)
  node.parent2 <- parent.node(node.parent, edge, node.ages)

  if ((is.finite(node.parent) == FALSE) ||
      (is.finite(node.parent2) == FALSE) ||
      !(node.parent %in% migration.nodes)){
    # || necessary here when node.parent is root. || exits logical test without checking second condition if first condition is TRUE. | returns logical(0) if node.parent is root as node.parent2 = numeric()
    # Require node.parent to be migration node for edge <selected.node, node.parent2> to exist in the tree
    return("REJECT")
  } else{
    #Identify maximal subtree containing edge <selected.node, node.parent>
    subtree <- list()
    subtree$root <- node.parent2

    while (! subtree$root %in% c(migration.nodes, root.node)){
      #Sweep upwards through tree until subtree root is either the root of the tree, or a migration node
      subtree$root <- parent.node(subtree$root, edge, node.ages)
    }

    subtree$nodes <- subtree$root
    subtree$edges <- integer()
    subtree$leaves <- integer()
    active.nodes <- subtree$root

    while (FALSE %in% (active.nodes %in% c(migration.nodes[! migration.nodes %in% c(subtree$root, node.parent)], leaf.nodes))){
      #Sweep downwards from subtree root until hitting migration or leaf nodes. Ignore node.parent as a migration node
      new.active <- integer()
      for (i in active.nodes){
        i.children <- child.nodes(i, edge, node.ages)
        subtree$nodes <- c(subtree$nodes, i.children)

        for (j in i.children){
          subtree$edges <- c(subtree$edges, get.edge.id(i,j,edge))
          if (j %in% c(migration.nodes[migration.nodes != node.parent], leaf.nodes)){
            subtree$leaves <- c(subtree$leaves, j)
          } else{
            new.active <- c(new.active, j)
          }
        }
        active.nodes <- new.active
      }
    }

    if (TRUE %in% (subtree$nodes %in% leaf.nodes)){
      return("REJECT")
    } else{
      #Update subtree deme labels
      exterior.deme <- edge.deme[get.edge.id(subtree$root, parent.node(subtree$root, edge, node.ages), edge)]
      for (i in subtree$leaves){
        terminal.child <- child.nodes(i,edge, node.ages)  #Identify child of each subtree leaf outside of subtree
        external.edge <- get.edge.id(i,terminal.child, edge)  #Edge ID of
        exterior.deme <- c(exterior.deme, edge.deme[external.edge])
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
