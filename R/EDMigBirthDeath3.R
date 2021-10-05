#### Redoing mig birth
ed.mig.birth.3 <- function(ED, n.deme){
  root.node <- ED[is.na(ED[,2]), 1]
  coalescence.nodes <- ED[!is.na(ED[,4]),1]
  coalescence.nodes <- coalescence.nodes[coalescence.nodes != root.node]
  leaf.nodes <- ED[(is.na(ED[,3])) & (is.na(ED[,4])), 1]
  migration.nodes <- ED[!ED[,1] %in% c(root.node, coalescence.nodes, leaf.nodes), 1]

  edge.length <- numeric(dim(ED)[1])
  for (i in c(coalescence.nodes, leaf.nodes, migration.nodes)){
    node.row <- which(ED[,1] == i)
    parent.row <- which(ED[,1] == ED[node.row, 2])
    edge.length[node.row] <- ED[node.row, 6] - ED[parent.row, 6]
  }

  tree.length <- sum(edge.length)  #Total tree length
  new.location <- runif(1, 0, tree.length)

  child.row <- min(which(cumsum(edge.length) >= new.location))
  child.node <- ED[child.row, 1]

  parent.node <- ED[child.row, 2]
  parent.row <- which(ED[,1] == parent.node)

  #Calculating subtree containing edge <new.node, child.node> terminating in migration or leaf nodes
  subtree.nodes <- child.node
  active.nodes <- child.node
  while (TRUE %in% (active.nodes %in% coalescence.nodes)){
    active.nodes <- active.nodes[active.nodes %in% coalescence.nodes]

    for (i in active.nodes){
      i.row <- which(ED[,1] == i)
      subtree.nodes <- c(subtree.nodes, ED[i.row, 3:4])  #Add children of i to subtree
      active.nodes <- c(active.nodes[active.nodes != i], ED[i.row, 3:4])  #Remove i from active nodes, add children
    }
  }

  if (TRUE %in% (subtree.nodes %in% leaf.nodes)){
    # REJECT
    return(list(ED = ED, prop.ratio = 0))
  } else{
    #Continue proposal
    subtree.leaves <- subtree.nodes[subtree.nodes %in% migration.nodes]
    forbidden.demes <- ED[child.row, 5]  #Cannot propose deme already on interior of subtree
    old.deme <- ED[child.row,5]
    for (i in subtree.leaves){  #Cannot propose deme already adjacent to subtree
      i.row <- which(ED[,1] == i)
      i.child <- ED[i.row, 3]
      i.child.row <- which(ED[,1] == i.child)
      forbidden.demes <- c(forbidden.demes, ED[i.child.row, 5])  #Deme below node i
    }
    forbidden.demes <- unique(forbidden.demes)

    if (length(forbidden.demes) == n.deme){
      #No valid deme for proposal
      #REJECT
      return(list(ED = ED, prop.ratio = 0))
    }

    proposal.deme <- sample.vector((1:n.deme)[-forbidden.demes], 1)
    #Update deme across subtree
    for (i in subtree.nodes){
      row <- which(ED[,1] == i)
      ED[row, 5] <- proposal.deme
    }

    new.node <- max(ED[,1]) + 1

    #Update child.row with new parent
    ED[child.row, 2] <- new.node

    #Update parent row with new child
    which.child <- which(ED[parent.row, 3:4] == child.node)
    ED[parent.row, 2 + which.child] <- new.node  #which.child either 1 or 2 so add 2 to get to col no.

    #Add row for new.node
    new.age <- ED[child.row, 6] - (new.location - max(cumsum(edge.length)[cumsum(edge.length) < new.location])) #Age of new.node (currently have distance along tree from new.location)

    ED <- rbind(ED, c(new.node, parent.node, child.node, NA, old.deme, new.age))

    return(list(ED= ED, prop.ratio = 1))
  }
}


#### Redoing Mig death
ed.mig.death.3 <- function(ED, n.deme){
  root.node <- ED[is.na(ED[,2]), 1]
  coalescence.nodes <- ED[!is.na(ED[,4]),1]
  coalescence.nodes <- coalescence.nodes[coalescence.nodes != root.node]
  leaf.nodes <- ED[(is.na(ED[,3])) & (is.na(ED[,4])), 1]
  migration.nodes <- ED[!ED[,1] %in% c(root.node, coalescence.nodes, leaf.nodes), 1]

  selected.node <- sample.vector(migration.nodes, 1)
  selected.row <- which(ED[,1] == selected.node)

  child.node <- ED[selected.row, 3]
  child.row <- which(ED[,1] == child.node)

  #Calculating subtree containing edge <new.node, child.node> terminating in migration or leaf nodes
  subtree.nodes <- child.node
  active.nodes <- child.node
  while (TRUE %in% (active.nodes %in% coalescence.nodes)){
    active.nodes <- active.nodes[active.nodes %in% coalescence.nodes]

    for (i in active.nodes){
      i.row <- which(ED[,1] == i)
      subtree.nodes <- c(subtree.nodes, ED[i.row, 3:4])  #Add children of i to subtree
      active.nodes <- c(active.nodes[active.nodes != i], ED[i.row, 3:4])  #Remove i from active nodes, add children
    }
  }

  if (TRUE %in% (subtree.nodes %in% leaf.nodes)){
    # REJECT
    return(list(ED = ED, prop.ratio = 0))
  } else{
    #Continue proposal
    subtree.leaves <- subtree.nodes[subtree.nodes %in% migration.nodes]
    forbidden.demes <- numeric(0)
    for (i in subtree.leaves){  #Cannot change to deme adjacent to subtree
      i.row <- which(ED[,1] == i)
      i.child <- ED[i.row, 3]
      i.child.row <- which(ED[,1] == i.child)
      forbidden.demes <- c(forbidden.demes, ED[i.child.row, 5])  #Deme below node i
    }
    forbidden.demes <- unique(forbidden.demes)

    proposal.deme <- ED[selected.row, 5]  #Propose deme immediately above subtree
    if (proposal.deme %in% forbidden.demes){
      #Cannot change deme
      #REJECT
      return(list(ED = ED, prop.ratio = 0))
    }

    parent.node <- ED[selected.row, 2]
    parent.row <- which(ED[,1] == parent.node)

    #Update deme across subtree
    for (i in subtree.nodes){
      i.row <- which(ED[,1] == i)
      ED[i.row, 5] <- proposal.deme
    }

    ED[child.row, 2] <- parent.node  #Update child.node parent

    #Update parent.node child
    which.child <- which(ED[parent.row, 3:4] == selected.node)
    ED[parent.row, 2 + which.child] <- child.node

    ED <- ED[-selected.row,]

    return(list(ED = ED, prop.ratio = 1))
  }
}
