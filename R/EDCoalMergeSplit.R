### Coalescent node Merge Proposal

ed.coal.merge <- function(ED, n.deme){
  root.node <- ED[is.na(ED[,2]), 1]
  coalescence.nodes <- ED[!is.na(ED[,4]),1]
  migration.nodes <- ED[ is.na(ED[,4]) & (!is.na(ED[,3])) ,1]

  selected.node <- sample.vector(coalescence.nodes, 1)
  selected.row <- which(ED[,1] == selected.node)

  child.node.1 <- ED[selected.row, 3]
  child.node.2 <- ED[selected.row, 4]

  if ((!child.node.1 %in% migration.nodes) || (! child.node.2 %in% migration.nodes)){
    #REJECT as cannot delete coalescent node
    return(list(ED = ED, prop.ratio = 0, log.prop.ratio = -Inf))
  }

  child.row.1 <- which(ED[,1] == child.node.1)
  child.child.1 <- ED[child.row.1, 3]
  child.child.1.row <- which(ED[,1] == child.child.1)
  exterior.deme.1 <- ED[child.child.1.row,5] #Deme below child 1
  child.row.2 <- which(ED[,1] == child.node.2)
  child.child.2 <- ED[child.row.2, 3]
  child.child.2.row <- which(ED[,1] == child.child.2)
  exterior.deme.2 <- ED[child.child.2.row,5] #Deme below child 2

  if (exterior.deme.1 != exterior.deme.2){
    #REJECT as exterior demes are not the same
    return(list(ED = ED, prop.ratio = 0, log.prop.ratio = -Inf))
  }

  if (selected.node == root.node){
    prop.ratio <- 1 / ((n.deme - 1)* (ED[child.child.1.row, 6] - ED[selected.row, 6]) * (ED[child.child.2.row, 6] - ED[selected.row, 6]) )
    log.prop.ratio <- -log(n.deme - 1) - log(abs(ED[child.child.1.row, 6] - ED[selected.row, 6])) - log(abs(ED[child.child.2.row, 6] - ED[selected.row, 6]))

    ED[selected.row, c(3:5)] <- c(child.child.1, child.child.2, exterior.deme.1)
    ED[c(child.child.1.row, child.child.2.row), 2] <- selected.row
    ED <- ED[!ED[,1] %in% c(child.node.1, child.node.2), ] #Remove child migration nodes

    return(list(ED = ED, prop.ratio = prop.ratio, log.prop.ratio = log.prop.ratio))
  } else{
    parent.node <- ED[selected.row, 2]
    parent.row <- which(ED[,1] == parent.node)

    new.node <- max(ED[,1]) + 1
    new.node.time <- runif(1, ED[parent.row, 6], ED[selected.row, 6])
    old.deme <- ED[selected.row, 5]

    ED[selected.row, c(2:5)] <- c(new.node, child.child.1, child.child.2, exterior.deme.1)
    ED[c(child.child.1.row, child.child.2.row), 2] <- selected.node

    which.child <- which(ED[parent.row, 3:4] == selected.node)
    ED[parent.row, 2 + which.child] <- new.node

    prop.ratio <- (ED[selected.row, 6] - ED[parent.row ,6]) / ((ED[child.child.1.row, 6] - ED[selected.row, 6]) * (ED[child.child.2.row, 6] - ED[selected.row, 6]))
    log.prop.ratio <- log(abs(ED[selected.row, 6] - ED[parent.row ,6])) - log(abs(ED[child.child.1.row, 6] - ED[selected.row, 6])) - log(abs(ED[child.child.2.row, 6] - ED[selected.row, 6]))


    ED <- ED[!ED[,1] %in% c(child.node.1, child.node.2), ] #Remove child migration nodes
    ED <- rbind(ED, c(new.node, parent.node, selected.node, NA, old.deme, new.node.time)) #Add new parent migration node

    return(list(ED = ED, prop.ratio = prop.ratio, log.prop.ratio = log.prop.ratio))
  }
}


### Coalescent node split proposal

ed.coal.split <- function(ED, n.deme){
  root.node <- ED[is.na(ED[,2]), 1]
  coalescence.nodes <- ED[!is.na(ED[,4]),1]
  migration.nodes <- ED[ is.na(ED[,4]) & (!is.na(ED[,3])) ,1]

  selected.node <- sample.vector(coalescence.nodes, 1)
  selected.row <- which(ED[,1] == selected.node)

  if (selected.node == root.node){
    child.node.1 <- ED[selected.row, 3]
    child.row.1 <- which(ED[,1] == child.node.1)
    child.node.2 <- ED[selected.row, 4]
    child.row.2 <- which(ED[,1] == child.node.2)

    new.nodes <- max(ED[,1]) + c(1,2)
    new.node.time.1 <- runif(1, ED[selected.row, 6], ED[child.row.1, 6])
    new.node.time.2 <- runif(1, ED[selected.row, 6], ED[child.row.2, 6])

    proposal.deme <- sample.vector((1:n.deme)[- ED[selected.row, 5]], 1)  #Sample new deme not equal to current deme of the root

    prop.ratio <- (n.deme - 1)* (ED[child.row.1, 6] - ED[selected.row, 6]) * (ED[child.row.2, 6] - ED[selected.row, 6])
    log.prop.ratio <- log(n.deme - 1) + log(abs(ED[child.row.1, 6] - ED[selected.row, 6])) + log(abs(ED[child.row.2, 6] - ED[selected.row, 6]))
    ED[selected.row, 3:5] <- c(new.nodes, proposal.deme)
    ED[child.row.1, 2] <- new.nodes[1]
    ED[child.row.2, 2] <- new.nodes[2]
    ED <- rbind(ED,
                c(new.nodes[1], selected.node, child.node.1, NA, proposal.deme, new.node.time.1),
                c(new.nodes[2], selected.node, child.node.2, NA, proposal.deme, new.node.time.2))
    return(list(ED = ED, prop.ratio = prop.ratio, log.prop.ratio = log.prop.ratio))
  } else{
    parent.node <- ED[selected.row, 2]

    if (!parent.node %in% migration.nodes){
      #REJECT
      return(list(ED = ED, prop.ratio = 0, log.prop.ratio = -Inf))
    }

    parent.row <- which(ED[,1] == parent.node)
    child.node.1 <- ED[selected.row, 3]
    child.row.1 <- which(ED[,1] == child.node.1)
    child.node.2 <- ED[selected.row, 4]
    child.row.2 <- which(ED[,1] == child.node.2)

    parent.parent.node <- ED[parent.row, 2]
    parent.parent.row <- which(ED[,1] == parent.parent.node)
    new.nodes <- max(ED[,1]) + c(1,2)
    new.node.time.1 <- runif(1, ED[selected.row, 6], ED[child.row.1, 6])
    new.node.time.2 <- runif(1, ED[selected.row, 6], ED[child.row.2, 6])

    prop.ratio <- ((ED[child.row.1, 6] - ED[selected.row, 6]) * (ED[child.row.2, 6] - ED[selected.row, 6])) / (ED[selected.row, 6] - ED[parent.parent.row ,6])
    log.prop.ratio <- log(abs(ED[child.row.1, 6] - ED[selected.row, 6])) + log(abs(ED[child.row.2, 6] - ED[selected.row, 6])) - log(abs(ED[selected.row, 6] - ED[parent.parent.row ,6]))

    which.child <- which(ED[parent.parent.row, 3:4] == parent.node)
    ED[parent.parent.row, 2 + which.child] <- selected.node
    ED[child.row.1, 2] <- new.nodes[1]
    ED[child.row.2, 2] <- new.nodes[2]
    ED[selected.row, 2:5] <- c(parent.parent.node, new.nodes, ED[parent.row, 5])
    ED <- ED[-parent.row, ]
    ED <- rbind(ED,
                c(new.nodes[1], selected.node, child.node.1, NA, ED[selected.node, 5], new.node.time.1),
                c(new.nodes[2], selected.node, child.node.2, NA, ED[selected.node, 5], new.node.time.2))

    return(list(ED = ED, prop.ratio = prop.ratio, log.prop.ratio = log.prop.ratio))
  }
}
