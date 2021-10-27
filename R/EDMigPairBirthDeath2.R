## New migration pair birth

ed.mig.pair.birth.2 <- function(ED, n.deme){
  root.node <- ED[is.na(ED[,2]), 1]
  migration.nodes <- ED[ED[ is.na(ED[,4]) & (!is.na(ED[,3])) ,1]]

  #Sample branch in genealogy proportional to length
  edge.length <- numeric(dim(ED)[1])
  non.root.nodes <- ED[ED[,1] != root.node, 1]
  non.root.nodes <- non.root.nodes[!is.na(non.root.nodes)]
  for (j in non.root.nodes){
    node.row <- which(ED[,1] == j)
    parent.row <- which(ED[,1] == ED[node.row, 2])
    edge.length[node.row] <- ED[node.row, 6] - ED[parent.row, 6]
  }
  tree.length <- sum(edge.length)  #Total tree length

  genealogy.loc <- runif(1, 0, tree.length)
  #Two nodes above and below genealogy.loc
  below.row <- min(which(cumsum(edge.length) >= genealogy.loc))
  below.node <- ED[below.row, 1]

  above.node <- ED[below.row, 2]
  above.row <- which(ED[,1] == above.node)

  branch.nodes <- c(below.node, above.node)
  branch.rows <- c(below.row, above.row)

  active.node <- above.node
  active.row <- above.row
  while (active.node %in% migration.nodes){
    active.node <- ED[active.row, 2]
    active.row <- which(ED[,1] == active.node)
    branch.nodes <- c(branch.nodes, active.node)
    branch.rows <- c(branch.rows, active.row)
  }
  branch.root <- active.node  #Root end of genealogy branch
  active.node <- below.node
  active.row <- below.row
  while (active.node %in% migration.nodes){
    active.node <- ED[active.row, 3]
    active.row <- which(ED[,1] == active.node)
    branch.nodes <- c(branch.nodes, active.node)
    branch.rows <- c(branch.rows, active.row)
  }
  branch.leaf <- active.node  #Leaf end of genealogy branch
  branch.nodes <- branch.nodes[order(ED[branch.rows, 6], decreasing = TRUE)]  #Ordering branch nodes/rows nicely (unecessary but easier for checking)
  branch.rows <- branch.rows[order(ED[branch.rows, 6], decreasing = TRUE)]

  branch.length <- ED[which(ED[,1] == branch.leaf), 6] - ED[which(ED[,1] == branch.root), 6]

  new.node.locs <- sort(runif(2, 0, branch.length), decreasing = TRUE) #Distance along branch from branch "leaf"

  #Identify nodes surrounding the new nodes
  branch.edge.length <- edge.length[branch.rows]
  new.child.1.row <- branch.rows[min(which(cumsum(branch.edge.length) >= new.node.locs[1]))]
  new.child.1 <- ED[new.child.1.row, 1]

  new.parent.1 <- ED[new.child.1.row, 2]
  new.parent.1.row <- which(ED[,1] == new.parent.1)

  new.child.2.row <- branch.rows[min(which(cumsum(branch.edge.length) >= new.node.locs[2]))]
  new.child.2 <- ED[new.child.2.row, 1]

  new.parent.2 <- ED[new.child.2.row, 2]
  new.parent.2.row <- which(ED[,1] == new.parent.2)

  #Proposed demes; cannot propose the same colour as the current edge
  proposal.deme.1 <- sample.vector((1:n.deme)[-ED[new.child.1.row, 5]], 1)
  proposal.deme.2 <- sample.vector((1:n.deme)[-ED[new.child.2.row, 5]], 1)

  if (new.child.1 == new.child.2){
    #Special case - selected locations fall on same segment of branch
    if (proposal.deme.1 == proposal.deme.2){
      new.nodes <- max(ED[,1]) + c(1,2)
      which.child <- which(ED[new.parent.1.row, 3:4] == new.child.1)
      ED[new.parent.1.row, 2 + which.child] <- new.nodes[1]
      ED[new.child.1.row, 2] <- new.nodes[2]

      new.node.ages <- ED[which(ED[,1] == branch.leaf), 6] - new.node.locs

      ED <- rbind(ED,
                  c(new.nodes[1], new.parent.1, new.nodes[2], NA, ED[new.child.1.row, 5], new.node.ages[1]),
                  c(new.nodes[2], new.nodes[1], new.child.1, NA, proposal.deme.1, new.node.ages[2]) )

      M.branch <- length(branch.nodes) - 2
      prop.ratio <- branch.length^2 / ((M.branch + 1) * (M.branch + 2))
      return(list(ED = ED, prop.ratio = prop.ratio))

    } else{
      #REJECT in likelihood due to incompatible deme proposal
      return(list(ED = ED, prop.ratio = 0))
    }
  }

  if ((new.child.1 == new.parent.2) && (proposal.deme.1 == proposal.deme.2)){
    #Special case - left-moving deme update and right-moving deme update meet at a migration node and make incompatible deme labelling
    return(list(ED = ED, prop.ratio = 0))
  }

  if ((ED[ED[new.child.1.row, 3],6] == proposal.deme.1)  #Deme below right-moving deme update
      || (ED[new.parent.2.row, 6]== proposal.deme.2)){ #Deme above left-moving deme update
    #REJECT in likelihood due to incompatible deme proposal
    return(list(ED = ED, prop.ratio = 0))
  }

  new.nodes <- max(ED[,1]) + c(1,2)

  which.child <- which(ED[new.parent.1.row, 3:4] == new.child.1)
  ED[new.parent.1.row, 2 + which.child] <- new.nodes[1]
  old.deme <- ED[new.child.1.row, 5]
  ED[new.child.1.row, c(2,5)] <- c(new.nodes[1], proposal.deme.1)

  which.child <- which(ED[new.parent.2.row, 3:4] == new.child.2)
  ED[new.parent.2.row, 2 + which.child] <- new.nodes[2]
  ED[new.child.2.row, 2] <- new.nodes[2]

  new.node.ages <- ED[which(ED[,1] == branch.leaf), 6] - new.node.locs

  ED <- rbind(ED,
              c(new.nodes[1], new.parent.1, new.child.1, NA, old.deme, new.node.ages[1]),
              c(new.nodes[2], new.parent.2, new.child.2, NA, proposal.deme.2, new.node.ages[2])
              )

  M.branch <- length(branch.nodes) - 2
  prop.ratio <- branch.length^2 / ((M.branch + 1) * (M.branch + 2))
  return(list(ED = ED, prop.ratio = prop.ratio))
}

## New migration pair death
ed.mig.pair.death.2 <- function(ED, n.deme){
  root.node <- ED[is.na(ED[,2]), 1]
  migration.nodes <- ED[ED[ is.na(ED[,4]) & (!is.na(ED[,3])) ,1]]

  #Sample branch in genealogy proportional to length
  edge.length <- numeric(dim(ED)[1])
  non.root.nodes <- ED[ED[,1] != root.node, 1]
  non.root.nodes <- non.root.nodes[!is.na(non.root.nodes)]
  for (j in non.root.nodes){
    node.row <- which(ED[,1] == j)
    parent.row <- which(ED[,1] == ED[node.row, 2])
    edge.length[node.row] <- ED[node.row, 6] - ED[parent.row, 6]
  }
  tree.length <- sum(edge.length)  #Total tree length

  genealogy.loc <- runif(1, 0, tree.length)
  #Two nodes above and below genealogy.loc
  below.row <- min(which(cumsum(edge.length) >= genealogy.loc))
  below.node <- ED[below.row, 1]

  above.node <- ED[below.row, 2]
  above.row <- which(ED[,1] == above.node)

  branch.nodes <- c(below.node, above.node)
  branch.rows <- c(below.row, above.row)

  active.node <- above.node
  active.row <- above.row
  while (active.node %in% migration.nodes){
    active.node <- ED[active.row, 2]
    active.row <- which(ED[,1] == active.node)
    branch.nodes <- c(branch.nodes, active.node)
    branch.rows <- c(branch.rows, active.row)
  }
  branch.root <- active.node  #Root end of genealogy branch
  active.node <- below.node
  active.row <- below.row
  while (active.node %in% migration.nodes){
    active.node <- ED[active.row, 3]
    active.row <- which(ED[,1] == active.node)
    branch.nodes <- c(branch.nodes, active.node)
    branch.rows <- c(branch.rows, active.row)
  }
  branch.leaf <- active.node  #Leaf end of genealogy branch
  branch.nodes <- branch.nodes[order(ED[branch.rows, 6], decreasing = TRUE)]  #Ordering branch nodes/rows (necessary for the death move)
  branch.rows <- branch.rows[order(ED[branch.rows, 6], decreasing = TRUE)]

  branch.length <- ED[which(ED[,1] == branch.leaf), 6] - ED[which(ED[,1] == branch.root), 6]

  M.branch <- length(branch.nodes) - 2  #Number of migration events on the branch

  if (M.branch < 2){
    #REJECT as insufficient number of migrations on branch
    return(list(ED = ED, prop.ratio = 0))
  }

  node.sample <- sample(branch.nodes[! branch.nodes %in% c(branch.leaf, branch.root)], 2)
  sample.indices <- which(branch.nodes %in% node.sample)
  sample.indices <- sort(sample.indices, decreasing = TRUE) #Sort sample indices so that sample.indices[1] is the sampled node furthest from the root

  upper.node <- branch.rows[sample.indices[1]] #Left-most sampled node
  upper.parent <- ED[which(ED[,1] == upper.node), 2] #Node above left-most sampled node
  upper.child <- ED[which(ED[,1] == upper.node), 3] #Node below left-most sampled node

  lower.node <- branch.rows[sample.indices[2]]
  lower.parent <- ED[which(ED[,1] == lower.node), 2] #Node above the right-most sampled node
  lower.child <- ED[which(ED[,1] == lower.node), 3] #Node below right-most sampled node

  if (upper.child == lower.node){
    #Special case - sampled nodes are adjacent
    if (ED[upper.node, 5] != ED[lower.child, 5]){
      #REJECT due to inconsistent deme labelling
      return(list(ED = ED, prop.ratio = 0))
    }

    ED[which(ED[,1] == lower.child), 2] <- upper.parent

    parent.row <- which(ED[,1] == upper.parent)
    which.child <- which(ED[parent.row, 3:4] == upper.node)
    ED[parent.row, 2 + which.child] <- lower.child

    ED <- ED[! ED[,1] %in% node.sample, ]

    prop.ratio <- M.branch * (M.branch - 1) / branch.length^2
    return(list(ED = ED, prop.ratio = prop.ratio))
  }

  upper.parent.row <- which(ED[,1] == upper.parent)
  which.child <- which(ED[upper.parent.row, 3:4] == upper.node)
  ED[upper.parent.row, 2 + which.child] <- upper.child

  upper.child.row <- which(ED[,1] == upper.child)
  ED[upper.child.row, 2] <- upper.parent
  ED[upper.child.row, 5] <- ED[which(ED[,1] == upper.node), 5] #Pull deme downwards

  ED[which(ED[,1] == lower.parent), 3] <- lower.child  #lower.parent must be migration node
  ED[which(ED[,1] == lower.child), 2] <- lower.parent

  ED <- ED[!ED[,1] %in% node.sample, ]

  prop.ratio <- M.branch * (M.branch - 1) / branch.length^2
  return(list(ED = ED, prop.ratio = prop.ratio))
}
