mig.rate.update.2 <- function(ED, migration.matrix, n.deme = NA, shape = 1, rate = 10){

  m <- ed.node.count(ED, n.deme)$m

  observed.demes <- unique(ED[,5])
  if (is.na(n.deme)){
    n.deme <- max(observed.demes)
  }
  proposal.matrix <- matrix(0, n.deme, n.deme)

  root.node <- ED[is.na(ED[,2]), 1]
  root.row <- which(ED[,1] == root.node)

  deme.length <- rep(NA, n.deme)

  for (d in 1 : n.deme){
    rows.in.deme <- which(ED[, 5] == d)
    rows.in.deme <- rows.in.deme[rows.in.deme != root.row]

    if (length(rows.in.deme) > 0){
      parents <- ED[rows.in.deme, 2]
      parents.rows <- numeric(length(parents))
      for (i in 1 : length(parents)){
        parents.rows[i] <- which(ED[,1] == parents[i])
      }
      deme.length[d] <- sum(ED[rows.in.deme, 6] - ED[parents.rows, 6])
    } else{
      deme.length[d] <- 0
    }
  }

  for (i in 1:n.deme){
    for (j in (1:n.deme)[-i]){
      proposal.matrix[i,j] <- rgamma(1, shape = shape + m[i,j], scale = 1 / (rate + deme.length[i]))
    }
  }
  return(proposal.matrix)
}


eff.pop.update.2 <- function(ED, effective.population, n.deme, shape = 1, rate = 1){
  c <- ed.node.count(ED, n.deme)$c

  observed.demes <- unique(ED[,5])
  if (is.na(n.deme)){
    n.deme <- length(observed.demes)
  }
  proposal.eff.pop <- rep(0, n.deme)

  event.times <- unique(sort(ED[,6]))
  time.increments <- diff(event.times)
  check.times <- event.times[1:(length(event.times) - 1)] + time.increments/2
  if (any(diff(check.times) == 0)){
    problem.indices <- which(diff(check.times) == 0)
    for (i in 1 : length(problem.indices)){
      check.times <- check.times[-problem.indices[i]]
      time.increments <- time.increments[-problem.indices[i]]
    }
  }

  root.node <- ED[is.na(ED[,2]), 1]
  coalescence.nodes <- ED[!is.na(ED[,4]),1]
  migration.nodes <- ED[ is.na(ED[,4]) & (!is.na(ED[,3])) ,1]
  leaf.nodes <- ED[(is.na(ED[,3])) & (is.na(ED[,4])), 1]

  root.row <- which(ED[,1] == root.node)

  counted <- rep(0, dim(ED)[1])
  counted[root.row] <- 1

  k <- matrix(0, length(check.times), n.deme)
  root.child <- ED[root.row, 3]
  root.child.row <- which(ED[,1] == root.child)

  k[1, ED[root.child.row, 5]] <- 2

  for (i in 2 : length(check.times)){
    current.rows <- which((ED[,6] <= check.times[i]) & (counted == 0))
    counted[current.rows] <- 1
    k[i,] <- k[i-1,]

    if (length(current.rows) == 0){
      if (i == length(check.times)){
        current.rows <- which(counted == 0)
      }
    } else{
      if (length(current.rows) > 1){
        #Add multiple leaves
        for (j in current.rows){
          k[i, ED[j, 5]] <- k[i, ED[j, 5]] - 1
        }
      } else{
        if (ED[current.rows,1] %in% migration.nodes){  #Migration event
          k[i, ED[current.rows, 5]] <- k[i, ED[current.rows, 5]] - 1
          current.child <- ED[current.rows, 3]
          current.child.row <- which(ED[,1] == current.child)
          k[i, ED[current.child.row, 5]] <- k[i, ED[current.child.row, 5]] + 1
        } else if (ED[current.rows,1] %in% coalescence.nodes){
          k[i, ED[current.rows, 5]] <- k[i, ED[current.rows, 5]] + 1
        } else if (ED[current.rows,1] %in% leaf.nodes){ #Single leaf added
          k[i, ED[current.rows, 5]] <- k[i, ED[current.rows, 5]] - 1
        }
      }
    }
  }
  rate.constants <- t(k * (k-1) / 2) %*% time.increments

  for (i in 1:n.deme){
    proposal.eff.pop[i] <- 1/rgamma(1, shape = shape + c[i], scale = 1 / (rate + rate.constants[i]))  #Proposals are inverse-gamma
  }
  return(proposal.eff.pop)
}
