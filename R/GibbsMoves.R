mig.rate.update <- function(ED, migration.matrix, n.deme = NA, node.indices, shape = 1, rate = 10){

  m <- ed.node.count(ED, n.deme, node.indices)$m

  observed.demes <- unique(ED[,5])
  if (is.na(n.deme)){
    n.deme <- max(observed.demes)
  }
  proposal.matrix <- matrix(0, n.deme, n.deme)

  root.node <- ED[is.na(ED[,2]), 1]
  root.row <- node.indices[root.node]

  deme.length <- rep(NA, n.deme)

  for (d in 1 : n.deme){
    rows.in.deme <- which(ED[, 5] == d)
    rows.in.deme <- rows.in.deme[rows.in.deme != root.row]

    if (length(rows.in.deme) > 0){
      parents <- ED[rows.in.deme, 2]
      parents.rows <- numeric(length(parents))
      for (i in 1 : length(parents)){
        parents.rows[i] <- node.indices[parents[i]]
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


eff.pop.update <- function(ED, effective.population, n.deme, node.indices, shape = 1, rate = 1){
  c <- ed.node.count(ED, n.deme, node.indices)$c

  observed.demes <- unique(ED[,5])
  if (is.na(n.deme)){
    n.deme <- length(observed.demes)
  }
  proposal.eff.pop <- rep(0, n.deme)

  event.times <- unique(sort(ED[,6]))
  time.increments <- diff(event.times)

  k <- deme.decomp(ED, n.deme, node.indices)
  rate.constants <- t(k * (k-1) / 2) %*% time.increments

  for (i in 1:n.deme){
    proposal.eff.pop[i] <- 1/rgamma(1, shape = shape + c[i], scale = 1 / (rate + rate.constants[i]))  #Proposals are inverse-gamma
  }
  return(proposal.eff.pop)
}
