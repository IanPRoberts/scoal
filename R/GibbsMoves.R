mig.rate.update <- function(ED, migration.matrix, n.deme = NA, node.indices, shape = 1, rate = 10){

  m <- NodeCountC(ED, n.deme, node.indices)$m #ed.node.count(ED, n.deme, node.indices)$m
  if (is.na(n.deme)){
    n.deme <- max(ED[,5])
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
      parent.rows <- node.indices[parents]
      deme.length[d] <- sum(ED[rows.in.deme, 6] - ED[parent.rows, 6])
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
  c <- NodeCountC(ED, n.deme, node.indices)$c

  if (is.na(n.deme)){
    n.deme <- max(ED[,5])
  }
  proposal.eff.pop <- rep(0, n.deme)

  DemeDecomp <- DemeDecompC(ED, n.deme, node.indices)
  k <- DemeDecomp$k
  rate.constants <- t(k * (k-1) / 2) %*% DemeDecomp$time.increments

  for (i in 1:n.deme){
    proposal.eff.pop[i] <- 1/rgamma(1, shape = shape + c[i], scale = 1 / (rate + rate.constants[i]))  #Proposals are inverse-gamma
  }
  return(proposal.eff.pop)
}
