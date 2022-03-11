dta.mig.rate.update <- function(ED, migration.matrix, n.deme = NA, node.indices, shape = 1, rate = 10){

  m <- ed.node.count(ED, n.deme, node.indices)$m

  observed.demes <- unique(ED[,5])
  if (is.na(n.deme)){
    n.deme <- dim(migration.matrix)[1]
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
      parents.rows <- node.indices[parents]
      deme.length[d] <- sum(ED[rows.in.deme, 6] - ED[parents.rows, 6])
    } else{
      deme.length[d] <- 0
    }
  }

  for (i in 1:n.deme){
    for (j in (1:n.deme)[-i]){
      proposal.matrix[i,j] <- rgamma(1, shape = shape + m[j,i], scale = 1 / (rate + deme.length[i]))
    }
  }
  return(proposal.matrix)
}
