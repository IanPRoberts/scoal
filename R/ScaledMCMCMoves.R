##### Scaled structured coalescent

rel.eff.pop.prop <- function(ED, timescale.factor, n.deme, node.indices, shape = 1, rate = 1){
  c <- NodeCountC(ED, n.deme, node.indices)$c
  DemeDecomp <- DemeDecompC(ED, n.deme, node.indices)
  k <- DemeDecomp$k
  rate.constants <- t(k * (k-1) / 2) %*% DemeDecomp$time.increments

  proposal.eff.pop <- numeric(n.deme)

  for (i in 1:n.deme){
    proposal.eff.pop[i] <- 1/rgamma(1, shape = shape + c[i] / timescale.factor, rate = (rate + rate.constants[i]))  #Proposals are inverse-gamma
  }
  return(proposal.eff.pop)
}

rel.mig.prop <- function(ED, timescale.factor, n.deme, node.indices, shape = 1, rate = 10){
  m <- NodeCountC(ED, n.deme, node.indices)$m
  DemeDecomp <- DemeDecompC(ED, n.deme, node.indices)
  k <- DemeDecomp$k
  deme.length <- as.vector(t(k) %*% DemeDecomp$time.increments)

  proposal.matrix <- matrix(0, n.deme, n.deme)
  for (i in 1:n.deme){
    for (j in (1:n.deme)[-i]){
      proposal.matrix[i,j] <- rgamma(1, shape = shape + m[i,j], rate = rate + timescale.factor * deme.length[i])
    }
  }
  return(proposal.matrix)
}

timescale.factor.prop <- function(timescale.factor, rel.eff.pop, rel.mig.mat, ED, n.deme, node.indices,
                                  prior.shape = 0.001, prior.rate = 0.001, proposal.sd = 0.1){
  proposal <- abs(rnorm(1, timescale.factor, proposal.sd))

  node_count <- NodeCountC(ED, n.deme, node.indices)
  c <- node_count$c
  m <- node_count$m
  M <- sum(as.vector(m)) #Number of migrations
  n <- (dim(ED)[1] - M + 1)/2 #Number of leaves

  deme_decomp <- DemeDecompC(ED, n.deme, node.indices)
  k <- deme_decomp$k
  time_increments <- deme_decomp$time.increments
  coal_consts <- colSums(k * (k-1) * time_increments) / 2
  deme_lengths <- colSums(k * time_increments)

  log_accept_prob <- min(0, (1/timescale.factor - 1/proposal) * sum(coal_consts * 1/rel.eff.pop) +
                           (M - n + prior.shape) * (log(proposal) - log(timescale.factor)) +
                           (timescale.factor - proposal) * (sum(deme_lengths * rowSums(rel.mig.mat)) + prior.rate))
  return(list(timescale.factor = proposal, prop.ratio = exp(log_accept_prob), log.prop.ratio = log_accept_prob))
}
