##### Scaled structured coalescent

rel.eff.pop.prop <- function(ED, timescale.factor, n.deme, node.indices, shape = 1, rate = 1){
  c <- NodeCountC(ED, n.deme, node.indices)$c
  DemeDecomp <- DemeDecompC(ED, n.deme, node.indices)
  k <- DemeDecomp$k
  rate.constants <- t(k * (k-1) / 2) %*% DemeDecomp$time.increments

  proposal.eff.pop <- numeric(n.deme)

  for (i in 1:n.deme){
    proposal.eff.pop[i] <- 1/rgamma(1, shape = shape + c[i], rate = rate + rate.constants[i] / timescale.factor)  #Proposals are inverse-gamma
  }
  return(proposal.eff.pop)
}

rel.ep.mh <- function(ED, rel.eff.pop, time.scale, n.deme, node.indices, prior.shape = 1, prior.rate = 1, proposal.sd = sqrt(10 * time.scale)){
  c <- NodeCountC(ED, n.deme, node.indices)$c
  DemeDecomp <- DemeDecompC(ED, n.deme, node.indices)
  k <- DemeDecomp$k
  rate.constants <- t(k * (k-1) / 2) %*% DemeDecomp$time.increments

  proposal <- abs(rnorm(n.deme, rel.eff.pop, proposal.sd))
  log_accept_prob <- min(0, - sum((1/proposal - 1/rel.eff.pop) * (prior.rate + rate.constants/time.scale)) + sum((c + prior.shape + 1) * (log(rel.eff.pop) - log(proposal))))

  return(list(rel.eff.pop = proposal, prop.ratio = exp(log_accept_prob), log.prop.ratio = log_accept_prob))
}


rel.mig.prop <- function(ED, timescale.factor, n.deme, node.indices, shape = 1, rate = 10){
  m <- NodeCountC(ED, n.deme, node.indices)$m
  DemeDecomp <- DemeDecompC(ED, n.deme, node.indices)
  k <- DemeDecomp$k
  deme.length <- as.vector(t(k) %*% DemeDecomp$time.increments)

  proposal.matrix <- matrix(0, n.deme, n.deme)
  for (i in 1:n.deme){
    for (j in (1:n.deme)[-i]){
      proposal.matrix[i,j] <- rgamma(1, shape = shape + m[i,j], rate = rate + deme.length[i] / timescale.factor)
    }
  }
  return(proposal.matrix)
}

rel.mm.mh <- function(ED, rel.mig.mat, time.scale, n.deme, node.indices, prior.shape = 1, prior.rate = 10, proposal.sd = 10/time.scale){
  m <- NodeCountC(ED, n.deme, node.indices)$m
  DemeDecomp <- DemeDecompC(ED, n.deme, node.indices)
  k <- DemeDecomp$k
  deme.length <- as.vector(t(k) %*% DemeDecomp$time.increments)

  proposal <- matrix(0, n.deme, n.deme)
  for (i in 1:n.deme){
    for (j in (1:n.deme)[-i]){
      proposal[i,j] <- abs(rnorm(1, rel.mig.mat[i,j], proposal.sd))
    }
  }

  diag(proposal) <- diag(rel.mig.mat) <- 1

  log_accept_prob <- min(0, - sum(deme.length * (rowSums(proposal) - rowSums(rel.mig.mat)))/time.scale - prior.rate * sum(proposal - rel.mig.mat) + sum((m+prior.shape -1) * (log(proposal) - log(rel.mig.mat))))
  diag(proposal) <- 0
  return(list(rel.mig.mat = proposal, prop.ratio = exp(log_accept_prob), log.prop.ratio = log_accept_prob))
}
# timescale.factor.prop <- function(timescale.factor, rel.eff.pop, rel.mig.mat, ED, n.deme, node.indices,
#                                   prior.shape = 0.001, prior.rate = 0.001, proposal.sd = 0.1){
#   proposal <- abs(rnorm(1, timescale.factor, proposal.sd))
#
#   node_count <- NodeCountC(ED, n.deme, node.indices)
#   c <- node_count$c
#   m <- node_count$m
#   M <- sum(as.vector(m)) #Number of migrations
#   n <- (dim(ED)[1] - M + 1)/2 #Number of leaves
#
#   deme_decomp <- DemeDecompC(ED, n.deme, node.indices)
#   k <- deme_decomp$k
#   time_increments <- deme_decomp$time.increments
#   coal_consts <- colSums(k * (k-1) * time_increments) / 2
#   deme_lengths <- colSums(k * time_increments)
#
#   log_accept_prob <- min(0, (1/timescale.factor - 1/proposal) * sum(coal_consts * 1/rel.eff.pop) +
#                            (M - n + prior.shape) * (log(proposal) - log(timescale.factor)) +
#                            (timescale.factor - proposal) * (sum(deme_lengths * rowSums(rel.mig.mat)) + prior.rate))
#   return(list(timescale.factor = proposal, prop.ratio = exp(log_accept_prob), log.prop.ratio = log_accept_prob))
# }

timescale.factor.prop <- function(timescale.factor, rel.eff.pop, rel.mig.mat, ED, n.deme, node.indices,
                                  prior.shape = 0.001, prior.rate = 0.001){
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

  proposal <- 1/rgamma(1, shape = prior.shape + M + n - 1, rate = prior.rate + sum(coal_consts / rel.eff.pop) + sum(deme_lengths * rowSums(rel.mig.mat)))

  return(proposal)
}
