eff.pop.mh <- function(effective.pop, ED, n.deme = length(effective.pop), node.indices,
                       prior.shape = 0.001, prior.rate = 0.001, proposal.sd = sqrt(0.2)){
  proposal <- rnorm(n.deme, effective.pop, proposal.sd)

  node_count <- NodeCountC(ED, n.deme, node.indices)
  c <- node_count$c

  deme_decomp <- DemeDecompC(ED, n.deme, node.indices)
  k <- deme_decomp$k
  time_increments <- deme_decomp$time.increments
  coal_consts <- colSums(k * (k-1) * time_increments) / 2

  log_accept_prob <- min(0, sum(-(prior.shape + c + 1) * (log(proposal) - log(effective.pop)) - (prior.rate + coal_consts) * (1/proposal - 1 / effective.pop)))

  if (log(runif(1)) < log_accept_prob){
    return(proposal)
  } else{
    return(effective.pop)
  }
}
