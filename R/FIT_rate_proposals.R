##### Forward-in-time Gibbs moves

fit_scaled_mm_gibbs <- function(ED, coal_rate, time_scale, n_deme, node_indices, prior_shape = matrix(1, n_deme, n_deme), prior_rate = matrix(1, n_deme, n_deme)){
  m <- NodeCountC(ED, n_deme, node_indices)$m
  deme_decomp <- DemeDecompC(ED, n_deme, node_indices)
  deme_length <- colSums(deme_decomp$k*deme_decomp$time.increments)

  proposal <- matrix(rgamma(n_deme^2, prior_shape + t(m), prior_rate + time_scale * coal_rate %*% t(deme_length / coal_rate)), n_deme, n_deme)
  diag(proposal) <- 0
  return(proposal)
}

fit_ts_gibbs <- function(ED, coal_rate, mig_mat, n_deme, node_indices, prior_shape = 0.001, prior_rate = 0.001){
  node_count <- NodeCountC(ED, n_deme, node_indices)

  deme_decomp <- DemeDecompC(ED, n_deme, node_indices)
  k <- deme_decomp$k
  time_increments <- deme_decomp$time.increments
  rate_consts <- colSums((k * (k-1) / 2) * deme_decomp$time.increments)

  proposal <- rgamma(1, shape = prior_shape + sum(node_count$m) + sum(node_count$c), prior_rate + sum(rate_consts * coal_rate) + t(coal_rate) %*% mig_mat %*% diag(1/coal_rate) %*% t(k) %*% time_increments)
  return(proposal)
}

##### Forward-in-time Metropolis-Hastings for coal rate
fit_scaled_cr <- function(ED, coal_rate, mig_mat, time_scale, n_deme, node_indices, prior_shape = rep(1, n_deme), prior_rate = rep(1,n_deme), proposal_sd = time_scale){
  NC <- NodeCountC(ED, n_deme, node_indices)
  DD <- DemeDecompC(ED, n_deme, node_indices)

  rate_consts <- t(choose(DD$k, 2)) %*% DD$time.increments
  deme_lengths <- t(DD$k) %*% DD$time.increments

  t1 <- (t(NC$m) - NC$m) %*% rep(1,n_deme) + prior_shape +NC$c - 1
  t2 <- prior_rate + time_scale * rate_consts + time_scale * mig_mat %*% diag(1/coal_rate) %*% deme_lengths
  t3 <- time_scale * diag(as.vector(coal_rate %*% mig_mat)) %*% deme_lengths

  proposal <- abs(rnorm(n_deme, coal_rate, proposal_sd))
  log_ar <- sapply(t1 * (log(proposal) - log(coal_rate)) - t2 * (proposal - coal_rate) - t3 * (1/proposal - 1/coal_rate), function(x) min(0,x))
  accepted <- (log(runif(n_deme)) < log_ar)

  proposal[!accepted] <- coal_rate[!accepted]

  return(proposal)
}
