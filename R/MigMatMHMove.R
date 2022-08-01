mig.mat.mh <- function(ED, mig.mat, n.deme, node.indices, prior.shape = 1, prior.rate = 10, proposal.sd = 1){
  # m <- NodeCountC(ED, n.deme, node.indices)$m
  # DemeDecomp <- DemeDecompC(ED, n.deme, node.indices)
  # deme.length <- as.vector(t(DemeDecomp$k) %*% DemeDecomp$time.increments)
  #
  # proposal <- matrix(0, n.deme, n.deme)
  # for (i in 1:n.deme){
  #   for (j in (1:n.deme)[-i]){
  #     proposal[i,j] <- abs(rnorm(1, mig.mat[i,j], proposal.sd))
  #   }
  # }
  #
  # diag(proposal) <- diag(mig.mat) <- 1
  #
  # log_accept_prob <- min(0, - sum((deme.length + prior.rate) * (rowSums(proposal) - rowSums(mig.mat))) + sum((m+prior.shape -1) * (log(proposal) - log(mig.mat))))
  # diag(proposal) <- diag(mig.mat) <- 0
  #
  # if (log(runif(1)) < log_accept_prob){
  #   return(proposal)
  # } else{
  #   return(mig.mat)
  # }

  m <- NodeCountC(ED, n.deme, node.indices)$m
  DemeDecomp <- DemeDecompC(ED, n.deme, node.indices)
  deme.length <- as.vector(t(DemeDecomp$k) %*% DemeDecomp$time.increments)

  for (i in 1:n.deme){
    for (j in (1:n.deme)[-i]){
      proposal <- abs(rnorm(1, mig.mat[i,j], proposal.sd))
      log_accept_prob <- min(0, -(deme.length[i] + prior.rate) * (proposal - mig.mat[i,j]) + (m[i,j] + prior.shape - 1) * (log(proposal) - log(mig.mat[i,j])))

      if (log(runif(1)) < log_accept_prob){
        mig.mat[i,j] <- proposal
      }
    }
  }

  return(mig.mat)
}
