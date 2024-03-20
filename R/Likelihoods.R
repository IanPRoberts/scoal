#' Kingman Coalescent Likelihood
#'
#' Returns the likelihood and log-likelihood of a (possibly heterochronous)
#' coalescent tree
#'
#' @param phylo object of class \code{phylo}
#' @param coal_rate Coalescent rate of the population (reciprocal of the product of effective population size and generation length)
#'
#' @return A vector (\code{log-likelihood}, \code{likelihood}) giving the log-likelihood and likelihood of \code{phylo}
#'
#' @export

phylo_likelihood <- function(phylo, coal_rate){
  n <- length(phylo$tip.label)

  node.heights <- node.depth.edgelength(phylo)  #Horizontal distance from root
  event.times <- sort(unique(node.heights))  #Times of events, root at time=0
  time.increments <- diff(event.times)

  #Augments phylo$edge to add edge "start" and "end" times
  aug.edge <- matrix(0,nrow = (2 * (n-1)), ncol = 4)
  aug.edge[,c(1,2)] <- phylo$edge

  for (i in 1 : (2 * (n-1))){
    aug.edge[i,c(3,4)] <- c(node.heights[aug.edge[i,c(1,2)]])
  }

  #Calculate number of edges in tree midway between event times
  check.time <- event.times[-length(event.times)] + 0.5 * time.increments
  k <- rep(0,length(time.increments))  #Number of edges in tree between events
  for (i in 1 : length(time.increments)){
    k[i] <- sum( (aug.edge[,3] < check.time[i]) & (aug.edge[,4] > check.time[i]) )
  }

  likelihood <- 0  #Initialise log likelihood at 0
  likelihood <- sum(- (k * ( k-1) / 2 * coal_rate) * time.increments) + (n-1) * log(coal_rate)
  return(c(likelihood, exp(likelihood)))  #Output (log-likelihood, likelihood)
}

#' Structured Coalescent Likelihood
#'
#' Computes the likelihood and log-likelihood of a structured genealogy under the
#' structured coalescent model
#'
#' @param ED Extended data representation of a structured genealogy
#' @param coal_rate Vector of coalescent rates
#' @param bit_mig_mat Matrix of backwards-in-time migration rates
#' @param ED_NI Node indices obtained from NodeIndices()
#'
#' @return A vector (\code{log-likelihood}, \code{likelihood}) giving the log-likelihood and likelihood of \code{phylo}
#'
#' @export

SC_likelihood <- function(ED, coal_rate, bit_mig_mat, ED_NI){
  n_deme <- nrow(bit_mig_mat)
  ED_DD <- DemeDecomp(ED, n_deme, ED_NI)
  ED_NC <- NodeCount(ED, n_deme, ED_NI)

  k <- ED_DD$k
  time_increments <- ED_DD$time.increments
  c <- ED_NC$c
  m <- ED_NC$m

  deme_lengths <- colSums(k * time_increments)
  coal_constants <- colSums(k * (k-1) * time_increments) / ( 2 * lambda)
  mm_row_sum <- rowSums(bit_mig_mat)
  log_mig_mat <- log(bit_mig_mat)
  diag(log_mig_mat) <- 0

  like <- sum(m * log_mig_mat) - sum(c * log(effective.pop)) - sum(coal_constants) - sum(deme_lengths * mm_row_sum)

  return(list(log.likelihood = like, likelihood = exp(like)))
}


#' Discrete Trait Analysis likelihood
#'
#' Returns the DTA likelihood and log-DTA-likelihood of a structured coalescent tree
#'
#' @inheritParams SC_likelihood
#' @param fit_mig_mat Matrix of forwards-in-time migration rates
#'
#' @return A vector (\code{log-likelihood}, \code{likelihood}) giving the log-likelihood and likelihood of \code{phylo}
#'
#' @export

DTA_likelihood <- function(ED, fit_mig_mat, ED_NI){
  root_row <- which(is.na(ED[,2]))
  fmm_rowsums <- rowSums(fit_mig_mat)
  log_likelihood <- 0

  for (i in (1 : dim(ED)[1])[-root_row]){
    node_parent_row <- ED_NI[ED[i,2]]
    time_increment <- ED[i, 6] - ED[node_parent_row, 6]

    parent_deme <- ED[i, 5]
    if (is.na(ED[i,3])){ #Leaf node
      node_deme <- parent_deme
    } else{
      node_deme <- ED[ED_NI[ED[i,3]], 5]
    }

    log_likelihood <- log_likelihood - fmm_rowsums[parent_deme] * time_increment

    if (parent_deme != node_deme){
      log_likelihood <- log_likelihood + log(f_mm[parent_deme, node_deme])
    }
  }
  return(list(log.likelihood = log_likelihood, likelihood = exp(log_likelihood)))
}

#' MultiTypeTree Probability Density
#'
#' Computes the relevant probability density to compute transition kernels for
#' the MultiTypeTree node retype move
#'
#' @inheritParams SC_likelihood
#' @param bit_rates Transition matrix corresponding to the backwards-in-time migration process
#' @param eigen_vals Vector of eigenvalues of \code{bit_rates}
#' @param eigen_vecs Matrix of eigenvectors of \code{bit_rates}
#' @param inverse_vecs Inverse matrix of \code{eigen_vecs}
#'
#' @return List consisting of the log.likelihood and likelihood of the structured genealogy
#'
#' @export

MTT_likelihood <- function(ED, bit_rates, ED_NI = NodeIndices(ED), eigen_vals = eigen(bit_rates)$values, eigen_vecs = eigen(bit_rates)$vectors, inverse_vecs = solve(eigen_vecs)){
  parent_rows <- ED_NI[ED[,2]]
  edge_lengths <- ED[,6] - ED[parent_rows, 6]
  log_like <- 0

  parent_demes <- ED[parent_rows, 5]
  for (i in 1 : nrow(ED)){
    if (!is.na(parent_rows[i])){
      current_deme <- ED[i,5]
      log_like <- log_like + edge_lengths[i] * bit_rates[current_deme, current_deme] #bit_rates[i,i] = rowSums(bit_mig_mat)[i]

      if (parent_demes[i] != current_deme){
        log_like <- log_like + log(bit_rates[current_deme, parent_demes[i]])
      }
    }
  }

  non_mig_nodes <- which((!is.na(ED[,4])) | (is.na(ED[,3])))

  for (row_id in non_mig_nodes){
    parent_row <- ED_NI[ED[row_id, 7]]
    if (!is.na(parent_row)){
      edge_length <- ED[row_id, 6] - ED[parent_row, 6]
      trans_mat <- eigen_vecs %*% diag(exp(edge_length * eigen_vals)) %*% inverse_vecs
      log_like <- log_like - log(trans_mat[ED[row_id, 5], ED[parent_row, 5]])
    }
  }

  return(list(log.likelihood = log_like, likelihood = exp(log_like)))
}
