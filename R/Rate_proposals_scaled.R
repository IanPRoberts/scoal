#' Coalescent Rate Gibbs Update
#'
#' Performs Gibbs updates for all parameters in a structured coalescent phylogeny
#'
#'
#' @param ED Extended data representation of a phylogeny including migration history
#' @param time_scale Time scale for the migration matrix
#' @param node_indices Vector identifying row of each node label
#' @param prior_shape Gamma shape parameter for prior
#' @param prior_rate Gamma rate parameter for prior
#'
#' @return Updated parameters
#'
#' @export

scaled_coal_rate_gibbs_update <- function(ED, time_scale, n_deme, node_indices, prior_shape = 1, prior_rate = 1){
  c <- NodeCountC(ED, n_deme, node_indices)$c
  DemeDecomp <- DemeDecompC(ED, n_deme, node_indices)
  k <- DemeDecomp$k
  rate_consts <- colSums((k * (k-1) / 2) * DemeDecomp$time.increments) #t(k * (k-1) / 2) %*% DemeDecomp$time.increments

  proposal <- rgamma(n_deme, prior_shape + c, prior_rate + time_scale * rate_consts)
  return(proposal)
}

#' Migration Rates Matrix Gibbs Update
#'
#' Performs Gibbs updates for all parameters in a structured coalescent phylogeny
#'
#'
#' @param ED Extended data representation of a phylogeny including migration history
#' @param time_scale Time scale for the migration matrix
#' @param n_deme Number of demes
#' @param node_indices Vector identifying row of each node label
#' @param prior_shape Gamma shape parameter for prior
#' @param prior_rate Gamma rate parameter for prior
#'
#' @return Updated parameters
#'
#' @export


scaled_mig_mat_gibbs_update <- function(ED, time_scale, n_deme, node_indices, prior_shape = 1, prior_rate = 1){
  m <- NodeCountC(ED, n_deme, node_indices)$m
  DemeDecomp <- DemeDecompC(ED, n_deme, node_indices)
  k <- DemeDecomp$k
  deme_length <- colSums(k*DemeDecomp$time.increments) #as.vector(t(k) %*% DemeDecomp$time.increments)

  proposal <- matrix(rgamma(n_deme^2, prior_shape + m, prior_rate + deme_length * time_scale), n_deme, n_deme)
  diag(proposal) <- 0
  return(proposal)
}

#' Time Scale Gibbs Update
#'
#' Performs Gibbs updates for all parameters in a structured coalescent phylogeny
#'
#'
#' @param ED Extended data representation of a phylogeny including migration history
#' @param coal_rate Vector of coalescent rates
#' @param mig_mat Matrix of migration rates
#' @param n_deme Number of demes
#' @param node_indices Vector identifying row of each node label
#' @param prior_shape Gamma shape parameter for prior
#' @param prior_rate Gamma rate parameter for prior
#'
#' @return Updated parameters
#'
#' @export

scaled_time_scale_gibbs_update <- function(coal_rate, mig_mat, ED, n_deme, node_indices, prior_shape = 0.001, prior_rate = 0.001){
  node_count <- NodeCountC(ED, n_deme, node_indices)
  c <- node_count$c
  m <- node_count$m

  deme_decomp <- DemeDecompC(ED, n_deme, node_indices)
  k <- deme_decomp$k
  time_increments <- deme_decomp$time.increments
  rate_consts <- colSums((k * (k-1) / 2) * deme_decomp$time.increments)
  deme_length <- colSums(k*deme_decomp$time.increments)

  proposal <- rgamma(1, shape = prior_shape + sum(m) + sum(c), rate = prior_rate + sum(rate_consts * coal_rate) + sum(deme_length * rowSums(mig_mat)))

  return(proposal)
}
