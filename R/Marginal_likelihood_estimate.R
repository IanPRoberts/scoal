#' Marginal likelihood estimation
#'
#' Estimates the marginal likelihood under the structured coalescent using either the harmonic mean, thermodynamic integration or stepping stones sampling
#'
#' @param N MCMC iterations to use in marginal likelihood estimation
#' @param N0 MCMC iterations of burn-in from the posterior distribution
#' @param ED Initial migration history in ED format (9 columns)
#' @param coal_rate Initial coalescent rates
#' @param bit_mig_mat Initial backwards-in-time migration rates
#' @param st_radius Initial radius for migration history proposals
#' @param adaptive Logical value indicating whether to adapt the radius during the burn-in phase (default TRUE)
#' @param cr_mode Coalescent rates prior mode
#' @param cr_var Coalescent rates prior variance
#' @param mm_mode Backward-in-time migration rates prior mode
#' @param mm_var Backward-in-time migration rates prior variance
#' @param adaptation_rate Adaptation rate (unused if adaptive = FALSE)
#' @param target_accept_rate Target acceptance rate for migration history updates
#' @param n_stones Number of stepping stones to use (unused in harmonic mean or thermodynamic integration estimates)
#'
#' @export

structured_coalescent_marginal_likelihood <- function(N, N0, method,
                                                      ED, coal_rate, bit_mig_mat,
                                                      st_radius, adaptive=TRUE,
                                                      cr_mode, cr_var,
                                                      mm_mode, mm_var,
                                                      adaptation_rate=0.6, target_accept_rate=0.234,
                                                      n_stones=100,
                                                      output_file=NULL){
  #Select correct setup depending on input method
  if (method %in% c('harmonic', 'Harmonic', 'h', 'H')){
    estimate_method <- 'Harmonic mean'
    batch_size <- N
    n_batches <- 1
    betas <- as.vector(1)
  } else if (method %in% c('thermodynamic', 'Thermodynamic', 't', 'T')) {
    estimate_method <- 'Thermodynamic integration'
    batch_size <- 1
    n_batches <- N
    betas <- N:0/N
  } else if (method %in% c('stepping_stones', 'Stepping_stones', 's', 'S')){
    estimate_method <- 'Stepping stones sampling'
    n_batches <- n_stones + 1
    alpha <- 0.25
    batch_size <- ceiling(N/(n_batches+1))
    betas <- (0:n_batches/n_batches)^(1/alpha)
  } else {
    stop('Invalid estimation method entered')
  }

  if (!is.null(output_file)){
    cat('Estimate method: ', estimate_method, '\n\n', file=output_file)
    cat('beta\tlog likelihood\n', file=output_file, append=TRUE)
  }

  #Gibbs move for power posterior
  power_posterior_gibbs_moves <- function(ED, ED_NI, coal_rate, bit_mig_mat, power){
    ED_NC <- NodeCountC(ED, n_deme, ED_NI)
    ED_DD <- DemeDecompC(ED, n_deme, ED_NI)

    rate_consts <- colSums((ED_DD$k * (ED_DD$k-1) / 2) * ED_DD$time.increments)
    deme_length <- colSums(ED_DD$k*ED_DD$time.increments)

    cr_post_shape <- cr_shape + power * ED_NC$c
    cr_post_rate <- cr_rate + power * rate_consts
    cr <- rgamma(n_deme, cr_post_shape, cr_post_rate)

    mm_post_shape <- mm_shape + power * ED_NC$m
    mm_post_rate <- mm_rate + power * deme_length
    bmm <- matrix(rgamma(n_deme^2, mm_post_shape, mm_post_rate), n_deme, n_deme)
    diag(bmm) <- 0

    return(list(coal_rate=cr, bit_mig_mat=bmm))
  }

  n_deme <- nrow(bit_mig_mat)

  # Convert priors to shape-rate parameterisation from mode-variance
  cr_rate <- (cr_mode + sqrt(cr_mode^2 + 4 * cr_var))/(2 * cr_var)
  cr_shape <- 1 + cr_mode * cr_rate
  mm_rate <- (mm_mode + sqrt(mm_mode^2 + 4 * mm_var))/(2 * mm_var)
  mm_shape <- 1 + mm_mode * mm_rate

  # Forward-in-time rates and eigen decomposition
  fit_mig_mat <- FitMigMatC(bit_mig_mat, coal_rate)
  fit_rates <- fit_mig_mat
  diag(fit_rates) <- - rowSums(fit_mig_mat)

  eigen_decomp <- eigen(fit_rates)
  eigen_vals <- eigen_decomp$values
  eigen_vecs <- eigen_decomp$vectors
  inverse_vecs <- solve(eigen_vecs)

  # Initial likelihoods and priors
  ED_NI <- NodeIndicesC(ED)
  ED_SC <- SC_like_C(ED, coal_rate, bit_mig_mat, ED_NI)

  log_likelihoods <- matrix(NA, batch_size, n_batches)

  cat('\n', estimate_method, '\n')

  pb <- txtProgressBar(0, batch_size * n_batches, style = 3)

  ##########################
  # Burn in from posterior #
  ##########################

  for (iter_id in 1 : N0){
    setTxtProgressBar(pb, iter_id)

    #Gibbs move updates
    power_post_update <- power_posterior_gibbs_moves(ED, ED_NI, coal_rate, bit_mig_mat, 1) #Update from posterior (i.e. power=1)
    coal_rate <- power_post_update$coal_rate
    bit_mig_mat <- power_post_update$bit_mig_mat

    ED_SC <- SC_like_C(ED, coal_rate, bit_mig_mat, ED_NI) #Update (log-)likelihood following parameter updates

    subtree <- st_centre_dist(ED, st_radius, ED_NI)
    proposal <- local_DTA_subtree_proposal(subtree$EED, subtree$st_labels, fit_rates,
                                           eigen_decomp = eigen_decomp, inverse_vecs = inverse_vecs)

    prop_NI <- NodeIndicesC(proposal$proposal)
    prop_SC <- SC_like_C(proposal$proposal, coal_rate, bit_mig_mat, prop_NI)

    log_AR <- min(0, prop_SC - ED_SC + #log difference in power posteriors
                    local_DTA_likelihood(subtree$st_labels, coal_rate, bit_mig_mat)$log.likelihood - #log(q(H | H'))
                    proposal$prop_prob) #log(q(H' | H))

    if (log(runif(1)) < log_AR){
      # Accept (No need to save ED_SC as recomputed at next iteration)
      ED <- proposal$proposal
      ED_NI <- prop_NI
    }

    if (adaptive){
      # Update proposal radius
      st_radius <- exp(log(st_radius) + iter_id^(-adaptation_rate) * (exp(log_AR) - target_accept_rate))
    }
  }

  ##################################
  # Marginal likelihood estimation #
  ##################################

  for (batch_id in 1 : n_batches){
    for (iter_id in 1 : batch_size){
      setTxtProgressBar(pb, (batch_id - 1) * batch_size + iter_id + N0)
      #Gibbs move updates
      power_post_update <- power_posterior_gibbs_moves(ED, ED_NI, coal_rate, bit_mig_mat, betas[batch_id])
      coal_rate <- power_post_update$coal_rate
      bit_mig_mat <- power_post_update$bit_mig_mat

      ED_SC <- SC_like_C(ED, coal_rate, bit_mig_mat, ED_NI) #Update (log-)likelihood following parameter updates

      subtree <- st_centre_dist(ED, st_radius, ED_NI)
      proposal <- local_DTA_subtree_proposal(subtree$EED, subtree$st_labels, fit_rates,
                                             eigen_decomp = eigen_decomp, inverse_vecs = inverse_vecs)

      prop_NI <- NodeIndicesC(proposal$proposal)
      prop_SC <- SC_like_C(proposal$proposal, coal_rate, bit_mig_mat, prop_NI)

      log_AR <- min(0, betas[batch_id] * (prop_SC - ED_SC) + #log difference in power posteriors
                      local_DTA_likelihood(subtree$st_labels, coal_rate, bit_mig_mat)$log.likelihood - #log(q(H | H'))
                      proposal$prop_prob) #log(q(H' | H))

      if (log(runif(1)) < log_AR){
        # Accept
        ED <- proposal$proposal
        ED_NI <- prop_NI
        ED_SC <- prop_SC
      }
      log_likelihoods[iter_id, batch_id] <- ED_SC

      if (!is.null(output_file)){
        cat(betas[batch_id], '\t', ED_SC, '\n', file=output_file, append = TRUE)
      }
    }
  }

  close(pb)

  if (estimate_method == 'Harmonic mean'){
    #Factor out minimum likelihood for numerical stability
    min_likelihood <- min(log_likelihoods)
    log_likelihoods <- log_likelihoods - min_likelihood
    marginal_likelihood <- log(batch_size * n_batches / sum(exp(-log_likelihoods))) + min_likelihood
  } else if (estimate_method == "Thermodynamic integration") {
    #Potential function given by log likelihood (already stored)
    log_likelihoods <- as.vector(log_likelihoods)
    marginal_likelihood <- (2 * sum(log_likelihoods) - log_likelihoods[1] - log_likelihoods[n_batches])/ (2 * n_batches)
  } else if (estimate_method == 'Stepping stones sampling'){
    #Estimate each ratio of normalising constants and take product to give (unlogged) marginal likelihood estimate
    ratio_estimates <- min_likelihood <- numeric(n_batches)
    for (batch_id in 1 : n_batches){
      #Factor out minimum likelihood for numerical stability - can use different minimum for each annealing level (beta value)
      min_likelihood[batch_id] <- min(log_likelihoods[,batch_id])
      log_likelihoods[,batch_id] <- log_likelihoods[,batch_id] - min_likelihood[batch_id]
      ratio_estimates[batch_id] <- log(mean(exp((betas[batch_id + 1] - betas[batch_id]) * log_likelihoods[,batch_id]))) + (betas[batch_id + 1] - betas[batch_id]) * min_likelihood[batch_id]
    }
    marginal_likelihood <- sum(ratio_estimates)
  }

  if (!is.null(output_file)){
    cat('\nlog(marginal likelihood) = ', marginal_likelihood, '\n', file=output_file, append = TRUE)
  }

  return(marginal_likelihood)
}

#' Partial marginal likelihood computation
#'
#' Computes the partial marginal likelihood obtained by integrating coalescent and migration rates out of the structured coalescent posterior
#'
#' @param ED Migration history in ED format (9 columns)
#' @param cr_mode Coalescent rates prior mode
#' @param cr_var Coalescent rates prior variance
#' @param mm_mode Backward-in-time migration rates prior mode
#' @param mm_var Backward-in-time migration rates prior variance
#' @param n_deme Number of demes
#'
#' @export

partial_marginal_likelihood <- function(ED, cr_mode, cr_var, mm_mode, mm_var, n_deme){
  # Convert priors to shape-rate parameterisation from mode-variance
  cr_rate <- (cr_mode + sqrt(cr_mode^2 + 4 * cr_var))/(2 * cr_var) #beta
  cr_shape <- 1 + cr_mode * cr_rate #alpha
  mm_rate <- (mm_mode + sqrt(mm_mode^2 + 4 * mm_var))/(2 * mm_var) #beta
  mm_shape <- 1 + mm_mode * mm_rate #alpha

  ED_NI <- NodeIndicesC(ED)
  ED_DD <- DemeDecompC(ED, n_deme, ED_NI)
  ED_NC <- NodeCountC(ED, n_deme, ED_NI)

  rate_consts <- colSums((ED_DD$k * (ED_DD$k-1) / 2) * ED_DD$time.increments)
  deme_length <- colSums(ED_DD$k*ED_DD$time.increments)

  log_pml <- cr_shape * n_deme * log(cr_rate) + mm_shape * n_deme * (n_deme - 1) * log(mm_rate) -
    n_deme * log(gamma(cr_shape)) - n_deme * (n_deme - 1) * log(gamma(mm_shape)) +
    sum(log(cr_shape + ED_NC$c) - (cr_shape + ED_NC$c) * log(cr_rate + rate_consts)) +
    sum(log(mm_shape + ED_NC$m) - (mm_shape + ED_NC$m) * log(mm_rate + deme_length)) -
    n_deme * log(mm_shape) + mm_shape * sum(log(mm_rate + deme_length)) #Correction for not removing diagonal entries on mm line

  return(list(likelihood=exp(log_pml), log.likelihood=log_pml))
}
