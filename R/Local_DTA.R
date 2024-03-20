#' Local DTA Sampling
#'
#' Updates the migration history associated with a subtree under a conditional DTA model
#'
#' @param ED Extended data representation of phylogenetic tree and initial migration history
#' @param st_labels Reduced extended data structure consisting of the rows of ED corresponding to coalescent nodes in the selected subtree only
#' @param bit_rates Transition matrix of a Markov process corresponding to a process with rates given by the backwards-in-time migration matrix
#' @param ED_NI Vector of row indices corresponding to which row of ED corresponds to each node label
#'
#' @return  List consisting of proposal, a structured phylogeny with updates made to the selected subtree, and prop_prob, the probability of the given update
#'
#' @export

local_DTA <- function(ED, st_labels, fit_rates, ED_NI = NodeIndices(ED),
                      eigen_decomp = eigen(fit_rates), inverse_vecs = solve(eigen_decomp$vectors)){
  n_deme <- nrow(fit_rates)

  #Identify migrations with child and parent inside subtree
  st_rm_rows <- is.na(st_labels[,4]) & (st_labels[,2] %in% st_labels[,1]) & (st_labels[,3] %in% st_labels[,1])
  ED_rm_rows <- ED_NI[st_labels[st_rm_rows, 1]]

  #Remove migrations from inside st_labels
  st_labels <- st_labels[!st_rm_rows,]

  ######################
  # Backwards messages #
  ######################
  #Reorder subtree nodes into decreasing node age
  st_order <- order(st_labels[,6], decreasing = TRUE)
  st_labels <- st_labels[st_order,]
  st_rows <- ED_NI[st_labels[,1]]


  n_nodes <- nrow(st_labels)
  st_root <- st_labels[n_nodes,1] #Label of st_root -> node with least node age
  st_root_row <- ED_NI[st_root]

  #Parent coal nodes within subtree
  st_parent_labels <- ED[st_rows, 7]
  st_parent_labels[st_parent_labels == st_parent_labels[n_nodes]] <- st_root
  st_parent_labels[n_nodes] <- st_root #In case st_root is global root
  st_parent_rows <- ED_NI[st_parent_labels]
  st_edge_lengths <- ED[st_rows, 6] - ED[st_parent_rows, 6]
  st_edge_lengths[is.na(st_edge_lengths)] <- 0


  #Child coal nodes within subtree
  st_children <- matrix(NA, n_nodes, 2) #unname(ED[st_rows, 8:9])
  st_child_rows <- matrix(NA, n_nodes, 2) #matrix(ED_NI[st_children], ncol = 2)

  for (st_id in 1 : n_nodes){
    for (child_id in 1 : 2){
      child_label <- ED[st_rows[st_id], 7 + child_id]
      child_row <- ED_NI[child_label]

      if (!is.na(child_row)){
        while (!(child_label %in% st_labels[,1])){
          child_label <- ED[child_row, 2]
          child_row <- ED_NI[child_label]
        }
      }
      st_children[st_id, child_id] <- child_label
      st_child_rows[st_id, child_id] <- child_row
    }
  }

  st_child_ids <- matrix(NA, n_nodes, 2)
  st_parent_ids <- numeric(n_nodes)
  for (st_id in 1 : n_nodes){
    st_parent_ids[st_parent_labels == st_labels[st_id]] <- st_id
    st_child_ids[st_children == st_labels[st_id]] <- st_id
  }

  #Transition matrices
  trans_mats <- array(0, c(n_deme, n_deme, n_nodes))
  for (st_id in 1:(n_nodes - 1)){ #Don't need to exponentiate to identity matrix for edge length 0 above root
    trans_mats[,,st_id] <- Re(eigen_decomp$vectors %*% diag(exp(st_edge_lengths[st_id] * eigen_decomp$values)) %*% inverse_vecs)
  }

  messages <- array(NA, dim = c(n_nodes, n_nodes, n_deme))

  for (st_id in 1 : (n_nodes - 1)){ #Skip subtree root (st_root corresponds to st_id == n_nodes)
    node_row <- st_rows[st_id]

    if ((is.na(st_labels[st_id, 3])) || st_children[st_id, 1] == ED[node_row, 1]){
      # Subtree leaf node contributes column of transition matrix ending at current deme
      current_deme <- ED[node_row, 5]
      messages[st_id, st_parent_ids[st_id], ] <- trans_mats[, current_deme, st_id]
    } else {
      # Other subtree nodes contribute transition matrix %*% incoming messages
      # messages[st_id, st_parent_ids[st_id], ] <- trans_mats[,, st_id] %*% apply(messages[st_child_ids[st_id,], st_id,], 2, prod)
      messages[st_id, st_parent_ids[st_id], ] <- trans_mats[,, st_id] %*% (messages[st_child_ids[st_id,1], st_id ,] * messages[st_child_ids[st_id,2], st_id, ])
      messages[st_id, st_parent_ids[st_id], ] <- messages[st_id, st_parent_ids[st_id], ] / sum(messages[st_id, st_parent_ids[st_id], ])
    }
  }

  #####################
  # Forwards sampling #
  #####################
  #When st_root is global root (a.s. the only situation where st_root may be coalescent node)
  #Need to sample deme at st_root as well as internal node!

  if (is.na(ED[st_rows[n_nodes], 2])){ #st_root == global_root
    #Deme distribution just product of incoming messages from below
    node_dist <- messages[st_child_ids[n_nodes, 1], n_nodes,] * messages[st_child_ids[n_nodes, 2], n_nodes,]
    ED[st_root_row, 5] <- sample.int(n_deme, 1, prob = node_dist)
  }

  for (st_id in (n_nodes - 1):1){ #Loop over non-subtree-leaf coalescent nodes
    node_row <- st_rows[st_id]
    if ((!is.na(st_child_rows[st_id, 1])) && (st_child_rows[st_id, 1] != node_row)){
      parent_deme <- ED[st_parent_rows[st_id], 5]

      #Deme distribution is product of transition from determined parent with product of backwards messages from node's children
      node_dist <- trans_mats[parent_deme,,st_id] * messages[st_child_ids[st_id,1], st_id,] * messages[st_child_ids[st_id,2], st_id,]
      ED[node_row, 5] <- sample.int(n_deme, 1, prob = node_dist)
    }
  }

  ##################
  # Subtree update #
  ##################
  max_label <- max(ED[,1])
  log_like <- 0 #Proposal log likelihood

  #Fill in parent edge of each node in st_labels except st_root
  for (st_id in 1 : (n_nodes - 1)){ #Skip subtree root (a.s. last in node_order)
    node_row <- st_rows[st_id]
    node_label <- st_labels[st_id]

    parent_deme <- ED[st_parent_rows[st_id], 5]
    parent_time <- ED[st_parent_rows[st_id], 6]

    # mig_path <- ECctmc::sample_path_unif(parent_deme, ED[node_row, 5], 0, st_edge_lengths[st_id], Q = fit_rates)
    mig_path <- ECctmc::sample_path_unif3(parent_deme, ED[node_row, 5], 0, st_edge_lengths[st_id],
                                          Q = fit_rates, P = trans_mats[,,st_id])


    ###### ECctmc::sample_path_unif returns t(mig_path) if 0 migrations added
    ###### BUG IN SOURCE CODE!!!
    if (nrow(mig_path) == ncol(mig_path)){
      if (mig_path[1,2] == st_edge_lengths[st_id]){
        mig_path <- t(mig_path)
      }
    }

    n_mig <- nrow(mig_path) - 2

    which_child <- which(st_children[st_parent_ids[st_id], ] == node_label)

    if (n_mig > 0){ #Add new migrations
      parent_node <- st_parent_labels[st_id]

      ED[st_parent_rows[st_id], 2 + which_child] <- max_label + 1 #Update child of node_parent
      for (k in 1 : n_mig){
        #log-probability of holding time until next migration i -> j
        log_like <- log_like +
          log(fit_rates[mig_path[k, 2], mig_path[k + 1, 2]]) + #log(f_ij)
          fit_rates[mig_path[k, 2], mig_path[k, 2]] * (mig_path[k+1, 1] - mig_path[k,1]) # - f_{i+} * (t_j - t_i)

        max_label <- max_label + 1
        ED <- rbind(ED,
                    c(max_label, #New migration ID
                      parent_node,
                      max_label + 1, NA, #Next migration ID
                      mig_path[k, 2], #New migration deme
                      parent_time + mig_path[k+1, 1], #Migration time increases leaf-wards
                      ED[node_row, 7], #Maintain same coal node parent as node_row
                      ED[st_parent_rows[st_id], 7 + which_child], NA))
        parent_node <- max_label
        ED_NI[max_label] <- nrow(ED)
      }
      ED[nrow(ED), 3] <- node_label #Update child of final migration added
      ED[node_row, 2] <- max_label #Update parent of current node
    } else{ #No migrations added
      ED[node_row, 2] <- st_parent_labels[st_id] #Update parent of current node
      ED[st_parent_rows[st_id], 2 + which_child] <- node_label #Update child of parent node
    }

    #log-probability of no further migrations between last two events on mig_path
    log_like <- log_like + fit_rates[mig_path[n_mig + 1, 2], mig_path[n_mig + 2, 2]] * (mig_path[n_mig+2, 1] - mig_path[n_mig + 1, 1])
  }

  #############################
  # Remove virtual migrations #
  ############################
  # If st_root is a migration, it is a.s. self-migration
  # Else st_root is the global root and no change needs to be made
  if (is.na(ED[st_root_row, 4])){
    parent_node <- ED[st_root_row, 2]
    parent_row <- ED_NI[parent_node]
    which_child <- which(ED[parent_row, 8:9] == ED[st_root_row, 8])
    ED[parent_row, 2 + which_child] <- ED[st_root_row, 3] #Child of parent is now child of st_root

    child_node <- ED[st_root_row, 3]
    child_row <- ED_NI[child_node]
    ED[child_row, 2] <- parent_node

    ED_rm_rows <- append(ED_rm_rows, st_root_row)
  }

  for (st_id in 1 : (n_nodes - 1)){
    current_row <- st_rows[st_id]

    if ((is.na(ED[current_row, 4])) & (!is.na(ED[current_row, 3]))){ #If current_row is migration then a.s. self-migration
      child_node <- ED[current_row, 3]
      child_row <- ED_NI[child_node]
      ED[child_row, 2] <- ED[current_row, 2] #Parent of child_node is parent(current_st_leaf)

      parent_node <- ED[current_row, 2]
      parent_row <- ED_NI[parent_node]
      which_child <- which(ED[parent_row, 3:4] == st_labels[st_id, 1])
      ED[parent_row, 2 + which_child] <- child_node

      ED_rm_rows <- append(ED_rm_rows, current_row)
    }
  }

  ED <- ED[-ED_rm_rows,]

  return(list(proposal=ED, prop_prob = log_like))
}


#' Gibbs Updates
#'
#' Performs Gibbs updates for structured coalescent parameters using conjugate gamma prior distributions
#'
#'
#' @param ED Extended data representation of a phylogeny including migration history
#' @param ED_NI Vector identifying row of each node label
#' @param n_deme Number of demes used in structured coalescent model
#' @param prior_shape Gamma shape parameter for prior
#' @param prior_rate Gamma rate parameter for prior
#'
#' @return Updated parameters
#'
#' @export

coal_rate_gibbs_update <- function(ED, ED_NI, n_deme, prior_shape = 1, prior_rate = 1){
  c <- NodeCount(ED, n_deme, ED_NI)$c
  ED_DD <- DemeDecomp(ED, n_deme, ED_NI)
  k <- ED_DD$k
  rate_consts <- colSums((k * (k-1) / 2) * ED_DD$time.increments) #t(k * (k-1) / 2) %*% DemeDecomp$time.increments

  proposal <- rgamma(n_deme, prior_shape + c, prior_rate + rate_consts)
  return(proposal)
}

#' @rdname coal_rate_gibbs_update
#'
#' @export

mig_mat_gibbs_update <- function(ED, n_deme, ED_NI, prior_shape = 1, prior_rate = 1){
  m <- NodeCount(ED, n_deme, ED_NI)$m
  ED_DD <- DemeDecomp(ED, n_deme, ED_NI)
  k <- ED_DD$k
  deme_length <- colSums(k*ED_DD$time.increments) #as.vector(t(k) %*% DemeDecomp$time.increments)

  proposal <- matrix(rgamma(n_deme^2, prior_shape + m, prior_rate + deme_length), n_deme, n_deme)
  diag(proposal) <- 0
  return(proposal)
}

#' Local DTA MCMC (Adaptive Radius)
#'
#' Runs an MCMC chain using the local DTA method with radius-based subtrees (adaptively set radius)
#'
#' @param N Total number of MCMC iterations to complete
#' @param ED Migration history to begin at in ED format (9 columns)
#' @param coal_rate Initial estimate of coalescent rates
#' @param bit_mig_mat Initial estimate of backward-in-time migration rates matrix
#' @param st_radius Initial radius for subtrees
#' @param adaptive Logical value indicating whether or not to adapt the radius within the run (default TRUE)
#' @param cr_mode Coalescent rates prior mode
#' @param cr_var Coalescent rates prior variance
#' @param mm_mode Backward-in-time migration rates prior mode
#' @param mm_var Backward-in-time migration rates prior variance
#' @param thin Thinning rate for continuous parameter posterior samples
#' @param save_trees Logical value indicating whether tree samples are saved (i.e. whether .trees file is created)
#' @param tree_thin Thinning rate for tree samples to be saved
#' @param proposal_rates Relative rates of migration history and continuous parameter updates (migration history : coalescent rates : migration rates)
#' @param adaptation_rate Rate at which adaptive MCMC varies subtree radius (unused if adaptive = FALSE)
#' @param target_accept_rate Target acceptance rate for migration history updates
#' @param output_dir Directory to output log files to
#' @param run_name Run name to save log files as
#'
#' @export

adaptive_radius_MCMC <- function(N, ED, coal_rate, bit_mig_mat,
                                 st_radius, adaptive = TRUE,
                                 cr_mode, cr_var,
                                 mm_mode, mm_var,
                                 thin = max(N/5e3, 1), save_trees = TRUE, tree_thin = thin,
                                 proposal_rates=c(1e3, 1, 1),
                                 adaptation_rate = 0.6, target_accept_rate = 0.234,
                                 output_dir, run_name = 'Local_DTA'){

  start_time <- Sys.time()

  # Convert priors to rate-shape parameterisation from mode-variance
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
  ED_NI <- NodeIndices(ED)
  ED_SC <- SC_like_C(ED, coal_rate, bit_mig_mat, ED_NI)
  mm_prior <- sum(dgamma(bit_mig_mat, mm_shape, mm_rate, log = TRUE)[-(1 + 0:(n_deme - 1) * (n_deme + 1))])
  cr_prior <- sum(dgamma(coal_rate, cr_shape, cr_rate, log = TRUE))


  # Write iteration number, likelihood, posterior and subtree radius to stdout()
  cat('Sample',
      'Likelihood',
      'Posterior',
      'Subtree radius\n',
      sep = '\t')

  cat(0, # Sample
      sprintf('%.03f', ED_SC), # Likelihood
      sprintf('%.03f', ED_SC + mm_prior + cr_prior), # posterior
      sprintf('%.03f\n', st_radius), #subtree radius (being adapted!)
      sep = '\t')

  # Set up .freq file to store move acceptance frequencies
  freq_file <- file.path(output_dir,
                         paste0(run_name, '.freq'))
  freq <- matrix(0, 2, 3,
                 dimnames = list(c("#accept", "#reject"), c("DTA", "CR", "MM")))

  # Set up .log file to store posterior continuous parameters
  log_file <- file.path(output_dir,
                        paste0(run_name, '.log'))
  cat("sample",
      "likelihood",
      "posterior",
      paste("coal_rate_", 1:n_deme, sep = "", collapse =","),
      paste("backward_migration_rate_", as.vector(outer(1:n_deme, 1:n_deme, paste, sep = "_")[-(1 + 0:(n_deme - 1) * (n_deme + 1))]), sep = ""),
      "st_radius",
      file = log_file, sep =",")

  cat(paste0("\n", 0), #sample
      ED_SC, # likelihood
      ED_SC + mm_prior + cr_prior, # posterior
      coal_rate, #coal_rate
      as.vector(bit_mig_mat)[-(1 + 0:(n_deme - 1) * (n_deme + 1))], #backward_migration_rate
      st_radius, #Subtree radius
      file = log_file, sep =",", append = TRUE)

  # Set up .trees file to store posterior sampled trees
  if (save_trees){
    tree_file <- file.path(output_dir,
                           paste0(run_name, '.trees'))
    phylo <- ed.to.phylo(ED)
    treedata <- treeio::as.treedata(phylo)
    treedata@data <- tidytree::tibble(type = paste0("\"", phylo$node.deme, "\""), node = 1:length(phylo$node.deme))

    header <- capture.output(write.beast(treedata, file=stdout(), translate=TRUE, tree.name='STATE_0')) #Generate full .trees file for initial tree - need all except final "END;"
    header[2] <- paste("[R-package scoal, ", date(), "]\n\n", sep = "") #Update package line of .trees file to scoal
    cat(header[-length(header)], file = tree_file, sep = "\n") #Save file with updated package line, omitting "END;" on final line
  }

  for (x in 1 : N){
    move_id <- sample(1:3, 1, prob = proposal_rates)

    if (move_id == 1){
      subtree <- radius_subtree(ED, ED_NI, st_radius)
      proposal <- local_DTA(subtree$EED, subtree$st_labels, fit_rates,
                            eigen_decomp = eigen_decomp, inverse_vecs = inverse_vecs)
      prop <- proposal$proposal

      if ((nrow(ED) == nrow(prop)) &&(all(na.omit(as.vector(ED == prop))))){ # Early acceptance if prop == EED
        freq[1, move_id] <- freq[1, move_id] + 1
        log_AR <- 0
      } else {
        prop_NI <- NodeIndices(prop)
        prop_SC <- SC_like_C(prop, coal_rate, bit_mig_mat, prop_NI)

        log_AR <- min(0, prop_SC - ED_SC +
                        DTA_local_likelihood(subtree$st_labels, coal_rate, bit_mig_mat)$log.likelihood -
                        proposal$prop_prob)

        if (log(runif(1)) < log_AR){
          # Accept
          freq[1, move_id] <- freq[1, move_id] + 1

          ED <- prop
          ED_NI <- prop_NI
          ED_SC <- prop_SC
        }
      }

      if (adaptive){
        # Update proposal radius
        st_radius <- exp(log(st_radius) + x^(-adaptation_rate) * (exp(log_AR) - target_accept_rate))
      }
    } else {
      if (move_id == 2){
        coal_rate <- coal_rate_gibbs_update(ED, ED_NI, n_deme, cr_shape, cr_rate)
        cr_prior <- sum(dgamma(coal_rate, cr_shape, cr_rate, log = TRUE))
      } else if (move_id == 3){
        bit_mig_mat <- mig_mat_gibbs_update(ED, ED_NI, n_deme, mm_shape, mm_rate)
        mm_prior <- sum(dgamma(bit_mig_mat, mm_shape, mm_rate, log = TRUE)[-(1 + 0:(n_deme - 1) * (n_deme + 1))])
      }

      freq[1, move_id] <- freq[1, move_id] + 1 #Gibbs move always accepted

      # Update forward-in-time rates, pre-computed likelihoods and eigendecomposition
      ED_SC <- SC_like_C(ED, coal_rate, bit_mig_mat, ED_NI)

      fit_mig_mat <- FitMigMatC(bit_mig_mat, coal_rate)
      fit_rates <- fit_mig_mat
      diag(fit_rates) <- - rowSums(fit_mig_mat)

      eigen_decomp <- eigen(fit_rates)
      eigen_vals <- eigen_decomp$values
      eigen_vecs <- eigen_decomp$vectors
      inverse_vecs <- solve(eigen_vecs)
    }

    freq[2, move_id] <- freq[2, move_id] + 1 # Increment proposal quantity for move_id


    if (x %% thin == 0){
      # Write iteration, likelihood, posterior and subtree radius to stdout()
      cat(x, # Sample
          sprintf('%.03f', ED_SC), # Likelihood
          sprintf('%.03f', ED_SC + mm_prior + cr_prior), # posterior
          sprintf('%.03f\n', st_radius), #subtree radius (being adapted!)
          sep = '\t')

      # Write continuous parameters to .log file
      cat(paste0("\n", x), #sample
          ED_SC, # likelihood
          ED_SC + mm_prior + cr_prior, # posterior
          coal_rate, #coal_rate
          as.vector(bit_mig_mat)[-(1 + 0:(n_deme - 1) * (n_deme + 1))], #Mig mat
          st_radius, #subtree radius (being adapted!)
          #ED[n_leaf + 1:(n_leaf-1), 'Deme'], #Coalescent node demes
          file = log_file, sep =",", append = TRUE)

      # Update .freq file (Overwrites existing file entirely)
      write.table(freq,
                  file = freq_file,
                  row.names = c('#ACCEPT', '#TOTAL'),
                  col.names = TRUE)
    }

    if (save_trees && (x %% tree_thin == 0)){
      # Write current tree to .trees file
      phylo <- ed.to.phylo(ED)
      treedata <- treeio::as.treedata(phylo)
      treedata@data <- tidytree::tibble(type = paste0("\"", phylo$node.deme, "\""), node = 1:length(phylo$node.deme))
      cat("\tTREE STATE_", x, " = ",
          treeio::write.beast.newick(treedata), "\n",
          file = tree_file, append = TRUE, sep = "")
    }
  }
  cat("END;", file = tree_file, append = TRUE, sep = "")

  end_time <- Sys.time()
  cat('\nElapsed time:', round(as.numeric(end_time - start_time), digits = 1), 'seconds')
}
