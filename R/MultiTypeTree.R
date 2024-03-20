#' MultiTypeTree Node Retype
#'
#' Generates a new migration history on a subtree drawn using MTT_st_coal_node()
#'
#' @param ED Extended data representation of phylogenetic tree and initial migration history
#' @param st_labels Reduced extended data structure consisting of the rows of ED corresponding to coalescent nodes in the selected subtree only
#' @param bit_rates Transition matrix of a Markov process corresponding to a process with rates given by the backwards-in-time migration matrix
#' @param ED_NI Vector of row indices corresponding to which row of ED corresponds to each node label
#'
#' @return  List consisting of proposal, a structured phylogeny with updates made to the selected subtree, and prop_prob, the probability of the given update
#'
#' @export

MTT_node_retype <- function(ED, st_labels, bit_rates, ED_NI = NodeIndices(ED), eigen_vals = eigen(bit_rates)$values, eigen_vecs = eigen(bit_rates)$vectors, inverse_vecs = solve(eigen_vecs)){
  st_root <- unname(st_labels[1,1]) #st_labels[1,] always corresponds to st_root by construction
  st_leaves <- st_labels[!(st_labels[,8] %in% st_labels[,1]),1]

  #Identify migrations with child and parent inside subtree
  rm_mig_pos <- (is.na(ED[,9])) & (ED[,8] %in% st_labels[,1]) & (ED[,7] %in% st_labels[,1])
  ED_rm_rows <- ED_NI[ED[rm_mig_pos,1]]

  internal_coal_nodes <- st_labels[!(st_labels[,1] %in% c(st_root, st_leaves)),1]

  if (is.na(st_labels[1,7])){ #If st_root is global root and not fully determined
    if (all(st_labels[1, 8:9] %in% st_labels[,1])){
      internal_coal_nodes <- append(internal_coal_nodes, st_root)
    }
  }

  int_coal_node_row <- ED_NI[internal_coal_nodes]

  #### Sample new deme at internal coalescent nodes
  n_deme <- nrow(bit_rates)
  ED[int_coal_node_row, 5] <- sample(1:n_deme, length(internal_coal_nodes), TRUE)

  ### Complete migration history
  max_label <- ED[nrow(ED), 1] #max(ED[,1]) is label of final row
  log_like <- 0 #Proposal log likelihood

  for (row_id in 2 : nrow(st_labels)){ #Row 1 = st_root
    node_row <- ED_NI[st_labels[row_id, 1]]
    node_deme <- ED[node_row, 5]

    parent_row <- ED_NI[st_labels[row_id, 7]]
    parent_deme <- ED[parent_row, 5]

    edge_length <- ED[node_row, 6] - ED[parent_row, 6]

    trans_mat <- Re(eigen_vecs %*% diag(exp(eigen_vals * edge_length))  %*% inverse_vecs)

    mig_path <- ECctmc::sample_path_unif3(node_deme, parent_deme,
                                          0, edge_length,
                                          Q = bit_rates, P = trans_mat)



    ###### ECctmc::sample_path_unif returns t(mig_path) if 0 migrations added
    ###### BUG IN SOURCE CODE!!!
    if (nrow(mig_path) == ncol(mig_path)){
      if (mig_path[1,2] == edge_length){
        mig_path <- t(mig_path)
      }
    }

    n_mig <- nrow(mig_path) - 2
    which_child <- which(ED[parent_row, 8:9] == ED[node_row, 1])

    #log(bit_rates[i,i]) = NaN; removed by sum(..., na.rm=TRUE)
    log_like <- log_like +
      suppressWarnings(sum(log(bit_rates[mig_path[-(1 + 0:n_mig), 2] + n_deme * (mig_path[-1, 2] - 1)]), na.rm=TRUE)) + #sum(log(\lambda_{ij}))
      sum(bit_rates[mig_path[-nrow(mig_path), 2] + n_deme * (mig_path[-nrow(mig_path), 2] - 1)] * #lambda_{i+}
            (mig_path[-nrow(mig_path),1] - mig_path[-1, 1])) - # * (t_j - t_i)
      # sum(bit_rates[mig_path[-(1 + 0:n_mig), 2] + n_deme * (mig_path[-(1 + 0:n_mig), 2] - 1)] * (mig_path[-1,1] - mig_path[-(1 + 0:n_mig), 1])) - # - sum(\lambda_{i+} * (t_j - t_i))
      log(trans_mat[node_deme, parent_deme]) #Dividing by conditional probability

    if (n_mig > 0){ #Add new migrations
      new_migrations <- cbind(max_label + 1:n_mig, #New migration IDs
                              max_label + 1 + 1:n_mig, #Parent ID = next new migration ID
                              max_label - 1 + 1:n_mig, NA, #Child ID = previous new migration ID
                              mig_path[1 + 1:n_mig, 2], #New demes come from mig_path col 2 (rows 2 : nrow() - 1)
                              ED[node_row, 6] - mig_path[1 + 1:n_mig, 1],  #New times come from mig_path col 1 (rows 2 : nrow() - 1)
                              ED[parent_row, 1], #Parent coal node of current node
                              ED[node_row, 1], NA) #Child coal node of all new migrations is current node

      new_migrations[1, 3] <- ED[node_row, 1] #Correct first migration child
      new_migrations[n_mig, 2] <- ED[parent_row, 1] #Correct last migration parent

      ED[node_row, 2] <- max_label + 1 #Correct parent of current node
      ED[parent_row, 2 + which_child] <- max_label + n_mig #Correct child of parent_node

      ED <- rbind(ED, new_migrations)
      max_label <- max_label + n_mig
    } else { #Correct current and parent nodes to remove migrations
      ED[node_row, 2] <- ED[parent_row, 1] #Parent of current row is parent_row
      ED[parent_row, 2 + which_child] <- ED[node_row, 1] #Child of parent_row is current_row
    }

  }
  ### Remove old migration events
  if (length(ED_rm_rows) > 0){
    ED <- ED[-ED_rm_rows, ]
  }

  return(list(proposal=ED, prop_prob = log_like))
}

#' MultiTypeTree MCMC
#'
#' Runs a MCMC chain drawing joint samples of migration histories, backwards-in-time migration rates and coalescent rates.
#' MCMC chain runs with the node retype operator from MultiTypeTree for migration history updates and Gamma Gibbs moves for migration and coalescent rates
#'
#' @param N Number of MCMC iterations
#' @param ED Extended data representation of phylogenetic tree and initial migration history
#' @param coal_rate Vector of initial coalescent rates for the MCMC chain
#' @param bit_mig_mat Matrix of initial backwards-in-time migration rates for the MCMC chain
#' @param st_depth Number of generations below a coalescent node to include in the subtree for the node retype move
#' @param cr_mode Mode for prior on coalescent rates
#' @param cr_var Variance for prior on coalescent rates
#' @param mm_mode Mode for prior on backwards-in-time migration rates
#' @param mm_var Variance for prior on backwards-in-time migration rates
#' @param thin Thinning increment for the MCMC chain. A sample is saved every \texttt{thin} iterations
#' @param prop_rates Relative rates for proposals of migration history updates to coalescent rate updates to backwards-in-time migration rates updates
#' @param output_dir Directory to output files to for (thinned) results from the MCMC chain
#' @param run_name Name for run in output_dir
#'
#' @return Output files output_dir/run_name.freq, output_dir/run_name.log and output_dir/run_name.trees giving the acceptance rates of each MCMC move, thinned samples of coalescent and migration rates, and thinned samples of trees respectively
#'
#' @export

MTT_node_retype_MCMC <- function(N = 1e6,
                                 ED, coal_rate, bit_mig_mat,
                                 st_depth = 1, cr_mode = 1, cr_var = 1, mm_mode = 0.05, mm_var = 0.5,
                                 thin = 1e3, prop_rates = c(100, 1, 1),
                                 output_dir = '~', run_name = 'MTT_R'){
  pb <- txtProgressBar(0, N, style = 3)

  # Prior parameters
  cr_rate <- (cr_mode + sqrt(cr_mode^2 + 4 * cr_var))/(2 * cr_var)
  cr_shape <- 1 + cr_mode * cr_rate
  mm_rate <- (mm_mode + sqrt(mm_mode^2 + 4 * mm_var))/(2 * mm_var)
  mm_shape <- 1 + mm_mode * mm_rate

  bit_rates <- bit_mig_mat
  diag(bit_rates) <- - rowSums(bit_mig_mat)
  n_deme <- nrow(bit_rates)

  eigen_decomp <- eigen(bit_rates)
  eigen_vecs <- eigen_decomp$vectors
  eigen_vals <- eigen_decomp$values
  inverse_vecs <- solve(eigen_vecs)

  # Initial likelihoods and priors
  ED_NI <- NodeIndices(ED)
  ED_SC <- SC_like_C(ED, coal_rate, bit_mig_mat, ED_NI)

  mm_prior <- sum(dgamma(bit_mig_mat, mm_shape, mm_rate, log = TRUE)[-(1 + 0:(n_deme - 1) * (n_deme + 1))])
  cr_prior <- sum(dgamma(coal_rate, cr_shape, cr_rate, log = TRUE))

  # Initialise output files
  freq_file <- file.path(output_dir,
                         paste0(run_name, ".freq"))
  freq <- matrix(0, 2, 3,
                 dimnames = list(c("#accept", "#reject"), c("Node_Retype", "CR", "MM")))

  # log_file stores mig rates and coal rates etc
  log_file <- file.path(output_dir,
                        paste0(run_name, ".log"))
  # Column names
  cat("sample",
      "likelihood",
      "posterior",
      paste("coal_rate_", 1:n_deme, sep = "", collapse =","),
      paste("backward_migration_rate_", as.vector(outer(1:n_deme, 1:n_deme, paste, sep = "_")[-(1 + 0:(n_deme - 1) * (n_deme + 1))]), sep = ""),
      file = log_file, sep =",")

  # First row
  cat(paste0("\n", 0), #sample
      ED_SC, # likelihood
      ED_SC + mm_prior + cr_prior, # posterior
      coal_rate, #coal_rate
      as.vector(bit_mig_mat)[-(1 + 0:(n_deme - 1) * (n_deme + 1))], #Mig mat
      file = log_file, sep =",", append = TRUE)

  # tree_file stores posterior sampled trees
  tree_file <- file.path(output_dir,
                         paste0(run_name, ".trees"))
  phylo <- ed.to.phylo(ED)
  treedata <- treeio::as.treedata(phylo)
  treedata@data <- tidytree::tibble(type = paste0("\"", phylo$node.deme, "\""), node = 1:length(phylo$node.deme))

  # .nexus file header
  header <- capture.output(write.nexus(phylo, file = stdout(), translate = TRUE))
  header[2] <- paste("[R-package scoal, ", date(), "]\n\n", sep = "")
  cat(header[-(length(header)-0:1)], file = tree_file, sep = "\n")

  #Initial tree
  cat("\tTREE STATE_0 = ",
      treeio::write.beast.newick(treedata), "\n",
      file = tree_file, append = TRUE, sep = "")

  for (x in 1 : N){
    move_id <- sample(1:3, 1, prob = prop_rates)

    if (move_id == 1){
      subtree <- coal_node_subtree(ED, st_depth)
      proposal <- MTT_node_retype(ED, subtree$st_labels, bit_rates, ED_NI, eigen_vals, eigen_vecs, inverse_vecs)

      # Early acceptance if prop == ED
      if ((nrow(ED) == nrow(proposal$proposal)) && (all(na.omit(as.vector(ED == proposal$proposal))))){
        freq[1, move_id] <- freq[1, move_id] + 1
      } else {
        prop_NI <- NodeIndices(proposal$proposal)
        prop_SC <- SC_like_C(proposal$proposal, coal_rate, bit_mig_mat, prop_NI)
        log_AR <- min(0, Re(prop_SC - ED_SC +
                              MTT_local_like(ED, subtree$st_labels, bit_rates, ED_NI, eigen_vals, eigen_vecs, inverse_vecs)$log.likelihood -
                              MTT_local_like(proposal$proposal, subtree$st_labels, bit_rates, prop_NI, eigen_vals, eigen_vecs, inverse_vecs)$log.likelihood
        ))

        if (log(runif(1)) < log_AR){
          # Accept
          freq[1, move_id] <- freq[1, move_id] + 1

          ED <- proposal$proposal
          ED_NI <- prop_NI
          ED_SC <- prop_SC
        }
      }
    } else {
      if (move_id == 2){
        coal_rate <- scaled_coal_rate_gibbs_update(ED, 1, n_deme, ED_NI, cr_shape, cr_rate)
        cr_prior <- sum(dgamma(coal_rate, cr_shape, cr_rate, log = TRUE))
      } else if (move_id == 3){
        bit_mig_mat <- scaled_mig_mat_gibbs_update(ED, 1, n_deme, ED_NI, mm_shape, mm_rate)

        bit_rates <- bit_mig_mat
        diag(bit_rates) <- - rowSums(bit_mig_mat)

        eigen_decomp <- eigen(bit_rates)
        eigen_vecs <- eigen_decomp$vectors
        eigen_vals <- eigen_decomp$values
        inverse_vecs <- solve(eigen_vecs)

        mm_prior <- sum(dgamma(bit_mig_mat, mm_shape, mm_rate, log = TRUE)[-(1 + 0:(n_deme - 1) * (n_deme + 1))])
      }

      ED_SC <- SC_like_C(ED, coal_rate, bit_mig_mat, ED_NI)
      freq[1, move_id] <- freq[1, move_id] + 1
    }

    freq[2, move_id] <- freq[2, move_id] + 1

    setTxtProgressBar(pb, x)
    if (x %% thin == 0){ #x (mod thin) = 0, i.e. thin-many iterations have passed
      cat(paste0("\n", x), #sample
          ED_SC, # likelihood
          ED_SC + mm_prior + cr_prior, # posterior
          coal_rate, #coal_rate
          as.vector(bit_mig_mat)[-(1 + 0:(n_deme - 1) * (n_deme + 1))], #Mig mat
          file = log_file, sep =",", append = TRUE)

      phylo <- ed.to.phylo(ED)
      treedata <- treeio::as.treedata(phylo)
      treedata@data <- tidytree::tibble(type = paste0("\"", phylo$node.deme, "\""), node = 1:length(phylo$node.deme))
      cat("\tTREE STATE_", x, " = ",
          treeio::write.beast.newick(treedata), "\n",
          file = tree_file, append = TRUE, sep = "")

      write.table(freq,
                  file = freq_file,
                  row.names = c('#ACCEPT', '#TOTAL'),
                  col.names = TRUE)
    }

  }
  close(pb)
  cat("END;",
      file = tree_file, append = TRUE, sep = "")
}
