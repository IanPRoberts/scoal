#' DTA Sampling
#'
#' Samples N migration histories under the DTA model deme fixed at all leaves using belief propagation
#'
#'
#' @param ED Extended data representation of a phylogeny including migration history
#' @param mig_mat Forward-in-time migration matrix for the phylogeny
#' @param time_scale Time scale for the migration matrix
#' @param N Number of migration histories to sample
#' @param parallel Logical, whether or not to generate proposals in parallel
#' @param mc.cores (optional) If parallel is TRUE, number of cores to parallelise over
#'
#' @return List of proposed DTA migration histories with the same leaves as input ED
#'
#' @export

DTA_sampling <- function(ED, mig_mat, time_scale = 1, N = 1, parallel = FALSE, mc.cores = NA){
  options(mc.cores = 1)
  if (N > 1){
    if (parallel){
      if (is.na(mc.cores)){
        options(mc.cores = parallel::detectCores() /2)
      }
    }
  }

  fit_rates <- mig_mat * time_scale
  diag(fit_rates) <- - rowSums(mig_mat)

  top_ED <- strip.history(ED)
  node_indices <- NodeIndicesC(top_ED)
  n <- (dim(top_ED)[1] + 1)/2
  n_deme <- dim(mig_mat)[1]

  ### Transition matrices
  eigen_decomp <- eigen(fit_rates)
  V <- eigen_decomp$vectors
  V_inv <- solve(V)
  edge_lengths <- top_ED[,6] - top_ED[node_indices[top_ED[,2]], 6] #edge_lengths[i] = time between node at row i and its parent; edge_lengths[root] = NA

  trans_mats <- array(sapply(edge_lengths, function(x) V %*% diag(exp(eigen_decomp$values * x)) %*% V_inv),
                      dim = c(n_deme, n_deme, 2*n-1))

  ### Backwards messages
  messages <- array(0, dim = c(2*n-1, 2*n-1, n_deme))
  node_order <- order(top_ED[,6], decreasing = TRUE) #Ordering of coalescent nodes from lowest to root. Only need do once for entire MCMC!

  for (i in 1 : (2*n-2)){ #Final entry is root node - no messages here
    node_row <- node_order[i]
    node_parent <- top_ED[node_row, 2]
    node_parent_row <- node_indices[node_parent]

    if (is.na(top_ED[node_row, 3])){ #Leaf node
      leaf_deme <- top_ED[node_row, 5]
      messages[node_row, node_parent_row, ] <- trans_mats[ ,leaf_deme, node_row] #Fully determined leaf nodes
      #messages[node_row, node_parent_row, ] <- trans_mats[,, node_row] %*% leaf_dists[,node_row] #Partially determined leaf nodes - need distributions in matrix leaf_dists
    } else{ #Coalescent node
      node_children <- top_ED[node_row, 3:4]
      node_children_rows <- node_indices[node_children]
      messages[node_row, node_parent_row, ] <- trans_mats[,, node_row] %*% apply(messages[node_children_rows, node_row,], 2, prod)
    }
  }

  proposal_func <- function(x){
    prop_ED <- top_ED #New proposal
    max_label <- max(prop_ED[,1])

    for (i in (2*n-1):1){
      node_row <- node_order[i]

      if (!is.na(prop_ED[node_row, 3])){ #Coalescent node
        node_children <- prop_ED[node_row, 3:4]
        node_children_rows <- node_indices[node_children]

        node_dist <- apply(messages[node_children_rows, node_row,], 2, prod) #Product of root-wards child messages

        if (!is.na(prop_ED[node_row, 2])){ #Non-root coalescent node
          parent_row <- node_indices[prop_ED[node_row,2]]
          parent_deme <- prop_ED[parent_row, 5]
          node_dist <- node_dist * trans_mats[parent_deme,,node_row] #Multiply by leaf-wards parent message
        }

        prop_ED[node_row, 5] <- sample(1:n_deme, 1, prob = node_dist) #Sample node distribution with probability proportional to node_dist
      }

      if (!is.na(prop_ED[node_row, 2])){ #Non-root nodes
        parent_row <- node_indices[prop_ED[node_row,2]] #Recalculate parent_row in case of leaf node
        mig_path <- ECctmc::sample_path(prop_ED[parent_row, 5], prop_ED[node_row, 5], 0, edge_lengths[node_row], Q = fit_rates)
        n_mig <- dim(mig_path)[1] - 2

        if (n_mig > 0){
          parent_node <- prop_ED[parent_row, 1]
          which_child <- which(prop_ED[parent_row, 3:4] == prop_ED[node_row])
          prop_ED[parent_row, 2 + which_child] <- max_label + 1

          for (k in 1 : n_mig){
            prop_ED <- rbind(prop_ED,
                             c(max_label + 1, parent_node, max_label + 2, NA, mig_path[k, 2], prop_ED[parent_row, 6] + mig_path[k+1, 1]))
            max_label <- max_label + 1
            parent_node <- max_label
          }
          prop_ED[dim(prop_ED)[1], 3] <- prop_ED[node_row, 1]
          prop_ED[node_row, 2] <- max_label
        }
      }
    }
    return(prop_ED)
  }

  if (parallel){
    proposals <- parallel::mclapply(1:N, proposal_func)
  } else {
    proposals <- lapply(1:N, proposal_func)
  }
  return(proposals)
}

#' Local DTA Update
#'
#' Updates a migration history under the DTA model locally around a single coalescent node
#'
#'
#' @param ED Extended data representation of a phylogeny including migration history
#' @param mig_mat Forward-in-time migration matrix for the phylogeny
#' @param time_scale Time scale for the migration matrix
#' @param selected_node Label for node proposal will be centered on (optional)
#'
#' @return Proposal ED
#'
#' @export

local_DTA <- function(ED, mig_mat, time_scale, selected_node = NA){
  #############################################
  #   Select subtree & sample internal deme   #
  #############################################

  n_deme <- dim(mig_mat)[1]
  fit_rates <- mig_mat * time_scale
  diag(fit_rates) <- - rowSums(fit_rates)
  top_ED <- strip.history(ED)
  node_indices <- NodeIndicesC(top_ED)
  ED_node_indices <- NodeIndicesC(ED)

  if (is.na(selected_node)){
    leaf_nodes <- top_ED[is.na(top_ED[,3]),1]
    coal_nodes <- top_ED[-node_indices[leaf_nodes],1]
    selected_node <- sample(coal_nodes, 1)
  }

  selected_row <- node_indices[selected_node]
  child_rows <- ED_node_indices[top_ED[selected_row, 3:4]]
  child_edge_lengths <- ED[child_rows, 6] - ED[selected_row, 6]

  trans_mats <- lapply(1:2, function(x) expm::expm(child_edge_lengths[x] * fit_rates))
  node_dist <- trans_mats[[1]][, ED[child_rows[[1]], 5]] * trans_mats[[2]][, ED[child_rows[[2]], 5]]
  node_dist <- node_dist / sum(node_dist)


  if (!is.na(top_ED[selected_row, 2])){ #Non-root node
    parent_row <- ED_node_indices[top_ED[selected_row, 2]]
    par_trans_mat <- expm::expm((ED[selected_row, 6] - ED[parent_row, 6]) * fit_rates)
    node_dist <- node_dist * par_trans_mat[ED[parent_row, 5],]
  }

  node_deme <- sample(1:n_deme, 1, prob = node_dist) #New deme at selected_row

  #######################################
  #   Construct new migration history   #
  #######################################

  prop_ED <- ED
  prop_node_indices <- NodeIndicesC(prop_ED)
  prop_selected_row <- prop_node_indices[selected_node]
  prop_child_rows <- prop_node_indices[top_ED[selected_row, 3:4]]
  max_label <- max(prop_ED[,1])

  prop_ED[prop_selected_row, 5] <- node_deme
  selected_time <- prop_ED[prop_selected_row, 6]

  rm_rows <- numeric(0)
  for (i in 1 : 2){
    active_row <- prop_child_rows[i]

    edge_length <- prop_ED[active_row, 6] - selected_time
    mig_path <- ECctmc::sample_path(node_deme, prop_ED[active_row, 5], 0, edge_length, Q = fit_rates)
    n_mig <- dim(mig_path)[1] - 2

    while (prop_ED[active_row, 2] != prop_ED[prop_selected_row, 1]){ # Identify migrations between selected and child[i]
      rm_rows <- c(rm_rows, prop_node_indices[prop_ED[active_row, 2]])
      active_row <- prop_node_indices[prop_ED[active_row, 2]]
    }

    if (n_mig > 0){ #Add new migrations
      parent_node <- selected_node
      prop_ED[prop_selected_row, 2 + i] <- max_label + 1

      for (k in 1 : n_mig){
        prop_ED <- rbind(prop_ED,
                         c(max_label + 1, parent_node, max_label + 2, NA, mig_path[k, 2], prop_ED[prop_selected_row, 6] + mig_path[k+1, 1]))
        max_label <- max_label + 1
        parent_node <- max_label
      }
      prop_ED[dim(prop_ED)[1], 3] <- top_ED[child_rows[i], 1]
      prop_ED[prop_child_rows[i], 2] <- max_label
    } else {
      prop_ED[prop_child_rows[i], 2] <- selected_node
      prop_ED[prop_selected_row, 2 + i] <- prop_child_rows[i]
    }
  }

  if (!is.na(top_ED[selected_row, 2])){ #Repeat for parent edge if not root node
    active_row <- prop_selected_row
    prop_parent_row <- prop_node_indices[top_ED[parent_row, 1]]

    edge_length <- selected_time - prop_ED[prop_parent_row, 6]
    mig_path <- ECctmc::sample_path(prop_ED[prop_parent_row, 5], node_deme, 0, edge_length, Q = fit_rates)
    n_mig <- dim(mig_path)[1] - 2

    while (prop_ED[active_row, 2] != prop_ED[prop_parent_row, 1]){ #Identify migrations between parent and node
      rm_rows <- c(rm_rows, prop_node_indices[prop_ED[active_row, 2]])
      active_row <- prop_node_indices[prop_ED[active_row, 2]]
    }

    if (n_mig > 0){ #Add new migrations
      parent_node <- prop_ED[prop_parent_row, 1]
      which.child <- which(prop_ED[prop_parent_row, 3:4] %in% c(selected_node, prop_ED[rm_rows, 1]))
      prop_ED[prop_parent_row, 2 + which.child] <- max_label + 1

      for (k in 1 : n_mig){
        prop_ED <- rbind(prop_ED,
                         c(max_label + 1, parent_node, max_label + 2, NA, mig_path[k, 2], prop_ED[prop_parent_row, 6] + mig_path[k+1, 1]))
        max_label <- max_label + 1
        parent_node <- max_label
      }
      prop_ED[dim(prop_ED)[1], 3] <- selected_node
      prop_ED[prop_selected_row, 2] <- max_label
    } else {
      prop_ED[prop_selected_row, 2] <- prop_ED[prop_parent_row, 1]
      which.child <- which(prop_ED[prop_parent_row, 3:4] %in% prop_ED[rm_rows, 1])
      prop_ED[prop_parent_row, 2 + which.child] <- selected_node
    }
  }
  if (length(rm_rows) > 0) prop_ED <- prop_ED[-rm_rows,]
  return(prop_ED)
}

#' Local DTA Update for EED data structure
#'
#' Updates a migration history using the DTA model for adjacent branches to a coalescent node via belief propagation
#'
#'
#' @param EED EED representation of a phylogeny including migration history
#' @param fit_mig_mat Forward-in-time (relative) migration matrix for the phylogeny
#' @param time_scale Time scale for the migration matrix
#' @param selected_node Label for node proposal will be centered on (optional)
#'
#' @return Proposal ED
#'
#' @export

EED_local_DTA <- function(EED, fit_mig_mat, time_scale, selected_node = NA){
  #############################################
  #   Select subtree & sample internal deme   #
  #############################################

  n_deme <- dim(fit_mig_mat)[1]
  fit_rates <- fit_mig_mat * time_scale
  diag(fit_rates) <- - rowSums(fit_rates)
  node_indices <- NodeIndicesC(EED)

  if (is.na(selected_node)){
    selected_node <- sample(EED[!is.na(EED[,4]), 1], 1)
  }

  selected_row <- node_indices[selected_node]
  child_rows <- node_indices[EED[selected_row, 8:9]]
  child_edge_lengths <- EED[child_rows, 6] - EED[selected_row, 6]

  trans_mats <- lapply(1:2, function(x) expm::expm(child_edge_lengths[x] * fit_rates))
  node_dist <- trans_mats[[1]][, EED[child_rows[[1]], 5]] * trans_mats[[2]][, EED[child_rows[[2]], 5]]
  node_dist <- node_dist / sum(node_dist)


  if (!is.na(EED[selected_row, 7])){ #Non-root node
    parent_row <- node_indices[EED[selected_row, 7]]
    par_trans_mat <- expm::expm((EED[selected_row, 6] - EED[parent_row, 6]) * fit_rates)
    node_dist <- node_dist * par_trans_mat[EED[parent_row, 5],]
  }

  node_deme <- sample(1:n_deme, 1, prob = node_dist) #New deme at selected_row

  #######################################
  #   Construct new migration history   #
  #######################################

  prop <- EED
  max_label <- max(prop[,1])

  prop[selected_row, 5] <- node_deme
  selected_time <- prop[selected_row, 6]

  rm_rows <- which(xor(prop[,8] %in% prop[selected_row, c(1,8:9)],
                       prop[,9] %in% prop[selected_row, c(1,8:9)]))
  if (!is.na(prop[selected_row, 7])){
    rm_rows <- rm_rows[prop[rm_rows, 1] != prop[selected_row, 7]]
  }

  for (i in 1 : 2){
    active_row <- child_rows[i]
    edge_length <- prop[active_row, 6] - selected_time
    mig_path <- ECctmc::sample_path(node_deme, prop[active_row, 5], 0, edge_length, Q = fit_rates)
    n_mig <- dim(mig_path)[1] - 2

    if (n_mig > 0){ #Add new migrations
      parent_node <- selected_node
      prop[selected_row, 2 + i] <- max_label + 1

      for (k in 1 : n_mig){
        prop <- rbind(prop,
                      c(max_label + 1, parent_node, max_label + 2, NA, mig_path[k, 2], prop[selected_row, 6] + mig_path[k+1, 1], selected_node, prop[selected_row, 7 + i], NA))
        max_label <- max_label + 1
        parent_node <- max_label
      }
      prop[dim(prop)[1], 3] <- EED[child_rows[i], 1]
      prop[child_rows[i], 2] <- max_label
    } else {
      prop[child_rows[i], 2] <- selected_node
      prop[selected_row, 2 + i] <- child_rows[i]
    }
  }

  if (!is.na(EED[selected_row, 2])){ #Repeat for parent edge if not root node
    edge_length <- selected_time - prop[parent_row, 6]
    mig_path <- ECctmc::sample_path(prop[parent_row, 5], node_deme, 0, edge_length, Q = fit_rates)
    n_mig <- dim(mig_path)[1] - 2
    if (prop[parent_row, 8] == selected_node) which.child <- 1 else which.child <- 2 #which.child <- which(prop[parent_row, 8:9] == selected_node)

    if (n_mig > 0){ #Add new migrations
      parent_node <- prop[parent_row, 1]
      prop[parent_row, 2 + which.child] <- max_label + 1

      for (k in 1 : n_mig){
        prop <- rbind(prop,
                      c(max_label + 1, parent_node, max_label + 2, NA, mig_path[k, 2], prop[parent_row, 6] + mig_path[k+1, 1], prop[selected_row, 7], selected_node, NA))
        max_label <- max_label + 1
        parent_node <- max_label
      }
      prop[dim(prop)[1], 3] <- selected_node
      prop[selected_row, 2] <- max_label
    } else {
      prop[selected_row, 2] <- prop[parent_row, 1]
      prop[parent_row, 2 + which.child] <- selected_node
    }
  }
  if (length(rm_rows) > 0) prop <- prop[-rm_rows,]
  return(list(proposal = prop, node_dist = node_dist))
}

#' Local DTA Update for EED data structure
#'
#' Updates a migration history using the DTA model for adjacent branches to a coalescent node via rejection sampling
#'
#'
#' @param EED EED representation of a phylogeny including migration history
#' @param fit_mig_mat Forward-in-time (relative) migration matrix for the phylogeny
#' @param time_scale Time scale for the migration matrix
#' @param selected_node Label for node proposal will be centered on (optional)
#'
#' @return Proposal ED
#'
#' @export

local_DTA_rejection <- function(EED, coal_node, fit_mig_mat, node_indices, max_attempts = 1e5){
  event_rates <- rowSums(fit_mig_mat)
  coal_row <- node_indices[coal_node]

  subtree_parent <- EED[coal_row, 7]
  sp_row <- node_indices[subtree_parent]

  subtree_children <- EED[coal_row, 8:9]
  sc_rows <- node_indices[subtree_children]

  coal_deme <- EED[coal_row, 5]
  parent_deme <- EED[sp_row, 5]
  child_demes <- EED[sc_rows, 5]

  coal_time <- EED[coal_row, 6]
  parent_time <- EED[sp_row, 6]
  child_times <- EED[sc_rows, 6]

  rm_rows <- which(xor(EED[,8] %in% EED[coal_row, c(1,8:9)],
                       EED[,9] %in% EED[coal_row, c(1,8:9)]))
  if (!is.na(EED[coal_row, 7])){
    rm_rows <- rm_rows[EED[rm_rows, 1] != EED[coal_row, 7]]
  }

  accept <- FALSE
  count <- 0

  while ((!accept) && (count < max_attempts)){
    prop <- EED
    old_nrow <- dim(prop)[1]
    max_label <- max(EED[,1])

    if (!is.na(subtree_parent)){
      current_node <- subtree_parent
      current_deme <- parent_deme
      current_time <- parent_time + rexp(1, event_rates[current_deme])

      if (current_time < coal_time){
        # Update sp_row child if migration is added
        which.child <- which(EED[sp_row, 8:9] == coal_node)
        prop[sp_row, 2 + which.child] <- max_label + 1
      }

      #Edge above coal_node
      while (current_time < coal_time){
        # Add migration node if not passed
        prop <- rbind(prop, c(max_label + 1, current_node, max_label + 2, NA, current_deme, current_time, subtree_parent, coal_node, NA))

        max_label <- max_label + 1
        current_node <- max_label
        current_deme <- sample.int(n_deme, 1, prob = fit_mig_mat[current_deme, ])
        current_time <- current_time + rexp(1, event_rates[current_deme])
      }
      nrow <- dim(prop)[1]
      if (nrow > old_nrow){
        #Update coal_row if migration(s) added
        prop[coal_row, 2] <- max_label
        prop[coal_row, 5] <- current_deme

        #Update final migration child
        prop[dim(prop)[1], 3] <- coal_node
        old_nrow <- nrow
      }else {
        prop[coal_row, 2] <- subtree_parent

        which.child <- which(prop[sp_row, 8:9] == coal_node)
        prop[sp_row, 2 + which.child] <- coal_node
      }
    }

    for (child_id in 1 : 2){
      current_deme <- prop[coal_row, 5]
      current_node <- coal_node
      current_time <- coal_time + rexp(1, event_rates[current_deme])

      if (current_time < child_times[child_id]){
        prop[coal_row, 2 + child_id] <- max_label + 1
      }

      #Edge below coal_node
      while (current_time < child_times[child_id]){
        # Add migration node if not past
        prop <- rbind(prop, c(max_label + 1, current_node, max_label + 2, NA, current_deme, current_time, coal_node, subtree_children[child_id], NA))

        max_label <- max_label + 1
        current_node <- max_label
        current_deme <- sample.int(n_deme, 1, prob = fit_mig_mat[current_deme, ])
        current_time <- current_time + rexp(1, event_rates[current_deme])
      }

      if (current_deme != child_demes[child_id]){
        accept <- FALSE
        break
      }
      accept <- TRUE

      nrow <- dim(prop)[1]
      if (nrow > old_nrow){
        # Update subtree_children[child_id]
        prop[sc_rows[child_id], 2] <- max_label
        prop[sc_rows[child_id], 5] <- current_deme

        # Update child of final migration added
        prop[nrow, 3] <- subtree_children[child_id]

        old_nrow <- nrow
      } else {
        prop[coal_row, 2 + child_id] <- subtree_children[child_id]
        prop[sc_rows[child_id], 2] <- coal_node
      }
    }
    count <- count + 1
  }
  if (length(rm_rows) > 0){
    prop <- prop[-rm_rows,]
  }

  if (count == max_attempts) stop("Attempts exceeded max_attempts")
  return(prop)
}

#' DTA Sampling
#'
#' Samples a migration history under the DTA model with deme fixed at all leaves using rejection sampling
#'
#'
#' @param ED Extended ED representation of a phylogeny including migration history
#' @param fit_mig_mat Forward-in-time migration matrix for the phylogeny
#' @param time_scale Time scale for the migration matrix
#'
#' @return List of proposed DTA migration histories with the same leaves as input ED
#'
#' @export

DTA_rejection <- function(ED, fit_mig_mat, time_scale, node_indices, max_attempts = 1e4){
  n_deme <- dim(fit_mig_mat)[1]
  event_rates <- rowSums(fit_mig_mat)
  accept <- FALSE

  top_EED <- ed.to.eed(strip.history(ED))
  root_row <- which(is.na(top_EED[,2]))
  leaf_rows <- which(is.na(top_EED[,3]))
  top_max_label <- max(top_EED[,1])

  for (x in 1 : max_attempts){
    prop <- top_EED
    active_rows <- root_row
    root_deme <- sample.int(n_deme, 1) #Uniform root deme
    prop[root_row, 5] <- root_deme

    max_label <- top_max_label

    while (length(active_rows) > 0){
      for (row in active_rows){
        child_rows <- node_indices[na.omit(prop[row, 8:9])]
        row_time <- prop[row, 6]
        which.child <- 1

        for (child_row in child_rows){
          current_node <- prop[row, 1]
          current_deme <- prop[row, 5]
          current_time <- row_time + rexp(1, event_rates[current_deme])
          child_time <- prop[child_row, 6]

          if (current_time < child_time){
            # Update row_children if migrations being added
            prop[row, 2 + which.child] <- max_label + 1

            while (current_time < child_time){
              # Add migration node if not past
              prop <- rbind(prop,
                            c(max_label + 1, current_node, max_label + 2, NA, current_deme, current_time, prop[row, 1], prop[child_row, 1], NA))

              max_label <- max_label + 1
              current_node <- max_label
              current_deme <- sample.int(n_deme, 1, prob = fit_mig_mat[current_deme, ])
              current_time <- current_time + rexp(1, event_rates[current_deme])
            }

            #Update child_row parent and deme if migration(s) added
            prop[child_row, 2] <- max_label
            prop[child_row, 5] <- current_deme

            #Update final migration child
            prop[nrow(prop), 3] <- prop[child_row, 1]
          } else {
            #Child row has current node as parent, current deme as deme
            prop[child_row, 2] <- current_node
            prop[child_row, 5] <- current_deme
            prop[row, 2 + which.child] <- prop[child_row, 1]
          }

          which.child <- which.child + 1
        }

      }
      active_rows <- node_indices[na.omit(prop[active_rows, 8:9])]
    }

    if (all(prop[leaf_rows, 5] == top_EED[leaf_rows, 5])){
      accept <- TRUE
      break
    }
  }

  if (accept){
    return(prop)
  } else {
    stop("Attempts exceeded max_attempts")
  }
}
