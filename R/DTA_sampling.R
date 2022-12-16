#' DTA Sampling
#'
#' Samples N migration histories under the DTA model deme fixed at all leaves
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

DTA_sampling <- function(ED, mig_mat, time_scale = 1, N = 1, parallel = TRUE, mc.cores = NA){
  if (parallel){
    if (is.na(mc.cores)){
      options(mc.cores = parallel::detectCores() /2)
    }
  } else {
    options(mc.cores = 1)
  }

  fit_rates <- mig_mat
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
#'
#' @return Proposal ED
#'
#' @export

local_DTA <- function(ED, mig_mat, time_scale){
  #############################################
  #   Select subtree & sample internal deme   #
  #############################################

  n_deme <- dim(mig_mat)[1]
  fit_rates <- mig_mat * time_scale
  diag(fit_rates) <- - rowSums(fit_rates)
  top_ED <- strip.history(ED)
  node_indices <- NodeIndicesC(top_ED)
  ED_node_indices <- NodeIndicesC(ED)

  leaf_nodes <- top_ED[is.na(top_ED[,3]),1]
  coal_nodes <- top_ED[-node_indices[leaf_nodes],1]

  coal_rows <- node_indices[coal_nodes]

  selected_row <- sample(coal_rows, 1)
  selected_node <- top_ED[selected_row, 1]
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
      active_row <- prop_ED[active_row, 2]
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
      active_row <- prop_ED[active_row, 2]
    }

    if (n_mig > 0){ #Add new migrations
      parent_node <- prop_ED[prop_parent_row, 1]
      which.child <- which(prop_ED[prop_parent_row, 3:4] %in% c(selected_node, prop_node_indices[rm_rows]))
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
      which.child <- which(prop_ED[prop_parent_row, 3:4] %in% prop_node_indices[rm_rows])
      prop_ED[prop_parent_row, 2 + which.child] <- selected_node
    }
  }
  if (length(rm_rows) > 0) prop_ED <- prop_ED[-rm_rows,]
  return(prop_ED)
}
