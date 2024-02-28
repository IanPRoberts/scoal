#' Summary tree
#'
#' Generates a summary tree with pie charts at selected locations
#'
#'
#' @param ED_sample Sample of migration histories on same topology
#' @param locations Vector of locations to compute pie charts
#' @param n_deme (optional) Number of demes
#'
#' @return Returns deme frequencies at each location
#'
#' @export

summary_tree <- function(ED_sample, locations, n_deme = NA){
  if (is.na(n_deme)){
    n_deme <- max(ED_sample[[1]][,5])
  }

  top_ED <- strip.history(ED_sample[[1]]) #topology ED - all mig nodes & demes reset

  root_row <- which(is.na(top_ED[,2]))
  node_indices <- NodeIndicesC(top_ED)
  cum_lengths <- c(0,cumsum(top_ED[-root_row,6] - top_ED[node_indices[top_ED[-root_row, 2]],6]))

  below_coals <- sapply(1:length(locations), function(x) max(which(cum_lengths <= locations[x])))
  excess_height <- locations - cum_lengths[below_coals] #Height of locations above below coal node

  deme_freq <- matrix(0, n_deme, length(locations))

  for (j in 1 : length(ED_sample)){
    ED <- ED_sample[[j]]
    node_indices <- NodeIndicesC(ED)

    for (i in 1 : length(locations)){
      height <- 0
      row <- below_coals[i]
      deme <- ED[row, 5]

      while (height < excess_height[i]){
        parent_row <- node_indices[ED[row, 2]]
        deme <- ED[parent_row, 5]
        height <- ED[row, 6] - ED[parent_row, 6]
        row <- parent_row
      }

      deme_freq[deme, i] <- deme_freq[deme, i] + 1
    }
  }
  return(deme_freq)
}


#' Consensus tree
#'
#' Generates a consensus tree by discretising branches and colouring each section if the consensus probability is sufficiently large
#'
#'
#' @param ED_sample Sample of migration histories on same topology
#' @param n_splits Number of splits per branch used to discretise the tree
#' @param consensus_prob Consensus probability to be exceeded to allow colouring
#' @param plot Logical value of whether to plot consensus tree
#'
#' @return Returns consensus tree in ED format
#'
#' @export

consensus_tree <- function(ED_list, n_splits = 5, consensus_prob = 0.6, plot =  TRUE){
  n_deme <- max(ED_list[[1]][,5])
  n_trees <- length(ED_list)
  n_leaf <- sum(is.na(ED_list[[1]][,3]))

  #Strip migration events from tree
  topology <- ED_list[[1]]
  topology <- topology[(is.na(topology[,3])) | (!is.na(topology[,4])),]
  topology[,2:4] <- topology[,7:9]
  topology[!is.na(topology[,3]), 5] <- 0 #Set non-leaf demes to 0, i.e. undetermined

  NI_list <- lapply(ED_list, NodeIndicesC) # Node indices for each ED in ED_list

  #n_coal_nodes * n_trees matrix giving ordering of coalescent nodes for each tree
  ordered_coal_nodes <- sapply(ED_list, function(x){
    coal_rows <- !is.na(x[,4])
    x[coal_rows, 1][order(x[coal_rows, 6])] #x[...,1] gives coalescent node labels; [order(...)] sorts into ascending node age
  })

  max_label <- max(topology[,1])

  # Edges ending at coalescent node - index by (a.s. unique) event times
  for (coal_node_id in 2 : nrow(ordered_coal_nodes)){ #Root node is a.s. first node in coal_node_order
    observed_demes <- matrix(NA, n_trees, n_splits + 1)

    coal_node_row <- NI_list[[1]][ordered_coal_nodes[coal_node_id, 1]]
    edge_length <- topology[coal_node_row, 6] - topology[NI_list[[1]][topology[coal_node_row, 7]], 6] #Branch length of current branch only
    split_width <- edge_length / (n_splits + 1)
    observation_times <- topology[coal_node_row, 6] - split_width * 1:n_splits #Age of splits along branch

    for (tree_id in n_trees : 1){ #Iterate in reverse order to leave final setup at tree_id = 1 to be same as topology
      ED <- ED_list[[tree_id]]
      NI <- NI_list[[tree_id]]

      coal_node_label <- ordered_coal_nodes[coal_node_id, tree_id]
      coal_node_row <- NI[coal_node_label]

      observed_demes[tree_id, 1] <- ED[coal_node_row, 5] #Observed demes at current coalescent node

      #Find all nodes on branch (including coal_node and parent coal of coal_node)
      coal_parent <- ED[coal_node_row, 7]
      branch_nodes <- unique(c(coal_parent, #Possibly add coal_parent twice if ED[coal_parent, 8] == coal_node_label
                               ED[sapply(ED[,8], identical, coal_node_label), 1], #sapply(..., identical) avoids 'NA==...' returning NA
                               coal_node_label)) #Add coal_node_label at end
      branch_rows <- NI[branch_nodes]
      branch_rows <- branch_rows[order(ED[branch_rows, 6])] #Ensure branch nodes are stored in ascending age order

      observed_demes[tree_id, 1 + 1:n_splits] <- ED[branch_rows, 5][1 + findInterval(observation_times, ED[branch_rows, 6])]
    }

    consensus_deme <- numeric(n_splits + 1)
    consensus_freq <- consensus_prob * n_trees
    for (split_id in 1 : (n_splits + 1)){
      tab <- table(observed_demes[,split_id])

      possible_demes <- names(tab[tab > consensus_freq])
      if (length(possible_demes) == 0){
        consensus_deme[split_id] <- 0
      } else {
        consensus_deme[split_id] <- as.numeric(sample(possible_demes, 1)) #Sample uniformly in case of tie over multiple options (only possible if consensus_prob <= 0.5)
      }
    }

    # Add nodes to topology at split points
    topology <- rbind(topology,
                      cbind(max_label + 1 : n_splits, #Node ID
                            c(max_label + 2 : n_splits, coal_parent), #Parent
                            c(coal_node_label, max_label + 2 : n_splits - 1), #Child 1
                            NA, #Child 2
                            consensus_deme[-1], #Deme
                            observation_times, #Node Age
                            coal_parent, #Parent coal
                            coal_node_label, #Child coal 1
                            NA)) #Child coal 2

    topology[coal_node_row, c(2, 5)] <- c(max_label + 1, consensus_deme[1]) #Update parent of coal_node

    #Update child of parent_coal_node
    coal_parent_row <- NI[coal_parent]
    which_child <- which(topology[coal_parent_row, 8:9] == coal_node_label)
    topology[NI[coal_parent], 2 + which_child] <- max_label + n_splits

    max_label <- max_label + n_splits
  }

  # Edges ending at migration node - index by tip label
  tip_labels <- topology[is.na(topology[,3]), 1]

  for (tip_id in 1 : length(tip_labels)){
    observed_demes <- matrix(NA, n_trees, n_splits + 1)

    tip_row <- NI_list[[1]][tip_labels[tip_id]]
    edge_length <- topology[tip_row, 6] - topology[NI_list[[1]][topology[tip_row, 7]], 6] #Branch length above tip
    split_width <- edge_length / (n_splits + 1)
    observation_times <- topology[tip_row, 6] - split_width * 1 : n_splits #Age of splits along branch

    for (tree_id in n_trees : 1){
      ED <- ED_list[[tree_id]]
      NI <- NI_list[[tree_id]]

      tip_node_label <- tip_labels[tip_id]
      tip_node_row <- NI[tip_node_label]

      observed_demes[tree_id, 1] <- ED[tip_node_row, 5] #Observed demes at current tip

      #Find all nodes on branch (including coal_node and parent coal of coal_node)
      coal_parent <- ED[tip_node_row, 7]
      branch_nodes <- unique(c(coal_parent, #Possibly add coal_parent twice if ED[coal_parent, 8] == tip_node_label
                               ED[sapply(ED[,8], identical, tip_node_label), 1], #sapply(..., identical) avoids 'NA==...' returning NA
                               tip_node_label)) #Add tip_node_label at end
      branch_rows <- NI[branch_nodes]
      branch_rows <- branch_rows[order(ED[branch_rows, 6])] #Ensure branch nodes are stored in ascending age order

      observed_demes[tree_id, 1 + 1:n_splits] <- ED[branch_rows, 5][1 + findInterval(observation_times, ED[branch_rows, 6])]
    }

    consensus_deme <- numeric(n_splits + 1)
    consensus_freq <- consensus_prob * n_trees
    for (split_id in 1 : (n_splits + 1)){
      tab <- table(observed_demes[,split_id])

      possible_demes <- names(tab[tab >= consensus_freq])
      if (length(possible_demes) == 0){
        consensus_deme[split_id] <- 0
      } else {
        consensus_deme[split_id] <- as.numeric(sample(possible_demes, 1)) #Sample uniformly in case of tie over multiple options (only possible if consensus_prob <= 0.5)
      }
    }

    # Add nodes to topology at split points
    topology <- rbind(topology,
                      cbind(max_label + 1 : n_splits, #Node ID
                            c(max_label + 2 : n_splits, coal_parent), #Parent
                            c(tip_node_label, max_label + 2 : n_splits - 1), #Child 1
                            NA, #Child 2
                            consensus_deme[-1], #Deme
                            observation_times, #Node Age
                            coal_parent, #Parent coal
                            tip_node_label, #Child coal 1
                            NA)) #Child coal 2

    topology[tip_node_row, c(2, 5)] <- c(max_label + 1, consensus_deme[1]) #Update parent of coal_node

    #Update child of parent_coal_node
    coal_parent_row <- NI[coal_parent]
    which_child <- which(topology[coal_parent_row, 8:9] == coal_node_label)
    topology[NI[coal_parent], 2 + which_child] <- max_label + n_splits

    max_label <- max_label + n_splits
  }

  # Remove unnecessary self-migration events
  NI <- NodeIndicesC(topology)
  rm_rows <- numeric(0)
  for (row_id in 1 : nrow(topology)){
    if ((is.na(topology[row_id, 4])) & (!is.na(topology[row_id, 3]))){ #If migration event
      parent_row <- NI[topology[row_id, 2]]
      parent_deme <- topology[parent_row, 5]

      if (parent_deme == topology[row_id, 5]){ #If current deme matches parent deme
        rm_rows <- append(rm_rows, row_id) #Add current row to be removed
        which_child <- which(topology[parent_row, 3:4] == topology[row_id, 1])
        topology[parent_row, 2 + which_child] <- topology[row_id, 3] #Child(parent_row) is now child(current_row)

        child_row <- NI[topology[row_id, 3]]
        topology[child_row, 2] <- topology[row_id, 2] #Parent of child(current_row) is now parent(current_row)
      }
    }
  }

  topology <- topology[-rm_rows,]


  if (plot){
    structured.plot(topology)
  }

  return(ED=topology)
}

consensus_tree_old <- function(ED_sample, n_splits = 5, consensus_prob = 0.8, n_deme = NA, plot = TRUE){
  if (is.na(n_deme)){
    n_deme <- max(ED_sample[[1]][,5])
  }

  top_ED <- strip.history(ED_sample[[1]])
  root_node <- top_ED[is.na(top_ED[,2]),1]
  nodes <-top_ED[,1:2]
  nodes <- nodes[nodes[,1] != root_node,]

  deme_freq <- matrix(0, n_splits * dim(nodes)[1], n_deme)
  node_freq <- matrix(0, dim(nodes)[1], n_deme)

  pb <- txtProgressBar(min = 0, max = length(ED_sample), initial = 0, style = 3)

  for (i in 1 : length(ED_sample)){
    ED <- ED_sample[[i]]
    node_indices <- NodeIndicesC(ED)
    rows <- matrix(node_indices[nodes], ncol = 2)
    row_heights <- matrix(ED[rows, 6], ncol = 2)

    split_heights <- matrix(NA, dim(nodes)[1], n_splits)

    for (j in 1 : dim(nodes)[1]){
      node_freq[j, ED[rows[j,1], 5]] <- node_freq[j, ED[rows[j,1], 5]] + 1
      split_heights[j,] <- seq(row_heights[j,1], row_heights[j,2], length.out = n_splits + 2)[1 : n_splits + 1]
      height <- row_heights[j,1]
      node_row <- rows[j,1]

      for (k in 1 : n_splits){
        while (height > split_heights[j,k]){
          node_row <- node_indices[ED[node_row, 2]]
          height <- ED[node_row, 6]
        }
        deme_freq[(j - 1) * n_splits + k, ED[node_row, 5]] <- deme_freq[(j - 1) * n_splits + k, ED[node_row, 5]] + 1
      }
    }

    setTxtProgressBar(pb, i)
  }

  close(pb)

  consensus_freq <- consensus_prob * length(ED_sample)
  split_deme <- sapply(1 : dim(deme_freq)[1],
                       function(x){
                         if (any(deme_freq[x,] > consensus_freq)){
                           return(which(deme_freq[x,] == max(deme_freq[x,])))
                         } else{
                           return(0)
                         }
                       })

  node_deme <- sapply(1 : dim(node_freq)[1],
                      function(x){
                        if (any(node_freq[x,] > consensus_freq)){
                          return(which(node_freq[x,] == max(node_freq[x,])))
                        } else{
                          return(0)
                        }
                      })

  max_filled_row <- dim(top_ED)[1]
  consensus_ED <- matrix(0, max_filled_row + length(split_deme), 6)
  consensus_ED[1 : max_filled_row,] <- top_ED
  consensus_ED[(1 : max_filled_row)[-which(top_ED[,1] == root_node)], 5] <- node_deme
  node_indices <- NodeIndicesC(top_ED)
  rows <- matrix(node_indices[nodes], ncol = 2)
  max_label <- max(top_ED[,1])

  for (j in 1 : dim(rows)[1]){
    consensus_ED[rows[j,1], 2] <- max_label + 1
    for (k in 1 : n_splits){
      max_filled_row <- max_filled_row + 1
      max_label <- max_label + 1

      if (k == n_splits){
        parent <- top_ED[rows[j, 1], 2]
      } else {
        parent <- max_label + 1
      }

      consensus_ED[max_filled_row,] <- c(max_label, parent, max_label - 1, NA, split_deme[(j-1) * n_splits + k], split_heights[j,k])

    }
  }

  if (plot){
    structured.plot(consensus_ED)
  }

  return(list(ED = consensus_ED, deme_freq = deme_freq))
}
