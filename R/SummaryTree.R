#' Coalescent Node Pie Charts
#'
#' Tabulates the frequency of each deme at each leaf and coalescent node of a phylogenetic tree
#'
#'
#' @param ED_sample Sample of migration histories on same topology (ED with 9 columns)
#' @param plot Logical value whether to plot tree with superimposed pie charts
#' @param plot_ED (optional) ED to superimpose with pie charts (clear migration history if left blank)
#'
#' @return Returns deme frequencies at each location
#'
#' @export

coalescent_node_pie_charts <- function(ED_list, plot = TRUE, plot_ED = matrix(NA, 0, 9), cex = 0.5){
  n_trees <- length(ED_list)
  n_deme <- max(ED_list[[1]][,5])

  #Strip migration events from tree
  topology <- ED_list[[1]]
  topology <- topology[(is.na(topology[,3])) | (!is.na(topology[,4])),]
  topology[,2:4] <- topology[,7:9]
  topology[!is.na(topology[,3]), 5] <- 0

  #Construct deme frequency at each coalescent node in ascending node age
  deme_freq <- matrix(0, nrow(topology), n_deme)

  for (tree_id in n_trees : 1){ #Loop in reverse tree order to leave node_order as order(topology[,6]) after final iteration
    ED <- ED_list[[tree_id]]
    ED <- ED[(is.na(ED[,3])) | (!is.na(ED[,4])),]
    ED[,2:4] <- ED[,7:9]

    node_order <- order(ED[,6]) #Store entries of deme_freq in ascending age (root = 0, newest leaf = max(ED[,6]))

    for (row_id in 1 : nrow(ED)){ #Loop robust only for non-simultaneous events, i.e. coalescent nodes
      deme_freq[row_id, ED[node_order[row_id], 5]] <- deme_freq[row_id, ED[node_order[row_id], 5]] + 1
    }
  }

  #Remove leaf deme frequencies from deme_freq (should be a.s. one deme throughout run)
  deme_freq <- deme_freq[!is.na(topology[node_order, 4]),]

  if (nrow(plot_ED) == 0){
    plot_ED <- topology
    plot_ED[,5] <- 0
  }

  plot_ED_coal_nodes <- !is.na(plot_ED[,4])
  plot_ED_coal_node_order <- order(plot_ED[plot_ED_coal_nodes, 6])


  rownames(deme_freq) <- plot_ED[plot_ED_coal_nodes, 1][plot_ED_coal_node_order]

  if (plot){
    structured.plot(plot_ED)
    pie_plot <- rowSums(deme_freq == 0) < n_deme - 1 #Logical on whether 100% same deme observed (in which case no pie chart plotted!)
    nodelabels(node = as.numeric(rownames(deme_freq))[pie_plot],
               pie = deme_freq[pie_plot,]/rowSums(deme_freq[pie_plot,]),
               cex = cex)
  }

  return(list(ED = topology, node_freq = deme_freq))
}


#' Consensus tree
#'
#' Generates an approximate consensus tree by discretising branches and colouring each section if the consensus probability is sufficiently large
#'
#' @inheritParams coalescent_node_pie_charts
#'
#' @param n_splits Number of splits per branch used to discretise the tree
#' @param consensus_prob Consensus probability to be exceeded to allow colouring
#'
#' @return Returns consensus tree in ED format
#'
#' @export

approximate_consensus_tree <- function(ED_list, n_splits = 5, consensus_prob = 0.5, plot =  TRUE){
  n_deme <- max(ED_list[[1]][,5])
  n_trees <- length(ED_list)
  n_leaf <- sum(is.na(ED_list[[1]][,3]))

  #Strip migration events from tree
  topology <- ED_list[[1]]
  topology <- topology[(is.na(topology[,3])) | (!is.na(topology[,4])),]
  topology[,2:4] <- topology[,7:9]
  topology[!is.na(topology[,3]), 5] <- 0 #Set non-leaf demes to 0, i.e. undetermined

  NI_list <- lapply(ED_list, NodeIndices) # Node indices for each ED in ED_list

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
  NI <- NodeIndices(topology)
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
