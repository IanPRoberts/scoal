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
#' @param n_deme (optional) Number of demes
#'
#' @return Returns deme frequencies at each location
#'
#' @export

consensus_tree <- function(ED_sample, n_splits = 5, consensus_prob = 0.8, n_deme = NA){
  if (is.na(n_deme)){
    n_deme <- max(ED_sample[[1]][,5])
  }

  top_ED <- strip.history(ED_sample[[1]])
  root_node <- top_ED[is.na(top_ED[,2]),1]
  nodes <-top_ED[,1:2] #rbind(top_ED[is.na(top_ED[,3]), 1:2], top_ED[!is.na(top_ED[,4]), 1:2]) #Coalescent and leaf nodes with parents in second col
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

  phylo <- ed.to.phylo(consensus_ED)
  edge <- phylo$edge
  color.palette <- c("black", rainbow(n_deme))
  edge.color <- rep(NA,dim(edge)[1])
  for (i in 1 : dim(edge)[1]){
    edge.color[i] <- color.palette[phylo$node.deme[edge[i,2]] + 1]
  }

  plot(phylo, edge.color = edge.color, no.margin = TRUE, edge.width = 2, show.tip.label = FALSE)
  par(mar = 0.1 + c(5,4,4,2))

  return(consensus_ED)
}
