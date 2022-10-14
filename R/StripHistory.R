#' strip.history
#'
#' Strips the migration history from a phylogenetic tree
#'
#'
#' @param ED Extended data representation of a phylogeny including migration history
#'
#' @return An object of class \code{ED} with only coalescent nodes remaining and all located in deme 1
#'
#' @export

strip.history <- function(ED, node_indices = NA){
  if (length(node_indices) == 1){
    node_indices <- NodeIndicesC(ED)
  }

  coal_rows <- which(!is.na(ED[,4]))
  coal_nodes <- ED[coal_rows, 1]

  leaf_nodes <- ED[is.na(ED[,3]), 1]

  for (i in 1 : length(coal_rows)){
    row <- coal_rows[i]
    for (which_child in 3:4){
      child <- ED[row, which_child]
      while (!(child %in% c(coal_nodes, leaf_nodes))){
        child <- ED[node_indices[child], 3]
      }
      ED[row, which_child] <- child
      ED[node_indices[child], 2] <- ED[row, 1]
    }
  }

  ED[coal_rows, 5] <- 0
  ED <- ED[c(node_indices[leaf_nodes], coal_rows),]
  return(ED)
}
