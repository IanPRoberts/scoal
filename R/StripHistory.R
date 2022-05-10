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

strip.history <- function(ED){
  coalescence.nodes <- ED[(!is.na(ED[,3])) & (!is.na(ED[,4])), 1]
  leaf.nodes <- ED[(is.na(ED[,3])) & (is.na(ED[,4])), 1]
  root.node <- ED[is.na(ED[,2]), 1]

  for (node in c(leaf.nodes, coalescence.nodes)[-root.node]){
    node.parent <- ED[node, 2]
    while (! node.parent %in% coalescence.nodes){
      node.parent <- ED[node.parent, 2]
      ED[node, 2] <- node.parent
    }
  }

  migration.nodes <- ED[! ED[,1] %in% c(leaf.nodes, coalescence.nodes), 1]
  ED <- ED[! ED[,1] %in% migration.nodes,]
  ED[! ED[,1] %in% leaf.nodes, 5] <- 1

  return(ED)
}
