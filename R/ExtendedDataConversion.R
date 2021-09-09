#' Phylo Object to Extended Data
#'
#' Converts a phylo object for a structured coalescent tree into an extended
#' data format containing the parent and children of each node alongside each
#' node's deme
#'
#' @param phylo phylo object augmented with the deme of each node
#'
#' @return A matrix with each row containing the parent, children, deme and age of
#' each node in the tree.
#'
#' @export

phylo.to.ed <- function(phylo){
  nodes <- sort(unique(as.vector(phylo$edge)))
  n.nodes <- length(nodes)
  n.edges <- dim(phylo$edge)[1]

  Out <- matrix(NA, n.nodes, 6, dimnames = list(NULL, c("Node ID", "Parent", "Child 1", "Child 2", "Deme", "Node Age")))
  Out[,c(1,5,6)] <- c(nodes, phylo$node.deme, node.depth.edgelength(phylo))

  for (i in (1 : n.nodes)[-(length(phylo$tip.label)+1)]){
    Out[i,2] <- phylo$edge[which(phylo$edge[,2] == Out[i,1]),1]
  }

  for (i in (length(phylo$tip.label)+1):n.nodes){
    children <- phylo$edge[which(phylo$edge[,1] == Out[i,1]),2]
    if (length(children) == 1){
      Out[i,3] <- children
    } else {
      Out[i,c(3,4)] <- children
    }
  }
  return(Out)
}

#' Extended Data to Phylo Object
#'
#' Converts an extended data object into a phylo object augmented with the deme
#' of each node
#'
#' @param ED Extended data object to convert; matrix with columns Node ID, parent, child 1, child 2, deme, node age
#'
#' @return A phylo object augmented with node demes.
#'
#' @export
