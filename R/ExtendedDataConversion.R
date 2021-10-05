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

ed.to.phylo <- function(ED){
  n.nodes <- dim(ED)[1]
  n.tips <- sum(is.na(ED[,3]))

  #Check no gaps in node labelling scheme to allow phylo object to be generated
  if (max(ED[,1]) > n.nodes){
    missing.labels <- (1:n.nodes)[! (1:n.nodes) %in% ED[,1]]
    extra.labels <- as.vector(ED[ED[,1] > n.nodes,1])  #as.vector needed in case only 1 extra label has appeared

    node.label.mat <- ED[,1:4]
    count <- 1
    for (i in extra.labels){
      node.label.mat[node.label.mat == i] <- missing.labels[count]
      count <- count + 1
    }
    ED[, 1:4] <- node.label.mat
  }
  edge.list <- list()
  edge.length <- numeric(0)
  count <- 1
  for (i in (1:n.nodes)[-(n.tips + 1)]){
    row <- which(ED[,1] == i)
    parent.row <- which(ED[,1] == ED[row, 2])
    edge.list[[count]] <- c(ED[row, 2], i)
    edge.length[count] <- ED[row, 6] - ED[parent.row, 6]
    count <- count + 1
  }
  edge <- do.call(rbind,edge.list)  #construct edge matrix from edge.list

  phylo <- list()
  class(phylo) <- 'phylo'
  phylo$edge <- edge
  phylo$edge.length <- edge.length
  phylo$tip.label <- 1:n.tips #ED[1:n.tips,1]
  phylo$Nnode <- n.nodes - n.tips
  phylo$node.deme <- ED[order(ED[,1]),5]  #Order supplies the ordering of the rows in ED to get node demes in correct order

  return(phylo)
}
