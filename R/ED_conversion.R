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

  ED <- matrix(NA, n.nodes, 6, dimnames = list(NULL, c("Node ID", "Parent", "Child 1", "Child 2", "Deme", "Node Age")))
  ED[,c(1,5,6)] <- c(nodes, phylo$node.deme, node.depth.edgelength(phylo))

  for (i in (1 : n.nodes)[-(length(phylo$tip.label)+1)]){
    ED[i,2] <- phylo$edge[which(phylo$edge[,2] == ED[i,1]),1]
  }

  for (i in (length(phylo$tip.label)+1):n.nodes){
    children <- phylo$edge[which(phylo$edge[,1] == ED[i,1]),2]
    if (length(children) == 1){
      ED[i,3] <- children
    } else {
      ED[i,c(3,4)] <- children
    }
  }

  #Label root node n+1, coalescence nodes (n+2):(2n-1), migration nodes (2n-1):(2n+M-1)

  migration.nodes <- ED[(!is.na(ED[,3])) & (is.na(ED[,4])),1]
  coalescence.nodes <- ED[!is.na(ED[,4]),1]
  n <- length(ED[(is.na(ED[,3])) & (is.na(ED[,4])),1])
  spare.label <- max(ED[,1]) + 1

  if (! all((n+1):(2*n-1) %in% coalescence.nodes)){
    missing.labels <-((n+1):(2*n-1))[! ((n+1):(2*n-1) %in% coalescence.nodes)]
    extra.labels <- coalescence.nodes[! (coalescence.nodes %in% (n+1):(2*n-1))]

    node.label.mat <- ED[,1:4]

    for (i in 1 : length(missing.labels)){
      mig.row <- which(ED[,1] == missing.labels[i])
      coal.row <- which(ED[,1] == extra.labels[i])

      node.label.mat[node.label.mat == missing.labels[i]] <- spare.label
      node.label.mat[node.label.mat == extra.labels[i]] <- missing.labels[i]
      node.label.mat[node.label.mat == spare.label] <- extra.labels[i]
    }

    ED[, 1:4] <- node.label.mat
    ED <- ED[order(ED[,1]),]
  }

  class(ED) <- 'ED'
  return(ED)
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

  phylo <- ladderize(phylo)

  return(phylo)
}

#' ED to EED Conversion
#'
#' Converts an EED object for a structured coalescent tree into an EED, adding
#' the parent and children coalescent nodes for each node
#'
#' @param ED phylo object augmented with the deme of each node
#' @param node_indices (optional) vector relating row labels with row number
#'
#' @return A matrix with each row containing the parent node, child nodes, deme,
#' age, parent coalescent node and child coalescent nodes of
#' each node in the tree.
#'
#' @export

ed.to.eed <- function(ED, node_indices = NA){
  if (is.na(node_indices)){
    node_indices <- NodeIndicesC(ED)
  }

  EED <- matrix(NA, dim(ED)[1], 9, dimnames = list(NULL, c("Node ID", "Parent", "Child 1", "Child 2", "Deme", "Node Age", "Parent coal", "Child coal 1", "Child coal 2")))
  EED[,1:6] <- ED

  for (x in 1 : dim(ED)[1]){
    if (!is.na(ED[x,2])){ #Exclude root
      EED[x, 7] <- ED[x,2] #Identify parent coalescent node
      while (is.na(ED[node_indices[EED[x, 7]], 4])) EED[x, 7] <- ED[node_indices[EED[x, 7]], 2]
    }

    if (!is.na(ED[x, 3])){ #Exclude leaves
      if (is.na(ED[x, 4])){ #Migration only
        EED[x, 8] <- ED[x, 3]
        while ((is.na(ED[node_indices[EED[x,8]], 4])) && (!is.na(ED[node_indices[EED[x,8]], 3]))) EED[x, 8] <- ED[node_indices[EED[x, 8]], 3]
      } else { #Coalescence only
        for (y in 1 : 2){
          EED[x, 7 + y] <- ED[x, 2 + y]
          while ((is.na(ED[node_indices[EED[x, 7 + y]], 4])) && (!is.na(ED[node_indices[EED[x, 7 + y]], 3]))) EED[x, 7 + y] <- ED[node_indices[EED[x, 7 + y]], 3]
        }
      }
    }
  }
  return(EED)
}
