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


#' Conversion Between Tree Objects
#'
#' @description  as.ED is a generic function converting a tree object into a tree of class "ED".
#' @description as.EED is a generic function converting a tree object into a tree of class "EED".
#' @description There are currently methods implemented to convert between objects of type "phylo", "ED", "EED"
#'
#' @param x object to be converted
#'
#' @return An object of class "ED"
#'
#' @export

as.ED <- function(x){
  if (identical(class(x), "ED")) return(x)
  UseMethod("as.ED")
}

#' @rdname as.ED
#' @export
as.EED <- function(x){
  if (identical(class(x), "EED")) return(x)
  UseMethod("as.EED")
}

#' @rdname as.ED
#' @export
as.ED.default <- function(x, ...){
  if (inherits(x, "ED")) return(x)
  stop('object does not inherit the class "ED": found no appropriate method to convert it')
}

#' @rdname as.ED
#' @export
as.EED.default <- function(x, ...){
  if (inherits(x, "EED")) return(x)
  stop('object does not inherit the class "EED": found no appropriate method to convert it')
}

#' @rdname as.ED
#' @export
as.ED.phylo <- function(x, ...){
  n_leaf <- length(x$tip.label)
  n_node <- x$Nnode + n_leaf
  n_edge <- nrow(x$edge)

  if (is.null(x$node.deme)) x$node.deme <- rep(0, n_node)

  ED <- matrix(NA, n_node, 6, dimnames = list(NULL, c("Node ID", "Parent", "Child 1", "Child 2", "Deme", "Node Age")))
  ED[,c(1,5,6)] <- c(1:n_node, x$node.deme, node.depth.edgelength(x))

  ED[x$edge[,2], 2] <- x$edge[,1] #Assign parent nodes - assume x$edge has edges in order (parent, child)

  for (i in 1 : n_node){
    if (is.na(ED[ED[i,2], 3])){
      ED[ED[i,2], 3] <- ED[i,1]
    } else {
      ED[ED[i,2], 4] <- ED[i,1]
    }
  }

  class(ED) <- 'ED'
  return(ED)
}

#' @rdname as.ED
#' @export
as.ED.EED <- function(x, ...){
  structure(x[,1:6], class = "ED")
}

#' @rdname as.ED
#' @export
as.EED.ED <- function(x, ...){
  node_indices <- NodeIndicesC(x)

  n_node <- nrow(x)

  EED <- matrix(NA, n_node, 9, dimnames = list(NULL, c("Node ID", "Parent", "Child 1", "Child 2", "Deme", "Node Age", "Parent coal", "Child coal 1", "Child coal 2")))
  EED[,1:6] <- x
  EED[,7] <- x[,2]

  for (i in 1 : n_node){
    if (!is.na(x[i,2])){ #Exclude root
      while (is.na(x[node_indices[EED[i, 7]], 4])){ #Loop until hit coalescent node
        EED[i, 7] <- x[node_indices[EED[i, 7]], 2]
      }
    }

    if (!is.na(x[i, 3])){ #Exclude leaves
      if (is.na(x[i, 4])){ #Migration only
        EED[i, 8] <- x[i, 3]
        while ((is.na(x[node_indices[EED[i,8]], 4])) && (!is.na(x[node_indices[EED[i,8]], 3]))) EED[i, 8] <- x[node_indices[EED[i, 8]], 3]
      } else { #Coalescence only
        for (y in 1 : 2){
          EED[i, 7 + y] <- x[i, 2 + y]
          while ((is.na(x[node_indices[EED[i, 7 + y]], 4])) && (!is.na(x[node_indices[EED[i, 7 + y]], 3]))) EED[x, 7 + y] <- x[node_indices[EED[i, 7 + y]], 3]
        }
      }
    }
  }
  class(EED) <- "EED"
  return(EED)
}

#' @rdname as.ED
#' @export
as.EED.phylo <- function(x, ...){
  as.EED(as.ED(x))
}

#' @rdname as.ED
#' @export
as.phylo.ED <- function(x, ...){
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
  class(phylo) <- c('str_phylo', 'phylo')
  phylo$edge <- edge
  phylo$edge.length <- edge.length
  phylo$tip.label <- 1:n.tips #ED[1:n.tips,1]
  phylo$Nnode <- n.nodes - n.tips
  phylo$node.deme <- ED[order(ED[,1]),5]  #Order supplies the ordering of the rows in ED to get node demes in correct order

  phylo <- ladderize(phylo)

  return(phylo)
}

#' @rdname as.ED
#' @export
as.phylo.EED <- function(x, ...){
  as.phylo(as.ED(x))
}
