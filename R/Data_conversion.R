################# To do
# 1. [DONE] overload plot.strphylo to structured.plot
# [No longer exporting structured.plot in NAMESPACE but maintain for now until remainder of code fixed]
# 2. consolidate references to EED into ED and update ED references to EED only (e.g. Ewing topology moves need cols 7:9 to be added)
# [PARTIAL - Ewing MCMC moves updated]
# [REMAINING - DTA_sampling.R, Initial_tree.R]
# (consolidate to a single extended data structure rather than ED and EED)
# 3. MCMC tree output as beast .trees file (nexus with type annotation at each node)
# [PARTIAL - written write.beast.multiple to convert list of strphylo info into .trees file]
# [PARTIAL - ~/Desktop/ewing_mcmc.R and ~/Desktop/dta_mcmc.R save to .trees using cat() but inefficient implementation ED -> strphylo -> treeio::treedata -> treeio::write.beast.newick]
#################

#' Conversion to Structured Phylo
#'
#' Conversion to structured phylo objects
#'
#' @param x object to be converted
#'
#' @return An object of class "strphylo"
#'
#' @export

as.strphylo <- function(x, ...){
  if (identical(class(x), c("strphylo", "phylo"))) return(x)
  UseMethod("as.strphylo")
}

#' @rdname as.strphylo
#' @export

as.strphylo.default <- function(x, ...){
  if (inherits(x, "strphylo")) return(x)
  stop('object does not inherit the class "strphylo": found no appropriate method to convert it')
}

#' @rdname as.strphylo
#' @export

as.strphylo.phylo <- function(x, node.deme, ...){
  x$node.deme <- node.deme
  structure(x, class = c("strphylo", "phylo"))
}


###### Internal function to convert ED or EED to strphylo
#' @rdname as.strphylo
#'@export
as.strphylo.matrix <- function(x, ...){
  n_node <- nrow(x)
  n_tip <- sum(is.na(x[,3]))

  #Check no gaps in node labelling scheme to allow phylo object to be generated
  if (max(x[,1]) > n_node){
    missing_labs <- (1:n_node)[! (1:n_node) %in% x[,1]]
    extra_labs <- as.vector(x[x[,1] > n_node,1])  #as.vector needed in case only 1 extra label has appeared

    node_labels <- x[,1:4]
    count <- 1
    for (i in extra_labs){
      node_labels[node_labels == i] <- missing.labels[count]
      count <- count + 1
    }
    x[, 1:4] <- node_labels
  }
  edge.list <- list()
  edge.length <- numeric(0)
  count <- 1
  for (i in (1:n_node)[-(n_tip + 1)]){
    row <- which(x[,1] == i)
    parent.row <- which(x[,1] == x[row, 2])
    edge.list[[count]] <- c(x[row, 2], i)
    edge.length[count] <- x[row, 6] - x[parent.row, 6]
    count <- count + 1
  }
  edge <- do.call(rbind,edge.list)  #construct edge matrix from edge.list

  phylo <- list()
  class(phylo) <- c('strphylo', 'phylo')
  phylo$edge <- edge
  phylo$edge.length <- edge.length
  phylo$tip.label <- 1:n_tip
  phylo$Nnode <- n_node - n_tip
  phylo$node.deme <- x[order(x[,1]),5]  #Order supplies the ordering of the rows in x to get node demes in correct order

  phylo <- ladderize(phylo)
  return(phylo)
}
