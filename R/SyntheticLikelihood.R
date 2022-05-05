#' Calculates synthetic likelihood for a coalescent tree with migration history
#'
#' Computes a synthetic likelihood for a coalescent tree with migration history
#' with a Poisson marginal distribution on the number of migrations
#'
#' @param ED Extended data representation of a structured coalescent process
#' @param rate Poisson rate for synthetic likelihood
#' @param effective.pop Effective population sizes
#' @param migration.matrix Migration matrix
#' @param node.indices
#'
#' @return A vector (\code{log-likelihood}, \code{likelihood}) giving the log-likelihood and likelihood of \code{phylo}
#'
#' @export

synth.likelihood <- function(ED, effective.pop, rate, migration.matrix, node.indices){
  root.row <- which(is.na(ED[,2]))
  n.deme <- length(effective.pop)
  tree.length <- 0

  for (i in (1 : dim(ED)[1])[-root.row]){
    node.parent <- ED[i,2]
    parent.row <- node.indices[node.parent]
    tree.length <- tree.length + ED[i,6] - ED[parent.row,6]
  }

  M <- dim(ED)[1] - root.row
  like <- M * (log(rate) - log(tree.length) - log(n.deme - 1)) - log(n.deme) - rate

  return(list(log.likelihood = like, likelihood = exp(like)))
}
