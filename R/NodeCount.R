#' Counting of Migration and Coalescent Events
#'
#' Counts the number of migration events between each pair of demes and the number
#' of coalescence events occurring in each deme
#'
#' @param phylo An object of class \code{phylo} augmented with a deme for each node in the tree (with the deme of an edge being given by the end closest to the tips of the tree)
#' @param n.deme Number of demes in the structured coalescent process (optional argument)
#'
#' @return Returns a list containing a matrix M and a vector C. M has entries {m_{ij}} giving the number of migrations beginning in deme i and ending in deme j backwards in time. C has entries c_i giving the number of coalescence events occurring in deme i.
#'
#' @export


node.count <- function(phylo, n.deme = NA){
  demes <- unique(phylo$node.deme)

  if (is.na(n.deme)){
    n.deme <- length(demes)  ##max(demes)???
  }

  n <- length(phylo$tip.label)

  leaf.nodes <- 1:n
  migration.nodes <- phylo$edge[!phylo$edge[,1] %in% phylo$edge[duplicated(phylo$edge[,1]),1],1]
  n.migrations <- length(migration.nodes)
  coalescence.nodes <- ((n+1):(n.migrations+2*n-1))[!((n+1):(n.migrations+2*n-1) %in% migration.nodes)]

  M <- matrix(0, n.deme, n.deme)
  C <- numeric(n.deme)

  for (i in 1:n.deme){
    end.in.i <- migration.nodes[phylo$node.deme[migration.nodes] == i]
    origin.deme <- numeric(length(end.in.i))
    count <- 1
    for (j in end.in.i){
      origin.deme[count] <- phylo$node.deme[phylo$edge[which(phylo$edge[,1] == j),2]]
      count <- count + 1
    }
    M[,i] <- summary(factor(origin.deme, 1:3))
  }
  C <- summary(factor(phylo$node.deme[coalescence.nodes], 1:3))
  out <- list()
  out$C <- C
  out$M <- M
  return(out)
}

#' Counting of Migration and Coalescent Events
#'
#' Counts the number of migration events between each pair of demes and the number
#' of coalescence events occurring in each deme
#'
#' @param ED Extended data object; matrix with columns Node ID, parent, child 1, child 2, deme, node age
#' @param n.deme Number of demes in the structured coalescent process (optional argument)
#'
#' @return Returns a list containing a matrix M and a vector C. M has entries {m_{ij}} giving the number of migrations beginning in deme i and ending in deme j backwards in time. C has entries c_i giving the number of coalescence events occurring in deme i.
#'
#' @export

ed.node.count <- function(ED, n.deme, node.indices){
  coalescence.nodes <- ED[!is.na(ED[,4]), 1]
  migration.nodes <- ED[ is.na(ED[,4]) & (!is.na(ED[,3])), 1]

  coalescence.rows <- node.indices[coalescence.nodes]
  migration.rows <- node.indices[migration.nodes]

  m <- matrix(0, n.deme, n.deme)
  c <- rep(0, n.deme)

  mig.target <- ED[migration.rows, 5]
  mig.child <- ED[migration.rows, 3]
  mig.child.row <- node.indices[mig.child]
  mig.origin <- ED[mig.child.row, 5]

  for (j in 1 : length(migration.nodes)){
    m[mig.origin[j], mig.target[j]] <- m[mig.origin[j], mig.target[j]] + 1
  }

  for (i in 1 : n.deme){
    c[i] <- sum(ED[coalescence.rows, 5] == i)
  }
  return(list(c = c, m = m))
}
