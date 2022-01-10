#' Migration Rates Update
#'
#' Performs an update of the migration rates matrix via a Gibbs move with Gamma prior and Gamma proposal distribution
#'
#' @param ED Extended data object; matrix with columns Node ID, parent, child 1, child 2, deme, node age
#' @param migration.matrix Current migration rates matrix for the structured coalescent process
#' @param n.deme Number of demes in the structured coalescent process (optional argument)
#' @param alpha Gamma shape parameter for the prior on migration rates
#' @param beta Gamma rate parameter for the prior on migration rates
#'
#' @return Updated migration rates matrix conditioned on the current migration history
#'
#' @export

mig.rate.update <- function(ED, migration.matrix, n.deme = NA, alpha = 1, beta = 10){
  m <- ed.node.count(ED, n.deme)$m

  observed.demes <- unique(ED[,5])
  if (is.na(n.deme)){
    n.deme <- length(observed.demes)
  }
  proposal.matrix <- matrix(0, n.deme, n.deme)

  node.heights <- ED[,6]
  event.times <- sort(unique(node.heights))  #Times of events, root at time=0
  time.increments <- diff(event.times)

  root.node <- ED[is.na(ED[,2]), 1]
  coalescence.nodes <- ED[!is.na(ED[,4]),1]
  migration.nodes <- ED[ is.na(ED[,4]) & (!is.na(ED[,3])) ,1]
  leaf.nodes <- ED[(is.na(ED[,3])) & (is.na(ED[,4])), 1]

  k <- matrix(0, nrow = length(event.times) - 1, ncol = n.deme)
  k[1,ED[which(ED[,1] == root.node),5]] <- 2
  for (i in 2 : (length(event.times) - 1)){
    current.rows <- which(ED[,6] == event.times[i])
    k[i,] <- k[i-1,]
    if (length(current.rows) > 1){ #Multiple leaves added simultaneously
      for (j in current.rows){
        k[i, ED[j, 5]] <- k[i, ED[j, 5]] + 1
      }
    } else{
      if (current.rows %in% migration.nodes){ #Migration event
        k[i, ED[current.rows, 5]] <- k[i, ED[current.rows, 5]] - 1
        current.child <- ED[current.rows, 3]
        k[i, ED[current.child, 5]] <- k[i, ED[current.child, 5]] + 1
      } else if (current.rows %in% coalescence.nodes){ #Coalescence event
        k[i, ED[current.rows, 5]] <- k[i, ED[current.rows, 5]] + 1
      } else{ #Single leaf added
        k[i, ED[current.rows, 5]] <- k[i, ED[current.rows, 5]] - 1
      }
    }
  }

  rate.constants <- t(k) %*% time.increments

  for (i in 1:n.deme){
    for (j in (1:n.deme)[-i]){
      proposal.matrix[i,j] <- rgamma(1, shape = alpha + m[i,j], scale = 1 / (beta + rate.constants[i]))
    }
  }
  return(proposal.matrix)
}
