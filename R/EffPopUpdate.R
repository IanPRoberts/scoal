#' Effective Population Size Update
#'
#' Performs an update of the effective population sizes via a Gibbs move with Inverse-gamma prior and Inverse-gamma proposal distribution
#'
#' @param ED Extended data object; matrix with columns Node ID, parent, child 1, child 2, deme, node age
#' @param effective.population Current migration rates matrix for the structured coalescent process
#' @param n.deme Number of demes in the structured coalescent process (optional argument)
#' @param alpha Inverse-gamma shape parameter for the prior on effective population sizes
#' @param beta Inverse-gamma rate parameter for the prior on effective population sizes
#'
#' @return Updated vector of effective population sizes conditioned on the current migration history
#'
#' @export

eff.pop.update <- function(ED, effective.population, n.deme, alpha = 1, beta = 1){
  c <- ed.node.count(ED, n.deme)$c

  observed.demes <- unique(ED[,5])
  if (is.na(n.deme)){
    n.deme <- length(observed.demes)
  }
  proposal.eff.pop <- rep(0, n.deme)

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
        k[i, ED[j, 5]] <- k[i, ED[j, 5]] - 1
      }
    } else{
      if (current.rows %in% migration.nodes){ #Migration event
        k[i, ED[current.rows, 5]] <- k[i, ED[current.rows, 5]] - 1
        current.child <- ED[current.rows, 3]
        current.child.row <- which(ED[,1] == current.child)
        k[i, ED[current.child.row, 5]] <- k[i, ED[current.child.row, 5]] + 1
      } else if (current.rows %in% coalescence.nodes){ #Coalescence event
        k[i, ED[current.rows, 5]] <- k[i, ED[current.rows, 5]] + 1
      } else{ #Single leaf added
        k[i, ED[current.rows, 5]] <- k[i, ED[current.rows, 5]] - 1
      }
    }
  }

  rate.constants <- t(k * (k-1) / 2) %*% time.increments

  for (i in 1:n.deme){
    proposal.eff.pop[i] <- 1/rgamma(1, shape = alpha + c[i], scale = 1 / (beta + rate.constants[i]))  #Proposals are inverse-gamma
  }
  return(proposal.eff.pop)
}
