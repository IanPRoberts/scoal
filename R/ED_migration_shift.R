#' Migration Shift MCMC Proposal
#'
#' Performs a migration shift MCMC proposal on a structured coalescent
#' migration history. A migration event is selected uniformly at random
#' and the time of the event is updated by resampling the event time between
#' the times of the preceding and following events.
#'
#' @param ED Extended data object; matrix with columns Node ID, parent, child 1, child 2, deme, node age
#' @param n.deme Number of distinct demes in the population
#'
#' @return Updated extended data object with the proposal from the migration shift proposal
#'
#' @export

ed.mig.shift.1 <- function(ED, n.deme){
  migration.nodes <- ED[ is.na(ED[,4]) & (!is.na(ED[,3])) ,1]
  selected.node <- sample.vector(migration.nodes, 1)
  selected.row <- which(ED[,1] == selected.node)

  child.node <- ED[selected.row, 3]
  child.row <- which(ED[,1] == child.node)
  child.time <- ED[child.row, 6]

  parent.node <- ED[selected.row, 2]
  parent.row <- which(ED[,1] == parent.node)
  parent.time <- ED[parent.row, 6]

  new.time <- runif(1, parent.time, child.time)
  ED[selected.row, 6] <- new.time

  return(list(ED = ED, prop.ratio = 1))
}


#' Migration Shift MCMC Proposal
#'
#' Performs a migration shift MCMC proposal on a structured coalescent
#' migration history. A migration event is selected uniformly at random
#' and the time of the event is updated by perturbing the event time by a
#' N(0, var) random variable with reflecting boundaries at the times of the
#' preceding and following events.
#'
#' @param ED Extended data object; matrix with columns Node ID, parent, child 1, child 2, deme, node age
#' @param n.deme Number of distinct demes in the population
#' @param sd Standard deviation of the perturbations
#'
#' @return Updated extended data object with the proposal from the migration shift proposal
#'
#' @export

ed.mig.shift.2 <- function(ED, n.deme, sd = 1){
  migration.nodes <- ED[ is.na(ED[,4]) & (!is.na(ED[,3])) ,1]
  selected.node <- sample.vector(migration.nodes, 1)
  selected.row <- which(ED[,1] == selected.node)
  selected.time <- ED[selected.row, 6]

  child.node <- ED[selected.row, 3]
  child.row <- which(ED[,1] == child.node)
  child.time <- ED[child.row, 6]

  parent.node <- ED[selected.row, 2]
  parent.row <- which(ED[,1] == parent.node)
  parent.time <- ED[parent.row, 6]

  new.time <- rnorm.reflect(1, selected.time, sd, parent.time, child.time)
  ED[selected.row, 6] <- new.time

  return(list(ED = ED, prop.ratio = 1))
}
