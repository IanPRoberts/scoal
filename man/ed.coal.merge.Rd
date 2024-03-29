% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ED_proposals.R, R/Ewing_ED_proposals.R
\name{ed.coal.merge}
\alias{ed.coal.merge}
\title{Coalescent Node Merge Proposal}
\usage{
ed.coal.merge(ED, n.deme, node.indices)

ed.coal.merge(ED, n.deme, node.indices)
}
\arguments{
\item{ED}{Extended data object; matrix with columns Node ID, parent, child 1, child 2, deme, node age}

\item{n.deme}{Number of distinct demes in the population}

\item{node.indices}{Vector giving row indices for node labels}
}
\value{
Updated extended data object with the proposal from the migration pair birth move

Updated extended data object with the proposal from the migration pair birth move
}
\description{
Performs a coalescent node merge move (Ewing et al. 2004). Merges two migration
nodes immediately below a coalescent node and place above coalescent node

Performs a coalescent node merge move (Ewing et al. 2004). Merges two migration
nodes immediately below a coalescent node and place above coalescent node
}
