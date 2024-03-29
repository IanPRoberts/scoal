% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Likelihoods.R
\name{structured.likelihood}
\alias{structured.likelihood}
\title{Calculates the likelihood of a structured coalescent tree}
\usage{
structured.likelihood(phylo, effective.pop, gen.length, migration.matrix)
}
\arguments{
\item{phylo}{object of class \code{phylo} augmented with deme labels for each node}

\item{effective.pop}{effective population sizes from which the sample is taken}

\item{gen.length}{mean generation length of the sampled individuals}

\item{migration.matrix}{matrix of migration rates between demes}
}
\value{
A vector (\code{log-likelihood}, \code{likelihood}) giving the log-likelihood and likelihood of \code{phylo}
}
\description{
Returns the likelihood and log-likelihood of a tree generated under the
structured coalescent process
}
