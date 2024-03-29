% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Local_likelihoods.R
\name{local_DTA_likelihood}
\alias{local_DTA_likelihood}
\title{Local DTA Likelihood}
\usage{
local_DTA_likelihood(st_labels, coal_rate, bit_mig_mat)
}
\arguments{
\item{st_labels}{subset of rows from an ED object giving all nodes in a subtree}

\item{coal_rate}{Vector of coalescent rates}

\item{bit_mig_mat}{Backward-in-time migration matrix for the phylogeny}
}
\value{
Likelihood & log_likelihood between events
}
\description{
Computes the DTA likelihood for a subtree
}
