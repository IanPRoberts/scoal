% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DTA_sampling.R
\name{DTA_rejection}
\alias{DTA_rejection}
\title{DTA Sampling}
\usage{
DTA_rejection(ED, fit_mig_mat, time_scale, node_indices, max_attempts = 10000)
}
\arguments{
\item{ED}{Extended ED representation of a phylogeny including migration history}

\item{fit_mig_mat}{Forward-in-time migration matrix for the phylogeny}

\item{time_scale}{Time scale for the migration matrix}
}
\value{
List of proposed DTA migration histories with the same leaves as input ED
}
\description{
Samples a migration history under the DTA model with deme fixed at all leaves using rejection sampling
}
