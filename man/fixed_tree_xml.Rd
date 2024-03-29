% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BEAST2_xml.R
\name{fixed_tree_xml}
\alias{fixed_tree_xml}
\title{Generates BEAST2 xml for MultiType Tree}
\usage{
fixed_tree_xml(
  strphylo,
  n_deme,
  coal_rate,
  bit_mig_mat,
  N = 1e+07,
  thin = 1000,
  tree_thin = thin,
  con = stdout(),
  run_name = "$(filebase)",
  BEAST2_package = "MTT",
  priors = "default",
  cr_shape = NULL,
  cr_rate = NULL,
  mm_shape = NULL,
  mm_rate = NULL
)
}
\arguments{
\item{strphylo}{structured phylo object giving initialisation condition for fixed tree run with MultiType Tree}

\item{coal_rate}{initial estimate for coalescent rates}

\item{N}{total number of MCMC iterations (including burn-in)}

\item{thin}{thinning increment for the MCMC}

\item{tree_thin}{Thinning rate for tree samples to be saved (Default value equal to thin)}

\item{con}{A connection object or a character string giving the location for the output xml file (stdout() prints to console)}

\item{run_name}{Name for logger files to be saved as}

\item{BEAST2_package}{Select package to prepare xml file for, either MTT (MultiTypeTree) or BASTA}

\item{priors}{either 'default' in which lognormal priors are used or 'gamma' in which gamma/inverse gamma priors are used (and prior parameters must be specified)}

\item{bit_mig_rate}{initial estimate for backward-in-time migration matrix}
}
\value{
output file or file content on screen
}
\description{
Generates a .xml file to run BEAST2 using either MultiTypeTree or BASTA with a fixed tree.
Initialises with all coalescent rates the same and all backward-in-time migration rates the same
}
