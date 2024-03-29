% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utils.R
\name{child_nodes}
\alias{child_nodes}
\title{Child nodes}
\usage{
child_nodes(node, edge, node.ages)
}
\arguments{
\item{node}{Node to find parent node of}

\item{edge}{Edge matrix from a \code{phylo} object}

\item{node.ages}{Time of each node in the \code{phylo} object}
}
\value{
Vector of node IDs of any child nodes. If the input node is a leaf, returns NA
}
\description{
Given a node from a \code{phylo} object, identifies the node IDs of any direct
child nodes
}
