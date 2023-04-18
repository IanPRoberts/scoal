#' Plot Structured Coalescent Tree
#'
#' Plots a structured coalescent tree input as a phylo object augmented with node demes
#'
#' @param phylo phylo object augmented with node demes
#'
#' @export

structured.plot <- function(x, n_deme = NA, time_axis = FALSE, root_time = NA, ...){
  if (! ('phylo' %in% class(x))){
    phylo <- ed.to.phylo(x)
  } else{
    phylo <- x
  }

  edge <- phylo$edge

  if (is.na(n_deme)){
    n_deme <- max(unique(phylo$node.deme))
  }

  if (any(phylo$node.deme == 0)){
    color.palette <- c(rainbow(n_deme), "black")
    phylo$node.deme[phylo$node.deme == 0] <- n_deme + 1
  } else{
    color.palette <- rainbow(n_deme)
  }

  edge.color <- color.palette[phylo$node.deme[edge[,2]]]

  plot(phylo, edge.color = edge.color, edge.width = 2, show.tip.label = FALSE, ...)

  if (time_axis){
    if (is.na(root_time)){
      root_time <- 0
    }
    axisPhylo(root.time = root_time, backward = FALSE)
  }
}


#' @export
plot.str_phylo <- function(x, n_deme = NULL, time_axis = FALSE, root_time = NULL, ...){
  if (is.null(n_deme)){
    n_deme <- max(x$node.deme)
  }

  color.palette <- c(rainbow(n_deme), "black")
  x$node.deme[x$node.deme == 0] <- n_deme + 1

  edge_color <- color.palette[x$node.deme[x$edge[,2]]]

  plot.phylo(x, edge.color = edge_color, edge.width = 2, ...)

  if (time_axis){
    if (is.null(root_time)){
      root_time <- 0
    }
    axisPhylo(root.time = root_time, backward = FALSE)
  }
}

#' @export
plot.ED <- function(x, n_deme = NULL, time_axis = FALSE, root_time = NULL, ...){
  plot(as.phylo(x), n_deme, time_axis, root_time, ...)
}

#' @export
plot.EED <- function(x, n_deme = NULL, time_axis = FALSE, root_time = NULL, ...){
  plot(as.phylo(x), n_deme, time_axis, root_time, ...)
}
