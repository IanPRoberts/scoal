#' ED Distance
#'
#' Computes the dissimilarity between a pair of migration histories on the same
#' tree via total branch length in different demes
#'
#' @param ED1 Extended data representation of first migration history
#' @param ED2 Extended data representation of second migration history
#'
#' @return Computed distance
#'
#' @export

ED.dist <- function(ED1, ED2, node_indices_1 = NA, node_indices_2 = NA){
  n_deme <- max(c(ED1[,5], ED2[,5]))
  event_times_1 <- ED1[,6]
  event_times_2 <- ED2[,6]
  pooled_event_times <- sort(unique(c(event_times_1, event_times_2)))
  pooled_time_incs <- diff(pooled_event_times)


  if (any(is.na(node_indices_1))){
    node_indices_1 <- NodeIndicesC(ED1)
  }

  if (any(is.na(node_indices_2))){
    node_indices_2 <- NodeIndicesC(ED2)
  }

  deme_decomp_1 <- DemeDecompC(ED1, n_deme, node_indices_1)
  deme_decomp_2 <- DemeDecompC(ED2, n_deme, node_indices_2)

  k_diff <- matrix(0, length(pooled_event_times)-1, n_deme)

  count_1 <- count_2 <- 1

  for (i in 2 : length(pooled_event_times) - 1){
    k_diff[i,] <- abs(deme_decomp_1$k[count_1,] - deme_decomp_2$k[count_2,])
    if (pooled_event_times[i+1] == deme_decomp_1$event.times[count_1 + 1]){
      count_1 <- count_1 + 1
    }

    if (pooled_event_times[i+1] == deme_decomp_2$event.times[count_2 + 1]){
      count_2 <- count_2 + 1
    }
  }

  tree_length <- sum(deme_decomp_1$k * deme_decomp_1$time.increments)


  score <- sum(k_diff * pooled_time_incs) / (2 * tree_length)  #k_diff double counts - remove from one deme and add to another

  return(distance = score)
}

#' Medoid ED
#'
#' Computes the medoid migration history from a sample (minimal mean distance
#' from all other sampled migration histories)
#'
#' @param ED_sample Sample of ED objects containing migration histories on the same topology
#' @param plot_medoid Logical whether to plot the medoid migration history
#'
#' @return Medoid ED
#'
#' @export

medoid.ED <- function(ED_sample, plot_medoid = FALSE){
  n_samples <- length(ED_sample)

  dists <- matrix(0, n_samples, n_samples)

  for (i in 1 : n_samples){
    ED1 <- ED_sample[[i]]
    node_indices_1 <- rep(0, max(ED1[,1]))
    for (j in 1 : dim(ED1)[1]){
      node_indices_1[ED1[j,1]] <- j
    }
    for (j in ((i+1) : n_samples)){
      ED2 <- ED_sample[[j]]
      dists[i,j] <- dists[j,i] <- ED.dist(ED1, ED2, node_indices_1)
    }
  }

  mean_dists <- rowMeans(dists)
  medoid_ED <- ED_sample[[min(which(mean_dists == min(mean_dists)))]]

  if (plot_medoid){
    structured.plot(medoid_ED)
  }

  return(medoid = medoid_ED)
}


#' ED Distance (coalescent node)
#'
#' Computes the dissimilarity between a pair of migration histories on the same
#' tree via number of coalescent nodes in different demes between the two histories
#'
#' @param ED1 Extended data representation of first migration history
#' @param ED2 Extended data representation of second migration history
#'
#' @return Computed distance
#'
#' @export

ED.dist.coal.node <- function(ED1, ED2, node_indices_1 = NA, node_indices_2 = NA){
  if (any(is.na(node_indices_1))){
    node_indices_1 <- NodeIndicesC(ED1)
  }

  if (any(is.na(node_indices_2))){
    node_indices_2 <- NodeIndicesC(ED2)
  }

  n_deme <- max(c(ED1[,5], ED2[,5]))

  coal_labels <- ED1[((!is.na(ED1[,3])) & (!is.na(ED1[,4]))),1]
  score <- sum(ED1[node_indices_1[coal_labels], 5] != ED2[node_indices_2[coal_labels], 5])

  return(score)
}
