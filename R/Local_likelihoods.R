#' @rdname SC_likelihood
#' @param st_labels Reduced extended data structure consisting of the rows of ED corresponding to nodes in the selected subtree only
#' @param ED_DD Deme decomposition obtained from DemeDecomp()
#' @param ED_NC Node count obtained from NodeCount()
#'
#' @export

SC_local_like <- function(ED, coal_rate, bit_mig_mat, st_labels,
                          ED_NI = NodeIndices(ED), ED_DD = DemeDecomp(ED, nrow(bit_mig_mat), ED_NI), ED_NC = NodeCount(ED, nrow(bit_mig_mat), ED_NI)){
  n_deme <- length(coal_rate)

  event_ids <- st_labels[,1]
  limit_rows <- ED_NI[event_ids]
  limit_times <- ED[limit_rows, 6]

  event_times <- ED_DD$event.times[-1]
  subset_rows <- (event_times > min(limit_times)) & (event_times <= max(limit_times))
  k_subset <- subset(ED_DD$k, subset_rows)
  ti_subset <- subset(ED_DD$time.increments, subset_rows)

  deme_lengths <- colSums(k_subset * ti_subset)
  coal_constants <- colSums(k_subset * (k_subset-1) * ti_subset) * coal_rate / 2
  mm_row_sum <- rowSums(bit_mig_mat)
  log_mig_mat <- log(bit_mig_mat)
  diag(log_mig_mat) <- 0

  like <- sum(ED_NC$m * log_mig_mat) + sum(ED_NC$c * log(coal_rate)) - sum(coal_constants) - sum(deme_lengths * mm_row_sum)
  return(list(log.likelihood = like, likelihood = exp(like)))
}

#' @rdname DTA_likelihood
#' @param st_labels Reduced extended data structure consisting of the rows of ED corresponding to nodes in the selected subtree only
#'
#' @export

DTA_local_likelihood <- function(st_labels, coal_rate, bit_mig_mat){
  # Local DTA log likelhood function evaluated over a subtree like subtree$st_labels

  is_root <- !(st_labels[,2] %in% st_labels[,1])
  is_leaf <- is.na(st_labels[,3]) | !(st_labels[,3] %in% c(st_labels[,1], NA))

  if (ncol(st_labels) == 9){
    st_labels[is_root, 7] <- NA
    st_labels[is_leaf, 8:9] <- NA
  }

  st_labels[is_root, 2] <- NA
  st_labels[is_leaf, 3:4] <- NA
  return(ScaledDTALikelihoodC(st_labels, coal_rate, 1, bit_mig_mat, NodeIndices(st_labels)))
}

#' @rdname MTT_likelihood
#' @param st_labels Reduced extended data structure consisting of the rows of ED corresponding to coalescent nodes in the selected subtree only
#'
#' @export

MTT_local_likelihood <- function(ED, st_labels, bit_rates, ED_NI = NodeIndices(ED),
                                 eigen_vals = eigen(bit_rates)$values, eigen_vecs = eigen(bit_rates)$vectors, inverse_vecs = solve(eigen_vecs)){
  #Modify st_labels into valid ED structure (no reference below leaves or above root)
  st_labels <- ED[ED_NI[st_labels[,1]],]
  st_labels[1,c(2,7)] <- NA #Parents of st_root (in row 1) = NA

  is_leaf <- is.na(st_labels[,3]) | #Leaf of ED
    !((st_labels[,8] %in% st_labels[,1]) | (st_labels[,9] %in% st_labels[,1])) #Neither child coalescent node is in st_labels
  st_labels[is_leaf, c(3:4, 8:9)] <- NA #Children of st_leaves = NA

  #Extract migration events in subtree and add to st_labels
  is_st_mig <- (ED[,7] %in% st_labels[,1]) & #Parent coalescent node in subtree
    ((ED[,8] %in% st_labels[,1]) | (ED[,9] %in% st_labels[,1])) & #Child coalescent node in subtree
    is.na(ED[,4]) #Is migration

  st_labels <- rbind(st_labels, ED[is_st_mig,])

  return(MTT_likelihood(st_labels, bit_rates, NodeIndices(st_labels), eigen_vals, eigen_vecs, inverse_vecs))
}
