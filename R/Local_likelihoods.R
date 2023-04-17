#' Local Structured Coalescent Likelihood
#'
#' Computes the structured coalescent likelihood locally between two specified events
#'
#'
#' @param ED Extended data representation of a phylogeny including migration history
#' @param coal_rate Vector of coalescent rates
#' @param bit_mig_mat Backward-in-time migration matrix for the phylogeny
#' @param event_ids Pair of event IDs to compute likelihood between
#' @param node_indices
#' @param deme_decomp
#' @return Likelihood & log_likelihood between events
#'
#' @export

SC_local_like <- function(ED, coal_rate, bit_mig_mat, event_ids, node_indices = NULL, deme_decomp = NULL, node_count = NULL){
    n_deme <- length(coal_rate)
    if (is.null(node_indices)) node_indices <- NodeIndicesC(ED)
    if (is.null(deme_decomp)) deme_decomp <- DemeDecompC(ED, n_deme, node_indices)
    if (is.null(node_count)) node_count <- NodeCountC(ED, n_deme, node_indices)

    limit_rows <- node_indices[event_ids]
    limit_times <- ED[limit_rows, 6]

    event_times <- deme_decomp$event.times[-1]
    subset_rows <- (event_times > min(limit_times)) & (event_times <= max(limit_times))
    k_subset <- subset(deme_decomp$k, subset_rows)
    ti_subset <- subset(deme_decomp$time.increments, subset_rows)

    deme_lengths <- colSums(k_subset * ti_subset)
    coal_constants <- colSums(k_subset * (k_subset-1) * ti_subset) * coal_rate / 2
    mm_row_sum <- rowSums(bit_mig_mat)
    log_mig_mat <- log(bit_mig_mat)
    diag(log_mig_mat) <- 0

    like <- sum(node_count$m * log_mig_mat) + sum(node_count$c * log(coal_rate)) - sum(coal_constants) - sum(deme_lengths * mm_row_sum)
    return(list(log.likelihood = like, likelihood = exp(like)))
}


