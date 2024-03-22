#' Fitch Algorithm
#'
#' Computes the minimum number of migration events required to explain a structured phylogeny
#'
#' @param ED Structured phylogeny giving structured tree with labelled tip demes
#' @param node_indices Output from NodeIndicesC()
#'
#' @export

fitch <- function(ED, node_indices){
  if (ncol(ED) == 9) ED <- ED[,1:6]
  n_deme <- max(ED[,5]) #Maximum number of required demes = max observed deme

  top_ED <- strip.history(ED, node_indices)
  top_ED <- cbind(top_ED, matrix(0, nrow(top_ED), n_deme)) #Add columns corresponding to optimal states to top_ED. Code 0 for FALSE, 1 for TRUE
  top_NI <- NodeIndicesC(top_ED)
  leaf_nodes <- top_ED[top_ED[,5] > 0, 1]
  active_rows <- top_NI[leaf_nodes]

  ####### Root-ward sweep
  while(length(active_rows) > 0){
    for (row in active_rows){
      if (is.na(top_ED[row, 3])){ #Leaf optimal state == sampled deme
        top_ED[row, 6 + top_ED[row, 5]] <- TRUE
      } else { #Coalescent node optimal states either intersection (if not empty) or union of child optimal states
        child_rows <- top_NI[top_ED[row, 3:4]]
        intersection <- (top_ED[child_rows[1], 6 + 1 : n_deme]) & (top_ED[child_rows[2], 6 + 1 : n_deme])
        if (any(intersection)){ #Intersection not empty
          top_ED[row, 6 + 1 : n_deme] <- intersection
        } else {
          top_ED[row, 6 + 1 : n_deme] <- (top_ED[child_rows[1], 6 + 1 : n_deme]) | (top_ED[child_rows[2], 6 + 1 : n_deme])
        }
      }
    }
    active_rows <- na.omit(top_NI[top_ED[active_rows, 2]])
  }

  ######## Leaf-ward sweep
  root_row <- which(is.na(ED[,2])) #root row
  top_ED[root_row, 5] <- sample.int(n_deme, 1, prob = top_ED[root_row, 6 + 1:n_deme])
  min_migs <- 0
  active_rows <- top_ED[root_row, 3:4]

  while (length(active_rows) > 0){
    for (row in active_rows){
      parent_row <- top_NI[top_ED[row, 2]]
      parent_deme <- top_ED[parent_row, 5]

      if (top_ED[row, 6 + parent_deme] == 1){
        top_ED[row, 5] <- parent_deme
      } else {
        top_ED[row, 5] <- sample.int(n_deme, 1, prob = top_ED[row, 6 + 1:n_deme])
        min_migs <- min_migs + 1
      }
    }
    active_rows <- na.omit(top_NI[top_ED[active_rows, 3:4]])
  }

  return(list(min_migs = min_migs, ED = top_ED))
}
