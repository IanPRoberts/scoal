#' Subtree centre sampling
#'
#' Methods to identify the subtree centre for MCMC
#'
#' @details
#' st_centre_location and st_centre_location2 sample a lineage at time centre_age uniformly.
#' st_centre_unif selects a subtree centre location uniformly on the tree
#'
#'
#' @param ED ED representation of a phylogeny including migration history
#' @param ED_NI NodeIndicesC(ED)
#' @param centre_age Age of the subtree centre (e.g. sampled uniformly between oldest tip and root)
#' @param root_row (optional) row of ED corresponding to the root of the tree
#'
#' @return List consisting of st_child, the label of the child event in ED and st_centre_loc giving the location in [0,1] along the branch segment above st_child where the subtree centre is located
#'
#' @export

# Direct computation of candidate nodes by looking at node ages ('Breadth-first search')

st_centre_location <- function(ED, ED_NI, centre_age, root_row = which(is.na(ED[,2]))){
  parent_ages <- ED[ED_NI[ED[,2]], 6] #Ages of parent nodes
  parent_ages[root_row] <- 0
  st_centre_children <- ED[(ED[,6] > centre_age) & (parent_ages < centre_age), 1]

  st_child <- st_centre_children[sample(length(st_centre_children), 1)] #Robust for edge case with long leaf branch leaving 1 lineage
  st_child_row <- ED_NI[st_child]
  st_child_age <- ED[st_child_row, 6]
  st_child_parent <- ED[st_child_row, 2]
  st_centre_loc <- (st_child_age - centre_age) / (st_child_age - ED[ED_NI[st_child_parent], 6])#Proportion along branch segment above st_child for subtree centre

  return(list(st_child=ED[st_child_row, 1], st_centre_loc=st_centre_loc, candidate_children=st_centre_children))
}

#' @rdname st_centre_location

#### Root-leaf traversal to identify candidate nodes in the tree ('Depth first search')
st_centre_location2 <- function(ED, ED_NI, centre_age, root_row = which(is.na(ED[,2]))){
  active_rows <- root_row
  st_centre_child_rows <- c() #Vector to store candidate subtree children

  while (length(active_rows) > 0){
    new_active_rows <- c()

    for (active_row in active_rows){
      for (child_id in 1:2){
        if (!is.na(ED[active_row, 2 + child_id])){ #If child is not NA
          child_row <- ED_NI[ED[active_row, 2 + child_id]]

          if (ED[child_row, 6] < centre_age){ #If child node is above subtree centre_age
            new_active_rows <- append(new_active_rows, child_row)
          } else {
            st_centre_child_rows <- append(st_centre_child_rows, child_row)
          }
        }
      }
    }

    active_rows <- new_active_rows
  }


  st_child_row <- sample(st_centre_child_rows, 1) #Sample child node uniformly from available options
  st_child_age <- ED[st_child_row, 6]
  st_child_parent <- ED[st_child_row, 2]
  st_centre_loc <- (st_child_age - centre_age) / (st_child_age - ED[ED_NI[st_child_parent], 6])#Proportion along branch segment above st_child for subtree centre

  return(list(st_child=ED[st_child_row, 1], st_centre_loc=st_centre_loc, candidate_children=ED[st_centre_child_rows,1]))
}


#' @rdname st_centre_location

st_centre_location_unif <- function(ED, ED_NI){
  parent_ages <- ED[ED_NI[ED[,2]], 6] #Ages of parent nodes
  parent_ages[root_row] <- 0

  st_child <- sample(ED[,1], 1, prob=ED[,6] - parent_ages)
  st_child_row <- ED_NI[st_child]
  st_child_age <- ED[st_child_row, 6]
  st_child_parent <- ED[st_child_row, 2]
  st_centre_loc <- (st_child_age - centre_age) / (st_child_age - ED[ED_NI[st_child_parent], 6])#Proportion along branch segment above st_child for subtree centre

  return(list(st_child=ED[st_child_row, 1], st_centre_loc=st_centre_loc))
}
