#' Subtree sampling
#'
#' Methods to sample subtrees and associated sub-migration histories from
#' structured genealogies
#'
#' @details
#' \code{radius_subtree} gives a radius-based subtree where points are included in the subtree up to a maximum distance from a fixed centrepoint.
#'
#' @param ED Extended data representation of a structured genealogy
#' @param ED_NI Node indices obtained from NodeIndices()
#' @param st_radius Maximum radius of the sampled subtree. If a leaf or the root are reached, radius may be smaller than st_radius
#' @param st_child (optional) Child node of the branch segment containing the subtree centre. If omitted, child node is selected uniformly at random
#' @param st_centre_loc (optional) Location of the subtree centre between st_child and st_child's parent. If omitted, sampled uniformly from (0,1)
#' @param edge_lengths (optional) Vector of lengths of branch segments between every node and its parent in the order of rows of ED
#'
#' @return List consisting of the structured genealogy with added self-migrations at the subtree root and leaves (if necessary) and ED object consisting only of subtree nodes (including any added self-migrations)
#'
#' @export


radius_subtree <- function(ED, ED_NI=NodeIndices(ED), st_radius, st_child = NULL, st_centre_loc = NULL, edge_lengths = NULL){
  if (is.null(st_child)){
    st_child <- sample(ED[,1], 1, prob = edge_lengths)
  }

  if (is.null(st_centre_loc)){
    st_centre_loc <- runif(1, 0, 1)
  }

  if (is.null(edge_lengths)){
    edge_lengths <- ED[,6] - ED[ED_NI[ED[,2]], 6]
    edge_lengths[is.na(edge_lengths)] <- 0
  }

  st_child_row <- ED_NI[st_child]
  st_root <- st_child
  st_root_row <- st_child_row

  st_centre_age <- ED[st_child_row, 6] - edge_lengths[st_child_row] * st_centre_loc

  root_node <- ED[is.na(ED[,2]), 1]
  root_row <- ED_NI[root_node]

  st_root_age <- max(ED[root_row, 6], st_centre_age - st_radius)

  #st_labels contains node_ID, parent, children and distance from st_centre
  st_labels <- matrix(NA, nrow = 0, ncol = 2,
                      dimnames = list(NULL, c("Node_ID", "Node_dist")))

  while (ED[st_root_row, 6] > st_root_age){
    st_labels <- rbind(st_labels,
                       ED[st_root_row, c(1, 6)])
    st_root <- ED[st_root_row, 2]
    st_root_row <- ED_NI[st_root]
  }

  if (ED[st_root_row, 6] < st_root_age){
    #Add virtual migration as st_root if nEDed
    max_label <- max(ED[,1]) + 1
    which_child <- which(ED[st_root_row, 3:4] == st_labels[nrow(st_labels), 1])
    st_root_child <- ED[st_root_row, 2 + which_child]

    if (is.na(ED[st_root_row, 4])){ #If st_root is a migration, parent coal is parent coal of st_root
      parent_coal <- ED[st_root_row, 7]
    } else { #Else parent coal is st_root
      parent_coal <- st_root
    }

    ED <- rbind(ED,
                 c(max_label, #Label
                   st_root, #Parent
                   ED[st_root_row, 2 + which_child], #Child 1
                   NA, #Child 2
                   ED[ED_NI[st_root_child], 5], #Deme
                   st_root_age, #Node age
                   parent_coal, #Parent coal
                   ED[st_root_row, 7 + which_child], #Child coal 1
                   NA #Child coal 2
                 ))

    ED[ED_NI[st_root_child], 2] <- max_label
    ED[st_root_row, 2 + which_child] <- max_label
    ED_NI[max_label] <- nrow(ED)

    st_labels <- rbind(st_labels,
                       ED[nrow(ED), c(1, 6)])
  } else { #No virtual migration required as st_root is global root
    max_label <- max(ED[,1])
    st_labels <- rbind(st_labels,
                       ED[st_root_row, c(1, 6)])
  }

  st_labels[, 2] <- st_centre_age - st_labels[, 2] #abs(st_centre_age - st_labels[, 2])

  current_parents <- st_labels[,1]

  # if (st_labels[1,2] > st_radius){ #If st_child age < st_leaf_age
  if (st_labels[1,2] < -st_radius){
    current_parents <- current_parents[-1]
  }

  while (length(current_parents) > 0){
    new_parents <- numeric(0)
    for (node in current_parents){
      node_row <- ED_NI[node]
      which_child <- which(!(na.omit(ED[node_row, 3:4]) %in% st_labels[,1]))
      node_dist <- st_labels[which(st_labels[,1] == node), 2]

      if (length(which_child) > 0){
        for (child_id in which_child){
          child_row <- ED_NI[ED[node_row, 2 + child_id]]

          #### NED to remove new_parents below st_leaf_age and global leaves
          # child_dist <- node_dist + edge_lengths[child_row]
          child_dist <- node_dist - edge_lengths[child_row]
          st_labels <- rbind(st_labels,
                             c(ED[child_row, 1], child_dist))

          if ((abs(child_dist) < st_radius) & (!is.na(ED[child_row, 3]))){
            new_parents <- append(new_parents, ED[child_row, 1])
          }
        }
      }
    }
    current_parents <- new_parents
  }

  for (row_id in 1 : nrow(st_labels)){
    # if (st_labels[row_id, 2] > st_radius){
    if (st_labels[row_id, 2] < - st_radius){
      #If node is greater than distance st_radius from st_centre add virtual migration
      max_label <- max_label + 1

      node_row <- ED_NI[st_labels[row_id, 1]]

      if ((is.na(ED[node_row, 4])) & (!is.na(ED[node_row, 3]))){
        #Current node is a migration -> child coal is child coal of current node
        child_coal <- ED[node_row, 8]
      } else {
        #Current node is a coalescent or leaf -> child coal is current node
        child_coal <- st_labels[row_id, 1]
      }

      ED <- rbind(ED,
                   c(max_label, #Label
                     ED[node_row, 2], #Parent
                     ED[node_row, 1], #Child 1
                     NA, #Child 2
                     ED[node_row, 5], #Deme
                     ED[node_row, 6] + st_radius - abs(st_labels[row_id, 2]), #Node age
                     ED[node_row, 7], #Parent coal
                     child_coal, #Child coal 1
                     NA #Child coal 2
                   ))

      parent_row <- ED_NI[ED[node_row, 2]]
      ED[node_row, 2] <- max_label
      which_child <- which(ED[parent_row, 3:4] == st_labels[row_id, 1])
      ED[parent_row, 2 + which_child] <- max_label

      ED_NI[max_label] <- nrow(ED)

      st_labels[row_id,] <- c(ED[nrow(ED), 1], - st_radius)
    }
  }

  return(list(ED = ED, st_labels = ED[ED_NI[st_labels[,1]],]))
}



#' @rdname radius_subtree
#'
#' @param st_depth Maximum number of descendant generations of coalescent nodes to include in the subtree (default 1)
#' @param selected_node (optional) Coalescent node to center subtree on. If omitted, a coalescent node is selected uniformly at random
#'
#' @details
#' \code{coal_node_subtree} gives a coalescent node-based subtree consisting of the parent branch and child branches up to \code{st_depth} generations of coalescent nodes below a central coalescent node
#'
#'


coal_node_subtree <- function(ED, ED_NI = NodeIndices(ED), st_depth = 1, selected_node = NULL){
  if (is.null(selected_node)){
    coal_nodes <- ED[!is.na(ED[,4]), 1]
    selected_node <- sample(coal_nodes, 1)
  }

  selected_row <- NI[selected_node]
  st_root <- ED[selected_row, 7]
  max_label <- max(ED[,1]) + 1

  if (is.na(st_root)){
    st_root <- selected_node
    st_labels <- st_root
  } else {
    st_labels <- c(st_root, selected_node)
  }

  active_rows <- selected_row

  for (gen in 1 : st_depth){
    child_coals <- unname(na.omit(ED[active_rows, 8:9]))
    st_labels <- append(st_labels,
                        child_coals)
    active_rows <- NI[child_coals]
  }


  return(list(ED = ED, st_labels = ED[NI[st_labels],]))
}
