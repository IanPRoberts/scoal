#' Coalescent Node Pie Charts
#'
#' Tabulates the frequency of each deme at each leaf and coalescent node of a phylogenetic tree
#'
#'
#' @param ED_sample Sample of migration histories on same topology (ED with 9 columns)
#' @param plot Logical value whether to plot tree with superimposed pie charts
#' @param plot_ED (optional) ED to superimpose with pie charts (clear migration history if left blank)
#'
#' @return Returns deme frequencies at each location
#'
#' @export

coalescent_node_pie_charts <- function(ED_list, plot = TRUE, plot_ED = matrix(NA, 0, 9), cex = 0.5, ...){
  n_trees <- length(ED_list)
  n_deme <- max(ED_list[[1]][,5])

  #Strip migration events from tree
  topology <- ED_list[[1]]
  topology <- topology[(is.na(topology[,3])) | (!is.na(topology[,4])),]
  topology[,2:4] <- topology[,7:9]
  topology[!is.na(topology[,3]), 5] <- 0

  #Construct deme frequency at each coalescent node in ascending node age
  deme_freq <- matrix(0, nrow(topology), n_deme)

  for (tree_id in n_trees : 1){ #Loop in reverse tree order to leave node_order as order(topology[,6]) after final iteration
    ED <- ED_list[[tree_id]]
    ED <- ED[(is.na(ED[,3])) | (!is.na(ED[,4])),]
    ED[,2:4] <- ED[,7:9]

    node_order <- order(ED[,6]) #Store entries of deme_freq in ascending age (root = 0, newest leaf = max(ED[,6]))

    for (row_id in 1 : nrow(ED)){ #Loop robust only for non-simultaneous events, i.e. coalescent nodes
      deme_freq[row_id, ED[node_order[row_id], 5]] <- deme_freq[row_id, ED[node_order[row_id], 5]] + 1
    }
  }

  #Remove leaf deme frequencies from deme_freq (should be a.s. one deme throughout run)
  deme_freq <- deme_freq[!is.na(topology[node_order, 4]),]

  if (nrow(plot_ED) == 0){
    plot_ED <- topology
    plot_ED[,5] <- 0
  }

  plot_ED_coal_nodes <- !is.na(plot_ED[,4])
  plot_ED_coal_node_order <- order(plot_ED[plot_ED_coal_nodes, 6])


  rownames(deme_freq) <- plot_ED[plot_ED_coal_nodes, 1][plot_ED_coal_node_order]

  if (plot){
    structured.plot(plot_ED)
    pie_plot <- rowSums(deme_freq == 0) < n_deme - 1 #Logical on whether 100% same deme observed (in which case no pie chart plotted!)
    nodelabels(node = as.numeric(rownames(deme_freq))[pie_plot],
               pie = deme_freq[pie_plot,]/rowSums(deme_freq[pie_plot,]),
               cex = cex)
  }

  return(list(ED = topology, node_freq = deme_freq))
}
