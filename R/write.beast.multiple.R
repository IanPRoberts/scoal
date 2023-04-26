#' write.beast.multiple
#'
#' Export multiple strphylo objects to BEAST NEXUS file
#'
#' @param ... either (i) a single object of class "strphylo"; (ii) multiple objects of class "strphylo" separated by commas; or (iii) a list of "strphylo" objects
#' @param file a file name specified either as a string; if left as default, outputs to standard output
#' @param translate logical; if TRUE, tip labels are replaced with tokens
#' @param thin integer from which the sample has been thinned
#' @param initial integer first index of sample
#'
#' @return None
#'
#' @export

write.beast.multiple <- function(..., file = stdout(), translate = TRUE, thin = 1, initial = 0){
  input <- list(...)
  n_tree <- length(input)

  header <- capture.output(write.nexus(input[[1]], file = stdout(), translate = TRUE))
  header[2] <- paste("[R-package scoal, ", date(), "]\n\n", sep = "")
  cat(header[-(length(header) - 0:1)], file = file, sep = "\n")

  for (x in 1 : n_tree){
    if (!inherits(input[[x]], "strphylo"))
    phylo <- input[[x]]
    class(phylo) <- "phylo"
    treedata <- treeio::as.treedata(phylo)
    treedata@data <- tidytree::tibble(type = paste0("\"", phylo$node.deme, "\""), node = 1:length(phylo$node.deme))
    newick_tree <- treeio::write.beast.newick(treedata)

    cat("tree STATE_", initial + thin * (x-1), " = ", newick_tree, "\n", file = file, append = TRUE, sep = "")
  }
  cat(header[length(header)])
}
