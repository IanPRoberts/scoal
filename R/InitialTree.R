#' Generates an initial tree for MCMC
#'
#' Adds migration nodes to a structured coalescent topology to create a full
#' structured coalescent tree. Distribution of demes at each coalescent node
#'
#' Generates an initial migration history for a structured coalescent tree given demes of every tip using ancestral character evolution (ACE)
#'
#' @param phylo object of class \code{phylo} giving coalescent node topology
#' @param leaf.deme vector of demes for each leaf in \code{phylo}
#'
#' @return An object of class \code{phylo} augmented with the deme of each node in the tree.
#'
#' @export

initial.tree <- function(phylo, leaf.deme){
  ACE <- ace(leaf.deme, phylo, type = "discrete")
  n <- length(leaf.deme)

  phylo$node.deme <- rep(NA, 2*n - 1)
  ED <- phylo.to.ed(phylo)
  ED[1:n, 5] <- leaf.deme

  observed.demes <- as.numeric(colnames(ACE$lik.anc))
  n.observed.deme <- dim(ACE$lik.anc)[2]

  for (i in 1 : (n-1)){
    ED[n + i, 5] <- observed.demes[sample.int(n.observed.deme, 1, prob = ACE$lik.anc[i,])]
  }

  for (i in (1 : (2 * n - 1))[-(n+1)]){
    parent.row <- which(ED[,1] == ED[i,2])

    if (ED[i,5] != ED[parent.row, 5]){  #If demes don't agree, add migration unif along edge
      new.node <- max(ED[,1]) + 1
      new.loc <- runif(1, ED[parent.row, 6], ED[i,6])
      ED <- rbind(ED, c(new.node, ED[parent.row,1], ED[i,1], NA, ED[parent.row, 5], new.loc))
      which.child <- which(ED[parent.row, 3:4] == ED[i,1])
      ED[parent.row, 2 + which.child] <- new.node
      ED[i, 2] <- new.node
    }
  }

  phylo <- ed.to.phylo(ED)
  return(phylo)
}
