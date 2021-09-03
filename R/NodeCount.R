node.count <- function(phylo){
  demes <- unique(phylo$node.deme)
  n.deme <- length(demes)
  n <- length(phylo$tip.label)

  leaf.nodes <- 1:n
  migration.nodes <- phylo$edge[!phylo$edge[,1] %in% phylo$edge[duplicated(phylo$edge[,1]),1],1]
  n.migrations <- length(migration.nodes)
  coalescence.nodes <- ((n+1):(n.migrations+2*n-1))[!((n+1):(n.migrations+2*n-1) %in% migration.nodes)]

  M <- matrix(0, n.deme, n.deme)
  C <- numeric(n.deme)

  for (i in 1:n.deme){
    end.in.i <- migration.nodes[phylo$node.deme[migration.nodes] == i]
    origin.deme <- numeric(length(end.in.i))
    count <- 1
    for (j in end.in.i){
      origin.deme[count] <- phylo$node.deme[phylo$edge[which(phylo$edge[,1] == j),2]]
      count <- count + 1
    }
    M[,i] <- summary(factor(origin.deme, 1:3))
  }
  return(M)
}
