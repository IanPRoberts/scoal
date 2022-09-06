f <- function(){

  phylo <- ape::read.nexus("data/didelot2015.nex")
  node.color <- read.table("data/didelot2015.nex", skip = 4, nrows = 260, comment.char = "")
  tip.label <- gsub('.{17}$', '', node.color[,1])
  tip.order <- sapply(1:260, function(x) which(tip.label == phylo$tip.label[x]))

  node.color <- as.factor(gsub('^.{15}|.$', '', node.color[,1]))
  leaf.deme <- sapply(1:260, function(x) which(unique(node.color) == node.color[x]))

  phylo$node.deme <- c(leaf.deme[tip.order], rep(0, 259))
  ED <- phylo.to.ed(didelot2015phylo)
  saveRDS(phylo, file = "data/didelot2015phylo.RDS")
  saveRDS(ED, "data/didelot2015ED.RDS")
}

g <- function(){
  phylo <- ape::read.tree("data/harris2010.nwk")
  leaf.data <- read.csv("data/harris2010.csv")

  demes <- unique(leaf.data$Location)
  n.deme <- 5
  leaf.deme <- sapply(1:58, function(x) which(demes == leaf.data))

  tip.order <- as.numeric(phylo$tip.label)

  phylo$node.deme <- c(leaf.deme[tip.order], rep(0, 57))
  ED <- phylo.to.ed(phylo)
  saveRDS(phylo, file = "data/harris2010phylo.RDS")
  saveRDS(ED, "data/harris2010ED.RDS")


  ### Convert to .nex file
  # library(tibble)
  # library(treeio)
  # tree=read.tree('data/harris2010.nwk')
  # table=read.table('data/harris2010.csv',sep=',',header = T)
  # ti=tibble(label=as.character(table[,1]),location=table[,3])
  # td=full_join(tree,ti,by='label')
  # write.beast(td,'data/harris2010.nex')
}
