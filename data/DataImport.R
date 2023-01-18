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
  leaf.deme <- sapply(1:58, function(x) which(demes == leaf.data$Location[x]))

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

h <- function(){
  phylo <- ape::read.tree("data/dudas2017.nwk")
  leaf.data <- read.csv("data/dudas2017.csv")

  demes <- unique(leaf.data$country)
  leaf.deme <- sapply(1:1610, function(x) which(demes == leaf.data$country[x]))
  tip.order <- sapply(1:1610, function(x) which(phylo$tip.label[x] == leaf.data$id))

  phylo$node.deme <- c(leaf.deme[tip.order], rep(0, 1609))
  ED <- phylo.to.ed(phylo)
  saveRDS(phylo, file = "data/dudas2017phylo.RDS")
  saveRDS(ED, "data/dudas2017ED.RDS")
}

k <- function(){
  phylo <- ape::read.nexus("data/MERS.nex")
  tip.data <- read.table("data/MERS.nex", skip = 5, nrows = 274, comment.char = "")
  leaf.data <- stringr::str_split(tip.data[,1], "\\|")

  host_data <- sapply(1:274, function(x) leaf.data[[x]][3])
  tip_label <- sapply(1:274, function(x) leaf.data[[x]][2])

  # leaf.deme <- numeric(274)
  # leaf.deme[grep("human", tip.data[,1])] <- 1
  # leaf.deme[grep("camel", tip.data[,1])] <- 2
  leaf.deme <- as.numeric(as.factor(host_data))

  phylo$node.deme <- c(leaf.deme, rep(0, 273))
  ED <- phylo.to.ed(phylo)
  structured.plot(ED)
  saveRDS(phylo, file = "data/MERSphylo.RDS")
  saveRDS(ED, "data/MERS_ED.RDS")
}
