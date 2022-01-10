f <- function(){

didelot2015phylo <- ape::read.nexus("data/didelot2015.nex")
node.color <- read.table("Data/didelot2015.nex", skip = 4, nrows = 260, comment.char = "")
node.color <- sub('...', '', gsub('[^0-9.]','',node.color[,]))
n.leaf.deme <- length(unique(node.color))
leaf.deme <- rep(0, 260)

for (i in 1:260){
  leaf.deme[i] <- (1:n.leaf.deme)[which(unique(node.color) == node.color[i])]
}
leaf.deme <- as.numeric(leaf.deme)

didelot2015phylo$node.deme <- c(node.color, rep(0, 259))
didelotED <- phylo.to.ed(didelot2015phylo)
View(didelotED)


save(didelot2015phylo, file = "data/didelot2015phylo.rda")
}

g <- function(){
  harris2010phylo <- ape::read.tree("data/harris2010.nwk")
  harris2010NodeData <- read.csv("data/harris2010.csv")
  demes <- levels(as.factor(harris2010NodeData$Location)) #levels(harris2010NodeData$Location)
  n.deme <- length(demes)

  n.leaf <- dim(harris2010NodeData)[1]
  node.deme <- rep(n.deme+1, 2 * n.leaf - 1)
  for (i in 1:n.leaf){
    node.deme[i] <- (1:n.deme)[which(demes == harris2010NodeData$Location[i])]
  }

  node.deme[1:n.leaf] <- node.deme[as.numeric(harris2010phylo$tip.label)]  #Reorder to satisfy weird ordering of tip.labels

  harris2010phylo$node.deme <- node.deme
  harrisED <- phylo.to.ed(harris2010phylo)



  ### Convert to .nex file
  library(tibble)
  library(treeio)
  tree=read.tree('data/harris2010.nwk')
  table=read.table('data/harris2010.csv',sep=',',header = T)
  ti=tibble(label=as.character(table[,1]),location=table[,3])
  td=full_join(tree,ti,by='label')
  write.beast(td,'data/harris2010.nex')
}
