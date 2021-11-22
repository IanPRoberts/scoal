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
