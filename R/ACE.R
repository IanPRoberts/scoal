h <- function(){

phylo <- harris2010phylo
n.leaf <- length(phylo$tip.label)
x <- phylo$node.deme[1:n.leaf]
n.deme <- 5
color.palette <- rainbow(5) #c("red", "green", "blue", "orange", "purple")

#### Ancestral Character Estimation
# Estimates probability of deme of each coal node using maximum likelihood

ACE <- ace(x, phylo, type = "discrete")
plot(phylo, type = "p", TRUE, label.offset = 1)
tiplabels(pch = 22, bg = color.palette[as.numeric(x)], frame = "circle", cex = 2, adj = 1)
nodelabels(pie = ACE$lik.anc, piecol = color.palette, cex = 0.75)


#Realisation of the ACE phylo with coloured edges - sample edge colour proportional to colour of coalescent node
ace.phylo <- phylo

for (i in 1 : dim(ACE$lik.anc)[1]){
  ace.phylo$node.deme[58 + i] <- sample(5, 1, prob = ACE$lik.anc[i,])
}
structured.plot(ace.phylo)
tiplabels(pch = 22, bg = color.palette[as.numeric(x)], frame = "circle", cex = 2, adj = 1)
nodelabels(pie = ACE$lik.anc, piecol = color.palette, cex = 0.75)

#### Most Parsimonious Reconstruction
# Estimates most parsimonious reconstruction for deme of each coal node
MPR <- MPR(x, unroot(phylo), 5)
MPR.deme <- phylo$node.deme
MPR.deme[(n.leaf+1):(length(MPR.deme))] <- MPR[,1]
MPR.deme
}
