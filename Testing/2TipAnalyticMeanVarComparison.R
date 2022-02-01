# Symmetric island model analytic mean time to coalescence from Section 13.2 of ST341 lnotes

n <- 2
N <- 1e5
n.deme <- 3
m <- 1/2 #Total migration rate

#Start in same deme
data.within <- matrix(c(1:2, rep(0,2), 1,1),nrow = n, ncol = 3)
data.between <- matrix(c(1:2, rep(0,2), 1,2),nrow = n, ncol = 3)

gen.length <- 1
effective.pop <- rep(1, n.deme)
migration.matrix <-matrix(m/(n.deme - 1), n.deme, n.deme)
diag(migration.matrix) <- 0

tree.height <- matrix(NA, N, 2)

for (i in 1 : N){
  phylo <- Structured.sim(data.within, effective.pop, gen.length, n.deme, migration.matrix, FALSE)
  ED <- phylo.to.ed(phylo)
  tree.height[i, 1] <- ED[1,6]
  phylo <- Structured.sim(data.between, effective.pop, gen.length, n.deme, migration.matrix, FALSE)
  ED <- phylo.to.ed(phylo)
  tree.height[i, 2] <- ED[1,6]
}

mean.tree.height.within <- n.deme
mean.tree.height.between <- n.deme + (n.deme - 1)/m

print(c(mean.tree.height.within, mean(tree.height[,1])), digits = 10)
print(c(mean.tree.height.between, mean(tree.height[,2])), digits = 10)
