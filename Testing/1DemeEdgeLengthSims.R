n <- 10
N <- 1e5
n.deme <- 1

data <- matrix(0,nrow = n, ncol = 3)
data[,1] <- 1:n
data[,2] <- rep(2022, n) #runif(n, min = 2010, max = 2022)
data[,3] <- rep(1, n)

gen.length <- 1
effective.pop <- rep(1, n.deme)  #1 / rgamma(n.deme, shape = eff.pop.alpha, scale = 1 /eff.pop.beta)
migration.matrix <-matrix(1/10, n.deme, n.deme) #matrix(rgamma(n.deme^2, shape = mig.mat.alpha, scale = 1/mig.mat.beta), n.deme, n.deme)
diag(migration.matrix) <- 0

edge.lengths <- matrix(NA, N, n-1)

for (i in 1 : N){
  phylo <- Structured.sim(data, effective.pop, gen.length, n.deme, migration.matrix, FALSE)
  event.times <- sort(unique(node.depth.edgelength(phylo)))[1:n]
  edge.lengths[i,] <- diff(event.times)
}

dim <- ceiling(sqrt(n-1))

png("Testing/1DemeEdgeLengthSims.png", width = 2000, height = 1500)
layout(matrix(1:dim^2, dim, dim, byrow = TRUE))
for (i in 1 : (n-1)){
  upper <- ceiling(10*max(edge.lengths[,i]))/10
  hist(edge.lengths[,i], freq = FALSE, breaks = seq(0, upper, length.out = 100), main = paste("Time interval", i), xlab = "Edge lengths")
  lines(seq(0, upper, length.out = 100), dexp(seq(0, upper, length.out = 100), choose(i+1,2)), lwd = 2, lty = 2, col = "red")
}
dev.off()
