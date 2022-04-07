set.seed(1)
N <- 1e6
N0 <- 1e3
lambda <- 20
n <- 10
n.deme <- 5


#############################################
data <- matrix(0,nrow = n, ncol = 3)
data[,1] <- 1:n
data[,2] <- rep(2022, n) #runif(n, min = 2010, max = 2022)
data[,3] <- sample.int(n.deme, n, replace = TRUE)

#Simulation Parameters
gen.length <- 1
effective.pop <- rep(1, n.deme)
migration.matrix <-matrix(1/10, n.deme, n.deme)
diag(migration.matrix) <- 0

phylo <- Structured.sim(data, effective.pop, gen.length, n.deme, migration.matrix, FALSE)
ED <- phylo.to.ed(phylo)
node.indices <- 1:dim(ED)[1]

M <- length(ED[(!is.na(ED[,3])) & (is.na(ED[,4])),1])

freq <- matrix(0, 2, 7)  #Row 1 no. of accepted proposals, row 2 no. of proposals

proposal.probs <- (1:4)/4 #Cumulative proposal probabilities for each reversible move (single birth/death : pair birth/death : merge/split : block recolour)

root.node <- which(is.na(ED[,2]))
non.root.nodes <- ED[ED[,1] != root.node, 1]
tree.length <- 0
for (i in non.root.nodes){
  i.row <- which(ED[,1] == i)
  parent.row <- which(ED[,1] == ED[i.row, 2])
  tree.length <- tree.length + (ED[i.row, 6] - ED[parent.row, 6])
}


n.stored.samples <- min(1e6, N)
Samples.to.store <- round(seq.int(1, N, length.out = n.stored.samples))
Sample <- numeric(n.stored.samples)
Count <- 1

prior.ratio <- lambda / ((n.deme - 1) * tree.length)

for (i in -N0 : N){
  U <- runif(1)
  V <- runif(1)
  W <- runif(1)

  if (U < proposal.probs[1]){
    if (V < 0.5){
      proposal <- ed.mig.birth(ED, n.deme, TRUE, node.indices)
      accept.prob <- min(1, lambda/(M+1))
      which.move <- 1
    } else{
      proposal <- ed.mig.death(ED, n.deme, TRUE, node.indices)
      accept.prob <- min(1, M / lambda)
      which.move <- 2
    }
  } else if (U < proposal.probs[2]){
    if (V < 0.5){
      proposal <- ed.mig.pair.birth(ED, n.deme, node.indices)
      accept.prob <- min(1, proposal$prop.ratio * prior.ratio^2)
      which.move <- 3
    } else{
      proposal <- ed.mig.pair.death(ED, n.deme, node.indices)
      accept.prob <- min(1, proposal$prop.ratio * prior.ratio^(-2))
      which.move <- 4
    }
  } else if (U < proposal.probs[3]){
    if (V < 0.5){
      proposal <- ed.coal.split(ED, n.deme, node.indices)
      proposal.M <- length(proposal$ED[(!is.na(proposal$ED[,3])) & is.na(proposal$ED[,4]), 1])
      dM <- proposal.M - M
      accept.prob <- min(1, proposal$prop.ratio * prior.ratio^dM)
      which.move <- 5
    } else{
      proposal <- ed.coal.merge(ED, n.deme, node.indices)
      proposal.M <- length(proposal$ED[(!is.na(proposal$ED[,3])) & is.na(proposal$ED[,4]), 1])
      dM <- proposal.M - M
      accept.prob <- min(1, proposal$prop.ratio * prior.ratio^dM)
      which.move <- 6
    }
  } else{
    proposal <- ed.block.recolour(ED, n.deme, TRUE, node.indices)
    which.move <- 7
    if (proposal$prop.ratio > 0){
      accept.prob <- 1
    } else{
      accept.prob <- 0
    }
  }

  if (W <= accept.prob){
    freq[1, which.move] <- freq[1, which.move] + 1
    ED <- proposal$ED
    node.indices <- proposal$node.indices
    M <- length(ED[(!is.na(ED[,3])) & is.na(ED[,4]), 1])

  }

  freq[2, which.move] <- freq[2, which.move] + 1

  if (i > 0){
    Sample[i] <- M
    #Count <- Count + 1
  }

  if (i %in% floor(0:100 * ((N+N0)/100))){
    print(i * 100 / (N+N0))
  }
}

lower <- min(Sample)
upper <- max(Sample)

hist(Sample, freq = FALSE, breaks = (lower:(upper + 1)) - 0.5)
lines(dpois(0:100, lambda), col = "red", lwd = 2, lty = 2)

png("C:/Users/ian_p/Desktop/18-MonthPanelFigs/MoveValidation.png", width = 16, height = 12, units = "cm", res = 144)
  layout(matrix(1:2, 1, 2))
  plot(Sample, type = 'l', main = "Trace Plot")
  hist(Sample, freq = FALSE, main = "Observed number of migrations", breaks = (lower:(upper + 1)) - 0.5)
  lines(dpois(0:100, lambda), col = "red", lwd = 2, lty = 2)
dev.off()
