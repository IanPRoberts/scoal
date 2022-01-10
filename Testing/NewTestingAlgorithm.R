##### Testing full MCMC algorithm using "dumb" likelihood - Poisson marginal on no. of migrations using Poi-likelihood ratio and proposal ratio instead of explicitly calculated versions

log.like <- function(ED, rate = 10, n.deme = NA){  #log-likelihood for Poisson marginal dist (up to additive const)
  root.node <- ED[is.na(ED[,2]), 1]
  coalescence.nodes <- ED[!is.na(ED[,4]),1]
  coalescence.nodes <- coalescence.nodes[coalescence.nodes != root.node]
  leaf.nodes <- ED[(is.na(ED[,3])) & (is.na(ED[,4])), 1]
  migration.nodes <- ED[ is.na(ED[,4]) & (!is.na(ED[,3])) ,1]

  edge.length <- numeric(dim(ED)[1])
  non.root.nodes <- c(coalescence.nodes, leaf.nodes, migration.nodes)
  non.root.nodes <- non.root.nodes[!is.na(non.root.nodes)]
  for (j in non.root.nodes){
    node.row <- which(ED[,1] == j)
    parent.row <- which(ED[,1] == ED[node.row, 2])
    edge.length[node.row] <- ED[node.row, 6] - ED[parent.row, 6]
  }

  tree.length <- sum(edge.length)  #Total tree length

  M <- length(migration.nodes)
  out <- M * (log(rate) - log(tree.length) - log(n.deme - 1))
}

n <- 10
n.deme <- 3

data <- matrix(0,nrow = n, ncol = 3)
data[,1] <- 1:n
data[,2] <- runif(n, min = 2010, max = 2022)
data[,3] <- sample.int(n.deme, n, replace = TRUE)


#Prior parameters
eff.pop.alpha <- 10
eff.pop.beta <- 2

mig.mat.alpha <- 1
mig.mat.beta <- 10

#Simulation parameters (drawn from priors)
gen.length <- 1
effective.pop <- 1 / rgamma(n.deme, shape = eff.pop.alpha, scale = 1 / eff.pop.beta)
migration.matrix <- matrix(rgamma(n.deme^2, shape = mig.mat.alpha, scale = 1/mig.mat.beta), n.deme, n.deme)
diag(migration.matrix) <- 0

phylo <- Structured.sim(data, effective.pop, gen.length, n.deme, migration.matrix, FALSE)

ED <- phylo.to.ed(phylo)
ED.like <- log.like(ED, lambda, n.deme)

Sample <- numeric(N)

N <- 1e6
N0 <- 1e3

freq <- matrix(0, 2, 7)  #Row 1 no. of accepted proposals, row 2 no. of proposals
M.freq = matrix(c(0:150, rep(0, 151)), 2, 151, byrow = TRUE)

migration.nodes <- ED[ is.na(ED[,4]) & (!is.na(ED[,3])) ,1]
M <- length(migration.nodes)


proposal.rates <- c(1, 4, 4, 1) #Relative rates of selecting each type of proposal mechanism
proposal.probs <- cumsum(proposal.rates/sum(proposal.rates)) #Cumulative proposal probabilities for each reversible move (single birth/death : pair birth/death : merge/split : block recolour)

lambda <- 50

for (i in -N0 : N){
  U <- runif(1)
  V <- runif(1)
  W <- runif(1)

  if (U < proposal.probs[1]){
    if (V < 0.5){
      proposal <- ed.mig.birth.4(ED, n.deme, FALSE)
      which.move <- 1
    } else{
      proposal <- ed.mig.death.4(ED, n.deme, FALSE)
      which.move <- 2
    }
  } else if (U < proposal.probs[2]){
    if (V < 0.5){
      proposal <- ed.mig.pair.birth(ED, n.deme)
      which.move <- 3
    } else{
      proposal <- ed.mig.pair.death(ED, n.deme)
      which.move <- 4
    }
  } else if (U < proposal.probs[3]){
    if (V < 0.5){
      proposal <- ed.coal.split(ED, n.deme)
      which.move <- 5
    } else{
      proposal <- ed.coal.merge(ED, n.deme)
      which.move <- 6
    }
  } else if (U < proposal.probs[4]){
    proposal <- ed.block.recolour(ED, n.deme, FALSE)
    which.move <- 7
  }

  freq[2, which.move] <- freq[2, which.move] + 1

  if (proposal$prop.ratio > 0){
    proposal.like <- log.like(proposal$ED, lambda, n.deme)
    log.accept.prob <- min(0, proposal.like - ED.like + proposal$log.prop.ratio)
    if (log(W) <= log.accept.prob){
      freq[1, which.move] <- freq[1, which.move] + 1
      ED <- proposal$ED
      ED.like <- proposal.like

      migration.nodes <- ED[ is.na(ED[,4]) & (!is.na(ED[,3])) ,1]
      M <- length(migration.nodes)
    }
  }

  if (i > 0){
    Sample[i] <- M
    M.freq[2, M+1] <- M.freq[2, M+1] + 1
  }

  #"Progress bar"
  if (i %in% floor(0:100 * ((N+N0)/100))){
    print(i * 100 / (N+N0))
  }
}

beepr::beep()  #Laptop dings on completion...

lower <- min(which(M.freq[2,] > 0))
upper <- lower + max(which(M.freq[2, (lower+1):151] > 0))
plot(M.freq[1,lower:upper], M.freq[2,lower:upper]/N, type = 'l', xlab = "M", ylab = "Density", main = paste(n, "leaves,", n.deme, "demes, lambda =", lambda, ", ", N, "iterations"))
lines(0:150, dpois(0:150, lambda), lty = 2, col = "red")

#plot(Sample, type = 'l')

hist(Sample, breaks = (min(Sample)-0.5):(max(Sample)+1), probability = TRUE)
lines(M.freq[1,lower:upper], M.freq[2,lower:upper]/N, type = 'l')
lines(0:150, dpois(0:150, lambda), lty = 2, col = "red")
