#set.seed(10)

n <- 10
n.deme <- 3

phylo <- Homochronous.Sim(1:n, 1, 1)
phylo$node.deme <- rep(1, 2 * n - 1)
#phylo <- ed.to.phylo(ED)

ED <- phylo.to.ed(phylo)
M <- length(ED[(!is.na(ED[,3])) & (is.na(ED[,4])),1])

N <- 1e5
N0 <- 1e3

lambda <- 10

freq <- matrix(0, 2, 7)  #Row 1 no. of accepted proposals, row 2 no. of proposals
M.freq = matrix(c(0:150, rep(0, 151)), 2, 151, byrow = TRUE)

proposal.probs <- c(0.1,0.5, 0.9, 1)  #Cumulative proposal probabilities for each reversible move (single birth/death : pair birth/death : merge/split : block recolour)

root.node <- which(is.na(ED[,2]))
non.root.nodes <- ED[ED[,1] != root.node, 1]
tree.length <- 0
for (i in non.root.nodes){
  i.row <- which(ED[,1] == i)
  parent.row <- which(ED[,1] == ED[i.row, 2])
  tree.length <- tree.length + (ED[i.row, 6] - ED[parent.row, 6])
}

Sample <- numeric(N)
prior.ratio <- lambda / ((n.deme - 1) * tree.length)

for (i in -N0 : N){
  U <- runif(1)
  V <- runif(1)
  W <- runif(1)

  if (U < proposal.probs[1]){
    if (V < 0.5){
      proposal <- ed.mig.birth.4(ED, n.deme, FALSE)
      accept.prob <- min(1, lambda/(M+1))
      which.move <- 1
    } else{
      proposal <- ed.mig.death.4(ED, n.deme, FALSE)
      accept.prob <- min(1, M / lambda)
      which.move <- 2
    }

    if (W <= accept.prob){
      ED <- proposal$ED
      M <- length(ED[(!is.na(ED[,3])) & (is.na(ED[,4])),1])
      if (i > 0){
        freq[1, which.move] <- freq[1, which.move] + 1
        M.freq[2, M + 1] <- M.freq[2, M + 1] + 1
      }
    }  else if (i > 0){
      M.freq[2, M + 1] <- M.freq[2, M + 1] + 1
    }
  } else if (U < proposal.probs[2]){
    if (V < 0.5){
      proposal <- ed.mig.pair.birth(ED, n.deme) #ed.pair.birth.3(ED, n.deme) #ed.mig.pair.birth(ED, n.deme)
      accept.prob <- min(1, proposal$prop.ratio * prior.ratio^2)
      which.move <- 3
    } else{
      proposal <- ed.mig.pair.death(ED, n.deme) #ed.pair.death.3(ED, n.deme) #ed.mig.pair.death(ED, n.deme)
      accept.prob <- min(1, proposal$prop.ratio * prior.ratio^(-2))
      which.move <- 4
    }

    if (W <= accept.prob){
      ED <- proposal$ED
      M <- length(ED[(!is.na(ED[,3])) & is.na(ED[,4]), 1])
      if (i > 0){
        freq[1, which.move] <- freq[1, which.move] + 1
        M.freq[2, M + 1] <- M.freq[2, M + 1] + 1
      }
    } else if (i > 0){
      M.freq[2, M + 1] <- M.freq[2, M + 1] + 1
    }
  } else if (U < proposal.probs[3]){
    if (V < 0.5){
      proposal <- ed.coal.split(ED, n.deme)
      proposal.M <- length(proposal$ED[(!is.na(proposal$ED[,3])) & is.na(proposal$ED[,4]), 1])
      dM <- proposal.M - M
      which.move <- 5
    } else{
      proposal <- ed.coal.merge(ED, n.deme)
      proposal.M <- length(proposal$ED[(!is.na(proposal$ED[,3])) & is.na(proposal$ED[,4]), 1])
      dM <- proposal.M - M
      which.move <- 6
    }
    accept.prob <- min(1, proposal$prop.ratio * (lambda / ((n.deme - 1) * tree.length))^dM)

    if (W <= accept.prob){
      ED <- proposal$ED
      M <- proposal.M
      if (i > 0){
        freq[1, which.move] <- freq[1, which.move] + 1
        M.freq[2, M + 1] <- M.freq[2, M + 1] + 1
      }
    } else if (i > 0){
      M.freq[2, M + 1] <- M.freq[2, M + 1] + 1
    }
  } else{
    proposal <- ed.block.recolour(ED, n.deme, FALSE)
    which.move <- 7
    if (proposal$prop.ratio > 0){
      ED <- proposal$ED
      if (i > 0){
        freq[1, which.move] <- freq[1, which.move] + 1
      }
    }
    if (i > 0){
      M.freq[2, M+1] <- M.freq[2, M+1] + 1
    }
  }

  freq[2, which.move] <- freq[2, which.move] + 1

  if (i > 0){
    Sample[i] <- M
  }

  if (i %in% floor(0:100 * ((N+N0)/100))){
    print(i * 100 / (N+N0))
  }
}

lower <- min(which(M.freq[2,] > 0))
upper <- lower + max(which(M.freq[2, (lower+1):151] > 0))
plot(M.freq[1,lower:upper], M.freq[2,lower:upper]/N, type = 'l', xlab = "M", ylab = "Density", main = paste(n, "leaves,", n.deme, "demes, lambda =", lambda, ", ", N, "iterations"))
lines(0:150, dpois(0:150, lambda), lty = 2, col = "red")

plot(Sample, type = 'l')

hist(Sample, breaks = (min(Sample)-0.5):(max(Sample)+1), probability = TRUE)
lines(M.freq[1,lower:upper], M.freq[2,lower:upper]/N, type = 'l')
lines(0:150, dpois(0:150, lambda), lty = 2, col = "red")

beepr::beep()  #Laptop dings on completion...
