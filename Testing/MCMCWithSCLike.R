set.seed(1)

n <- 10  #No. leaves in simulation
n.deme <- 3

#Simulated data for initial tree
data <- matrix(0,nrow = n, ncol = 3)
data[,1] <- 1:n
data[,2] <- runif(n, min = 2010, max = 2022)
data[,3] <- sample.int(n.deme, n, replace = TRUE)

#Prior parameters
eff.pop.alpha <- 10
eff.pop.beta <- 2

mig.mat.alpha <- 1
mig.mat.beta <- 10

#Draw gen.length, effective.pop and migration.matrix from priors
gen.length <- 1
effective.pop <- 1 / rgamma(n.deme, shape = eff.pop.alpha, scale = 1 / eff.pop.beta)
migration.matrix <- matrix(rgamma(n.deme^2, shape = mig.mat.alpha, scale = 1/mig.mat.beta), n.deme, n.deme)
diag(migration.matrix) <- 0

phylo <- Structured.sim(data, effective.pop, gen.length, n.deme, migration.matrix, TRUE) #Initial tree; "TRUTH"

ED <- phylo.to.ed(phylo)
ED.like <- ed.likelihood(ED, effective.pop, gen.length, migration.matrix)

max.post.like <- list(ED = ED, log.likelihood = ED.like$log.likelihood, likelihood = ED.like$likelihood, n = 0)

N <- 1e7
N0 <- 1e3

freq <- matrix(0, 2, 7)  #Row 1 no. of accepted proposals, row 2 no. of proposals

proposal.probs <- c(0.25,0.5, 0.75, 1)  #Cumulative proposal probabilities for each reversible move (single birth/death : pair birth/death : merge/split : block recolour)

for (i in -N0 : N){
  U <- runif(1)
  V <- runif(1)
  W <- runif(1)

  if (U < proposal.probs[1]){ #Single birth/death proposal
    if (V < 0.5){
      proposal <- ed.mig.birth.4(ED, n.deme, TRUE)
      which.move <- 1
    } else{
      proposal <- ed.mig.death.4(ED, n.deme, TRUE)
      which.move <- 2
    }
  } else if (U < proposal.probs[2]){ #Pair birth/death proposal
    if (V < 0.5){
      proposal <- ed.mig.pair.birth(ED, n.deme)
      which.move <- 3
    } else{
      proposal <- ed.mig.pair.death(ED, n.deme)
      which.move <- 4
    }
  } else if (U < proposal.probs[3]){ #Migration merge/split proposal
    if (V < 0.5){
      proposal <- ed.coal.split(ED, n.deme)
      which.move <- 5
    } else{
      proposal <- ed.coal.merge(ED, n.deme)
      which.move <- 6
    }
  } else if (U < proposal.probs[4]){ #Block recolour proposal
    proposal <- ed.block.recolour(ED, n.deme, TRUE)
    which.move <- 7
  }

  freq[2, which.move] <- freq[2, which.move] + 1

  if (proposal$prop.ratio > 0){
    proposal.like <- ed.likelihood(proposal$ED, effective.pop, gen.length, migration.matrix)
    log.accept.prob <- min(0, proposal.like$log.likelihood - ED.like$log.likelihood + proposal$log.prop.ratio)
    if (log(W) <= log.accept.prob){
      freq[1, which.move] <- freq[1, which.move] + 1
      ED <- proposal$ED
      ED.like <- proposal.like
    }

      if (ED.like$log.likelihood > max.post.like$log.likelihood){
        max.post.like <- list(ED = ED, log.likelihood = ED.like$log.likelihood, likelihood = ED.like$likelihood, n = i)
      }
  }

  #"Progress bar"
  if (i %in% floor(0:100 * ((N+N0)/100))){
    print(i * 100 / (N+N0))
  }
}

beepr::beep()  #Laptop dings on completion...
structured.plot(ed.to.phylo(ED))
cbind(ED[1:10,], phylo.to.ed(phylo)[1:10,])
