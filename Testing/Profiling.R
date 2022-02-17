require(profvis)
require(ape)
devtools::load_all()

profvis({
  max.move.fixed <- 6

  set.seed(1)
  N0 <- 1e4  #Burn in
  N <- 1e5  #Main MCMC run

  n <- 20
  n.deme <- 3
  proposal.rates <- c(rep(5, 4), 1, 1) #c(1, 4, 4, 1, 0.1, 0.1) #Relative rates of selecting each type of proposal mechanism

  data <- matrix(0,nrow = n, ncol = 3)
  data[,1] <- 1:n
  data[,2] <- rep(2022, n) #runif(n, min = 2010, max = 2022)
  data[,3] <- sample.int(n.deme, n, replace = TRUE)

  #Prior Parameters
  eff.pop.prior.mean <- 1
  eff.pop.prior.var <- 1
  mig.prior.mean <- 1/10
  mig.prior.var <- 1/20

  eff.pop.prior.shape <- eff.pop.prior.mean^2/eff.pop.prior.var + 2
  eff.pop.prior.rate <- eff.pop.prior.mean * (eff.pop.prior.mean^2/eff.pop.prior.var + 1)
  mig.prior.shape <- mig.prior.mean^2/mig.prior.var
  mig.prior.rate <- mig.prior.mean/mig.prior.var

  #Simulation Parameters
  n.deme <- 3
  gen.length <- 1
  effective.pop <- rep(1, n.deme)
  migration.matrix <-matrix(1, n.deme, n.deme)
  diag(migration.matrix) <- 0

  phylo <- Structured.sim(data, effective.pop, gen.length, n.deme, migration.matrix, FALSE)
  ED <- phylo.to.ed(phylo)

  max.label <- max(ED[,1])
  node.indices <- rep(0, max.label)
  for (j in 1 : dim(ED)[1]){
    node.indices[ED[j,1]] <- j
  }

  ED.like <- ed.likelihood(ED, effective.pop, gen.length, migration.matrix, node.indices)$log.likelihood
  freq <- matrix(0, 2, 9)  #Row 1 no. of accepted proposals, row 2 no. of proposals

  proposal.probs <- cumsum(proposal.rates/sum(proposal.rates)) #Cumulative proposal probabilities for each reversible move (single birth/death : pair birth/death : merge/split : block recolour)

  which.move <- max.move.fixed + 1

  for (i in -N0 : N){
    U <- runif(1)
    V <- runif(1)
    W <- runif(1)


    if (which.move > max.move.fixed){
      max.label <- max(max(ED[,1]), length(node.indices))
      node.indices <- rep(0, max.label)
      for (j in 1 : dim(ED)[1]){
        node.indices[ED[j,1]] <- j
      }
    }

    if (which.move %in% 5:6){
      max.label <- max(max(ED[,1]), length(node.indices))
      temp.indices <- rep(0, max.label)
      for (j in 1 : dim(ED)[1]){
        temp.indices[ED[j,1]] <- j
      }

      if (any(temp.indices - node.indices != 0)) browser()
    }

    if (U < proposal.probs[1]){
      if (V < 0.5){
        which.move <- 1
        proposal <- ed.mig.birth(ED, n.deme, TRUE, node.indices)
      } else{
        which.move <- 2
        proposal <- ed.mig.death(ED, n.deme, TRUE, node.indices)
      }
    } else if (U < proposal.probs[2]){
      if (V < 0.5){
        which.move <- 3
        proposal <- ed.mig.pair.birth(ED, n.deme, node.indices)
      } else{
        which.move <- 4
        proposal <- ed.mig.pair.death(ED, n.deme, node.indices)
      }
    } else if (U < proposal.probs[3]){
      if (V < 0.5){
        which.move <- 5
        proposal <- ed.coal.split(ED, n.deme, node.indices)
      } else{
        which.move <- 6
        proposal <- ed.coal.merge(ED, n.deme, node.indices)
      }
    } else if (U < proposal.probs[4]){
      which.move <- 7
      proposal <- ed.block.recolour(ED, n.deme, TRUE, node.indices)
    } else if (U < proposal.probs[5]){
      which.move <- 8
      effective.pop <- eff.pop.update(ED, effective.pop, n.deme, node.indices, shape = eff.pop.prior.shape, rate = eff.pop.prior.rate)
      ED.like <- ed.likelihood(ED, effective.pop, gen.length, migration.matrix, node.indices)$log.likelihood
    } else if (U < proposal.probs[6]){
      which.move <- 9
      migration.matrix <- mig.rate.update(ED, migration.matrix, n.deme, node.indices, shape = mig.prior.shape, rate = mig.prior.rate)
      ED.like <- ed.likelihood(ED, effective.pop, gen.length, migration.matrix, node.indices)$log.likelihood
    }

    freq[2, which.move] <- freq[2, which.move] + 1

    if ((which.move <= 7) && (proposal$prop.ratio > 0)){
      proposal.like <- ed.likelihood(proposal$ED, effective.pop, gen.length, migration.matrix, node.indices)$log.likelihood
      log.accept.prob <- min(0, proposal.like - ED.like + proposal$log.prop.ratio)
      if (log(W) <= log.accept.prob){
        freq[1, which.move] <- freq[1, which.move] + 1
        ED <- proposal$ED
        ED.like <- proposal.like

        if (which.move <= max.move.fixed){
          node.indices <- proposal$node.indices
        }
      }
    }

    #Progress bar
    if ((i + N0) %in% floor(0:100 * ((N+N0)/100))){
      print((i+N0) * 100 / (N+N0))
    }
  }
})

