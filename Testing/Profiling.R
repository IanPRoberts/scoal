require(profvis)
require(ape)
devtools::load_all()

profvis({
  set.seed(1)
  N0 <- 1e4  #Burn in
  N <- 1e5  #Main MCMC run

  n <- 20
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

    ED.like <- ed.likelihood(ED, effective.pop, gen.length, migration.matrix)$log.likelihood
    freq <- matrix(0, 2, 9)  #Row 1 no. of accepted proposals, row 2 no. of proposals

    proposal.probs <- cumsum(proposal.rates/sum(proposal.rates)) #Cumulative proposal probabilities for each reversible move (single birth/death : pair birth/death : merge/split : block recolour)


    for (i in -N0 : N){
      U <- runif(1)
      V <- runif(1)
      W <- runif(1)

      if (U < proposal.probs[1]){
        if (V < 0.5){
          proposal <- ed.mig.birth.4(ED, n.deme, TRUE)
          which.move <- 1
        } else{
          proposal <- ed.mig.death.4(ED, n.deme, TRUE)
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
        proposal <- ed.block.recolour(ED, n.deme, TRUE)
        which.move <- 7
      } else if (U < proposal.probs[5]){
        effective.pop <- eff.pop.update.2(ED, effective.pop, n.deme, shape = eff.pop.prior.shape, rate = eff.pop.prior.rate)
        which.move <- 8
        ED.like <- ed.likelihood(ED, effective.pop, gen.length, migration.matrix)$log.likelihood
      } else if (U < proposal.probs[6]){
        migration.matrix <- mig.rate.update.2(ED, migration.matrix, n.deme, shape = mig.prior.shape, rate = mig.prior.rate)
        which.move <- 9
        ED.like <- ed.likelihood(ED, effective.pop, gen.length, migration.matrix)$log.likelihood
      }

      freq[2, which.move] <- freq[2, which.move] + 1

      if ((which.move <= 7) && (proposal$prop.ratio > 0)){
        proposal.like <- ed.likelihood(proposal$ED, effective.pop, gen.length, migration.matrix)$log.likelihood
        log.accept.prob <- min(0, proposal.like - ED.like + proposal$log.prop.ratio)
        if (log(W) <= log.accept.prob){
          freq[1, which.move] <- freq[1, which.move] + 1
          ED <- proposal$ED
          ED.like <- proposal.like
        }
      }

      #Progress bar
      if ((i + N0) %in% floor(0:100 * ((N+N0)/100))){
        print((i+N0) * 100 / (N+N0))
      }
    }
})
