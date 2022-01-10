require(ape)
require(magick)
#devtools::load_all()

##### Testing full MCMC algorithm using structured coalescent likelihood
set.seed(10)
n <- 10
n.deme <- 3

data <- matrix(0,nrow = n, ncol = 3)
data[,1] <- 1:n
data[,2] <- rep(2022, n) #runif(n, min = 2010, max = 2022)
data[,3] <- rep(1, n) #sample.int(n.deme, n, replace = TRUE)


#Prior parameters
eff.pop.alpha <- 10
eff.pop.beta <- 2

mig.mat.alpha <- 1
mig.mat.beta <- 10

#Simulation parameters (drawn from priors)
gen.length <- 1
effective.pop <- rep(1, n.deme)  #1 / rgamma(n.deme, shape = eff.pop.alpha, scale = 1 /eff.pop.beta)
sim.migration.matrix <- matrix(1/10, n.deme, n.deme) #matrix(rgamma(n.deme^2, shape = mig.mat.alpha, scale = 1/mig.mat.beta), n.deme, n.deme)
diag(sim.migration.matrix) <- 0

migration.matrix <- sim.migration.matrix

#phylo <- Structured.sim(data, effective.pop, gen.length, n.deme, migration.matrix, FALSE)
phylo <- ape::rcoal(n)
phylo$node.deme <- rep(1, 2*n-1)

ED <- phylo.to.ed(phylo)
ED.like <- ed.likelihood(ED, effective.pop, gen.length, migration.matrix)$log.likelihood

N0 <- 0 #1e3
N <- 1e5

freq <- matrix(0, 2, 9)  #Row 1 no. of accepted proposals, row 2 no. of proposals

coalescence.nodes <- ED[!is.na(ED[,4]), 1]
n.coal <- length(coalescence.nodes)
coal.node.deme.freq <- matrix(0, n.coal, n.deme, dimnames = list(coalescence.nodes, NULL)) #Matrix to store freq of deme at each coalescence node during iteration

proposal.rates <- c(1, 4, 4, 1, 0.1, 0.1) #Relative rates of selecting each type of proposal mechanism
proposal.probs <- cumsum(proposal.rates/sum(proposal.rates)) #Cumulative proposal probabilities for each reversible move (single birth/death : pair birth/death : merge/split : block recolour)

max.posterior.sample <- list(ED = ED, effective.population = effective.pop, migration.matrix = migration.matrix, iteration = -N0, log.likelihood = ED.like)

######### Gif output of MCMC run initialisation

k <- 100  #Number of trees to store throughout MCMC run
n.zero <- floor(log10(k))
png(paste0("VideoOut/File", sprintf("%04d", 0),".png"))
structured.plot(ed.to.phylo(ED), n.deme)
dev.off()
video.count <- 1

for (i in -N0 : N){
#for (i in -N0 : (-938)){
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
    effective.pop <- eff.pop.update(ED, effective.pop, n.deme, 1, 1)  #Update as Gibbs move for effective population sizes with update parameters alpha = beta = 1
    which.move <- 8
    ED.like <- ed.likelihood(ED, effective.pop, gen.length, migration.matrix)$log.likelihood
  } else if (U < proposal.probs[6]){
    migration.matrix <- mig.rate.update(ED, migration.matrix, n.deme, 1, 10)  #Update as Gibbs move for migration rates with update parameters alpha = 1, beta = 10
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

  if (ED.like > max.posterior.sample$log.likelihood){
    max.posterior.sample <- list(ED = ED, effective.population = effective.pop, migration.matrix = migration.matrix, iteration = i, log.likelihood = ED.like)
  }

  if (i > 0){ #Record deme at coalescent nodes for sampled tree
    coalescence.nodes <- ED[!is.na(ED[,4]), 1]
    for (j in 1:n.coal){
      coal.row <- which(ED[,1] == coalescence.nodes[j])
      coal.node.deme.freq[j, ED[coal.row, 5]] <- coal.node.deme.freq[j, ED[coal.row, 5]] + 1
    }

    #### Record effective population sizes and migration matrix
  }

  if (i %in% floor((1:(k+1)) * N/(k+1))){  #Frames for gif output, record k plots spread evenly over the entire run
    png(paste0("VideoOut/File", sprintf("%04d", video.count), ".png"))
    structured.plot(ed.to.phylo(ED), n.deme)
    dev.off()

    video.count <- video.count + 1
  }

  #Progress bar
  if (i %in% floor(0:100 * ((N+N0)/100))){
    print(i * 100 / (N+N0))
  }
}

beepr::beep()  #Laptop dings on completion...

structured.plot(phylo, n.deme)  #Initial tree

#Convert ED to phylo and back to allow coalescence node labels to be identified
ED <- phylo.to.ed(ed.to.phylo(ED))
coalescence.nodes <- ED[!is.na(ED[,4]), 1]
node <- rep(0, n.coal)
for (j in 1 : n.coal){
  node[j] <- which(ED[,1] == coalescence.nodes[j])
}

#Produce gif of saved phylo plots
unlink("VideoOut/Vid.gif")  #Remove previous version of video
img_list <- sapply(paste0("VideoOut/File", sprintf("%04d", 0:(k+1)), ".png"), image_read)
img_joined <- image_join(img_list)
img_animated <- image_animate(img_joined, fps = 2)
image_write(image = img_animated, path = "VideoOut/Vid.gif")
file.show("VideoOut/Vid.gif")
#sapply(paste0("VideoOut/File", sprintf("%04d", 0:(video.count - 1)), ".png"), unlink)  #Wipe stored phylo plots



#Final sampled tree with coalescence node deme pie charts
structured.plot(ed.to.phylo(ED), n.deme)
color.palette <- rainbow(n.deme)
coal.node.rows <- numeric(n.coal)
coalescence.nodes <- ED[!is.na(ED[,4]), 1]
for (i in 1 : n.coal){
  coal.node.rows[i] <- which(ED[,1] == coalescence.nodes[i])
}
nodelabels(node = coal.node.rows,pie = coal.node.deme.freq/N, piecol = color.palette, cex = 0.75)
