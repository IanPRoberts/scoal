require(ape)
require(magick)
devtools::load_all()

set.seed(1)
N0 <- 1e4  #Burn in
N <- 1e5  #Main MCMC run

n <- 50
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

n.stored.samples <- 1e4
mig.eff.pop.sample <- array(0, c(n.deme, n.deme, n.stored.samples))
ED.sample <- list()
samples.to.store <- round(seq.int(1, N, length.out = n.stored.samples))
sample.count <- 1

coalescence.nodes <- ED[!is.na(ED[,4]), 1]
n.coal <- length(coalescence.nodes)
coal.node.deme.freq <- matrix(0, n.coal, n.deme, dimnames = list(coalescence.nodes, NULL)) #Matrix to store freq of deme at each coalescence node during iteration

mig.mat.prior <- dgamma(migration.matrix, shape = mig.prior.shape, rate = mig.prior.rate, log = TRUE)
diag(mig.mat.prior) <- 0
eff.pop.prior <- dgamma(1/effective.pop, shape = eff.pop.prior.shape, rate = eff.pop.prior.rate, log = TRUE)
log.prior <- sum(eff.pop.prior) + sum(mig.mat.prior)

max.posterior.sample <- list(ED = ED, effective.population = effective.pop, migration.matrix = migration.matrix, iteration = -N0, log.likelihood = ED.like, log.posterior = ED.like + log.prior)

######### Gif output of MCMC run initialisation

#Create folders to store images
new.directory <- file.path("./MCMC_Results", format(Sys.time(), "%F %H-%M"))
new.directory <- gsub(" ", "_", new.directory)
dir.create(new.directory)  #Create directory to store plots; directory name gives date and time
dir.create(paste0(new.directory,"/GifOut"))

#Initialise first image
k <- 100  #Number of trees to store throughout MCMC run
n.zero <- floor(log10(k))
png(paste0(new.directory,"/GifOut/File", sprintf("%04d", 0),".png"))
structured.plot(ed.to.phylo(ED), n.deme)
dev.off()
video.count <- 1

for (i in -N0 : N){
  U <- runif(1)
  V <- runif(1)
  W <- runif(1)

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
    eff.pop.prior <- dgamma(1/effective.pop, shape = eff.pop.prior.shape, rate = eff.pop.prior.rate, log = TRUE)
    ED.like <- ed.likelihood(ED, effective.pop, gen.length, migration.matrix, node.indices)$log.likelihood
  } else if (U < proposal.probs[6]){
    which.move <- 9
    migration.matrix <- mig.rate.update(ED, migration.matrix, n.deme, node.indices, shape = mig.prior.shape, rate = mig.prior.rate)
    mig.mat.prior <- dgamma(migration.matrix, shape = mig.prior.shape, rate = mig.prior.rate, log = TRUE)
    diag(mig.mat.prior) <- 0
    ED.like <- ed.likelihood(ED, effective.pop, gen.length, migration.matrix, node.indices)$log.likelihood
  }

  freq[2, which.move] <- freq[2, which.move] + 1

  if ((which.move <= 7) && (proposal$prop.ratio > 0)){
    proposal.like <- ed.likelihood(proposal$ED, effective.pop, gen.length, migration.matrix, proposal$node.indices)$log.likelihood
    log.accept.prob <- min(0, proposal.like - ED.like + proposal$log.prop.ratio)
    if (log(W) <= log.accept.prob){
      freq[1, which.move] <- freq[1, which.move] + 1
      ED <- proposal$ED
      ED.like <- proposal.like
      node.indices <- proposal$node.indices
    }
  }

  log.posterior <- sum(eff.pop.prior) + sum(mig.mat.prior) + ED.like
  if (log.posterior >= max.posterior.sample$log.posterior){
    max.posterior.sample <- list(ED = ED, effective.population = effective.pop, migration.matrix = migration.matrix, iteration = i, log.likelihood = ED.like, log.posterior = log.posterior)
  }

  if (i > 0){
    #Record deme at coalescent nodes for sampled tree
    coalescence.nodes <- ED[!is.na(ED[,4]), 1]
    for (j in 1:n.coal){
      coal.row <- which(ED[,1] == coalescence.nodes[j])
      coal.node.deme.freq[j, ED[coal.row, 5]] <- coal.node.deme.freq[j, ED[coal.row, 5]] + 1
    }
  }

  if (i %in% samples.to.store){
    ED.sample[[sample.count]] <- ED
    mig.eff.pop.sample[,,sample.count] <- migration.matrix
    diag(mig.eff.pop.sample[,,sample.count]) <- effective.pop
    sample.count <- sample.count + 1
  }

  if (i %in% floor((1:(k+1)) * N/(k+1))){  #Frames for gif output, record k plots spread evenly over the entire run
    png(paste0(new.directory,"/GifOut/File", sprintf("%04d", video.count), ".png"))
    structured.plot(ed.to.phylo(ED), n.deme)
    dev.off()

    video.count <- video.count + 1
  }

  #Progress bar
  if ((i + N0) %in% floor(0:100 * ((N+N0)/100))){
    print((i+N0) * 100 / (N+N0))
  }
}

img_list <- sapply(paste0(new.directory,"/GifOut/File", sprintf("%04d", 0:(k+1)), ".png"), image_read)
img_joined <- image_join(img_list)
img_animated <- image_animate(img_joined, fps = 2)
image_write(image = img_animated, path = paste0(new.directory,"/Vid.gif"))

# Plot matrix plot of histograms for migration rates and effective population sizes

png(paste0(new.directory, "/mig.eff.pop.sample.hist.png"), width = 2000, height = 1500)
layout(matrix(1:n.deme^2, n.deme, n.deme, byrow = TRUE))
for (i in 1:n.deme){
  for (j in 1:n.deme){
    if (i == j){
      upper <- ceiling(5*max(mig.eff.pop.sample[i, j, ]))/5
      hist(mig.eff.pop.sample[i, j, ], freq = FALSE, main = bquote(theta[.(i)]), xlab = bquote(theta[.(i)]), breaks = seq(0, upper, 1/5))
    } else{
      upper <- ceiling(2*max(mig.eff.pop.sample[i, j, ]))/2
      hist(mig.eff.pop.sample[i, j, ], freq = FALSE, main = bquote(lambda[paste("(", .(i), ",", .(j), ")")]), xlab = bquote(lambda[paste("(", .(i), ",", .(j), ")")]), breaks = seq(0, upper, 1/10), xlim = c(0,2))
    }
  }
}
dev.off()

# Matrix plot for trace plots of migration rates and effective population sizes
png(paste0(new.directory, "/mig.eff.pop.sample.trace.png"), width = 2000, height = 1500)
layout(matrix(1:n.deme^2, n.deme, n.deme, byrow = TRUE))
for (i in 1:n.deme){
  for (j in 1:n.deme){
    if (i == j){
      plot(mig.eff.pop.sample[i,j,], type = 'l', main = bquote(paste(theta[.(i)], ", mean =", .(round(mean(mig.eff.pop.sample[i,j,]), 4)))), ylab = bquote(theta[.(i)]) )
    } else{
      plot(mig.eff.pop.sample[i,j,], type = 'l', main = bquote(paste(lambda[paste("(", .(i), ",", .(j), ")")], ", mean =", .(round(mean(mig.eff.pop.sample[i,j,]), 4)))), ylab = bquote(lambda[paste("(", .(i), ",", .(j), ")")]))
    }
  }
}
dev.off()

# Matrix plot for histogram and trace plot of number of migrations

n.mig.sample <- numeric(n.stored.samples)
for (i in 1 : n.stored.samples){
  n.mig.sample[i] <- dim(ED.sample[[i]])[1] - 2*n + 1
}
png(paste0(new.directory, "/N.mig.plots.png"), width = 2000, height = 1500)
layout(matrix(1:2, 1, 2))
plot(n.mig.sample, type = 'l', main = "Trace Plot", xlab = "N.mig")
hist(n.mig.sample, freq = FALSE, main = "Observed number of migrations", xlab = "N.mig", breaks = 0:max(n.mig.sample + 1) - 0.5)
dev.off()


#Maximum posterior sampled tree with coalescence node deme pie charts
png(paste0(new.directory, "/max.post.sample.png"), width = 2000, height = 1500)
structured.plot(max.posterior.sample$ED, n.deme)
color.palette <- rainbow(n.deme)
coal.node.rows <- numeric(n.coal)
coalescence.nodes <- ED[!is.na(ED[,4]), 1]
for (i in 1 : n.coal){
  coal.node.rows[i] <- which(ED[,1] == coalescence.nodes[i])
}
nodelabels(node = coal.node.rows,pie = coal.node.deme.freq/N, piecol = color.palette, cex = 0.75)
dev.off()
