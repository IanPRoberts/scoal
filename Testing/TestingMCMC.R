devtools::load_all(".")
require(magick)
require(ape)

set.seed(100)

save_ED_sample <- TRUE #Boolean - whether to save ED_sample (Large array of large matrices = lots of memory :( )
N0 <- 1e4 #Burn in
N <- 1e6  #Main MCMC run

proposal.rates <- c(rep(5, 4), 1, 1, 1) #c(rep(5, 4), 1, 1, 1) #c(1, 4, 4, 1, 0.1, 0.1) #Relative rates of selecting each type of proposal mechanism

n <- 50
n_deme <- 3
data <- matrix(0,nrow = n, ncol = 3); data[,1] <- 1:n; data[,2] <- sample(2015:2020, n, replace = TRUE); data[,3] <- sample.int(n_deme, n, replace = TRUE)

#Simulation Parameters
time_scale <- 1
true_cr <- rep(1/20, n_deme)
true_mm <-matrix(0.01, n_deme, n_deme)
diag(true_mm) <- 0

true_phylo <- scaled_sim(data, true_cr, time_scale, n_deme, true_mm, FALSE)
true_ED <- phylo.to.ed(true_phylo)

ED <- initial.ed(strip.history(true_ED), data[,3])

#Prior Parameters
cr_prior_shape <- cr_prior_rate <- 1
mm_prior_shape <- mm_prior_rate <- 1

#Initial conditions
coal_rate <- rep(1, n_deme)
mig_mat <- matrix(1, n_deme, n_deme); diag(mig_mat) <- 0
time_scale <- 1

max_label <- max(ED[,1])
node_indices <- rep(0, max_label)
for (j in 1 : dim(ED)[1]){
  node_indices[ED[j,1]] <- j
}

ED_like <- ScaledLikelihoodC(ED, coal_rate, time_scale, mig_mat, node_indices)$log.likelihood


freq <- matrix(0, 2, 10)  #Row 1 no. of accepted proposals, row 2 no. of proposals
proposal.probs <- cumsum(proposal.rates/sum(proposal.rates)) #Cumulative proposal probabilities for each reversible move (single birth/death : pair birth/death : merge/split : block recolour)

n.stored.samples <- min(N, 1e5)
mig.eff.pop.sample <- array(0, c(n_deme, n_deme, n.stored.samples))
time_scale_sample <- numeric(n.stored.samples)
if (save_ED_sample == TRUE){
  ED_sample <- list()
}
samples.to.store <- round(seq.int(1, N, length.out = n.stored.samples))

sample.count <- 1

coalescence.nodes <- ED[!is.na(ED[,4]), 1]
n.coal <- length(coalescence.nodes)
coal.node.deme.freq <- matrix(0, n.coal, n_deme, dimnames = list(coalescence.nodes, NULL)) #Matrix to store freq of deme at each coalescence node during iteration

mig_mat_prior <- dgamma(mig_mat * time_scale, shape = mm_prior_shape, rate = mm_prior_rate, log = TRUE)
diag(mig_mat_prior) <- 0
coal_rate_prior <- dgamma(coal_rate, shape = cr_prior_shape, rate = cr_prior_rate, log = TRUE)
log.prior <- sum(coal_rate_prior) + sum(mig_mat_prior)

max.posterior.sample <- list(ED = ED, coalescence_rate = coal_rate, migration_matrix = mig_mat, iteration = -N0, log.likelihood = ED_like, log.posterior = ED_like + log.prior)

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
structured.plot(ed.to.phylo(ED), n_deme)
dev.off()
video.count <- 1

#Progress bar
pb <- txtProgressBar(min = 0, max = N0 + N, initial = 0, style = 3)

for (i in -N0 : N){
  U <- runif(1)
  V <- runif(1)
  W <- runif(1)

  if (U < proposal.probs[1]){
    if (V < 0.5){
      which.move <- 1
      proposal <- ed.mig.birth(ED, n_deme, TRUE, node_indices)
    } else{
      which.move <- 2
      proposal <- ed.mig.death(ED, n_deme, TRUE, node_indices)
    }
  } else if (U < proposal.probs[2]){
    if (V < 0.5){
      which.move <- 3
      proposal <- ed.mig.pair.birth(ED, n_deme, node_indices)
    } else{
      which.move <- 4
      proposal <- ed.mig.pair.death(ED, n_deme, node_indices)
    }
  } else if (U < proposal.probs[3]){
    if (V < 0.5){
      which.move <- 5
      proposal <- ed.coal.split(ED, n_deme, node_indices)
    } else{
      which.move <- 6
      proposal <- ed.coal.merge(ED, n_deme, node_indices)
    }
  } else if (U < proposal.probs[4]){
    which.move <- 7
    proposal <- ed.block.recolour(ED, n_deme, TRUE, node_indices)
  } else if (U < proposal.probs[5]){
    which.move <- 8
    coal_rate <- scaled_coal_rate_gibbs_update(ED, time_scale, n_deme, node_indices, cr_prior_shape, cr_prior_rate)
    coal_rate_prior <- dgamma(coal_rate, shape = cr_prior_shape, rate = cr_prior_rate, log = TRUE)
    ED_like <- ScaledLikelihoodC(ED, coal_rate, time_scale, mig_mat, node_indices)$log.likelihood
    freq[1, which.move] <- freq[1, which.move] + 1
  } else if (U < proposal.probs[6]){
    which.move <- 9
    mig_mat <- scaled_mig_mat_gibbs_update(ED, time_scale, n_deme, node_indices, mm_prior_shape, mm_prior_rate)
    mig_mat_prior <- dgamma(mig_mat, shape = mm_prior_shape, rate = mm_prior_rate, log = TRUE)
    diag(mig_mat_prior) <- 0
    ED_like <- ScaledLikelihoodC(ED, coal_rate, time_scale, mig_mat, node_indices)$log.likelihood
    freq[1, which.move] <- freq[1, which.move] + 1
  } else if (U < proposal.probs[7]){
    which.move <- 10
    time_scale <- scaled_time_scale_gibbs_update(coal_rate, mig_mat, ED, n_deme, node_indices, 0.001, 0.001)
    ED_like <- ScaledLikelihoodC(ED, coal_rate, time_scale, mig_mat, node_indices)$log.likelihood
    freq[1, which.move] <- freq[1, which.move] + 1
  }

  freq[2, which.move] <- freq[2, which.move] + 1

  if ((which.move <= 7) && (proposal$prop.ratio > 0)){
    proposal_like <- ScaledLikelihoodC(proposal$ED, coal_rate, time_scale, mig_mat, proposal$node.indices)$log.likelihood
    log.accept.prob <- min(0, proposal_like - ED_like + proposal$log.prop.ratio)
    if (log(W) <= log.accept.prob){
      freq[1, which.move] <- freq[1, which.move] + 1
      ED <- proposal$ED
      ED_like <- proposal_like
      node_indices <- proposal$node.indices
    }
  }

  log.posterior <- sum(coal_rate_prior) + sum(mig_mat_prior) + ED_like
  if (log.posterior >= max.posterior.sample$log.posterior){
    max.posterior.sample <- list(ED = ED, coalescence_rate = coal_rate, migration_matrix = mig_mat, iteration = i, log.likelihood = ED_like, log.posterior = log.posterior)
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
    if (save_ED_sample == TRUE){
      ED_sample[[sample.count]] <- ED
    }

    mig.eff.pop.sample[,,sample.count] <- mig_mat
    diag(mig.eff.pop.sample[,,sample.count]) <- coal_rate
    time_scale_sample[sample.count] <- time_scale
    sample.count <- sample.count + 1
  }

  if (i %in% floor((1:(k+1)) * N/(k+1))){  #Frames for gif output, record k plots spread evenly over the entire run
    png(paste0(new.directory,"/GifOut/File", sprintf("%04d", video.count), ".png"))
    structured.plot(ed.to.phylo(ED), n_deme)
    dev.off()

    video.count <- video.count + 1
  }

  #Progress bar
  setTxtProgressBar(pb, i + N0)
}

close(pb)

#Gif of sample progress

img_list <- sapply(paste0(new.directory,"/GifOut/File", sprintf("%04d", 0:(k+1)), ".png"), image_read)
img_joined <- image_join(img_list)
img_animated <- image_animate(img_joined, fps = 2)
image_write(image = img_animated, path = paste0(new.directory,"/Vid.gif"))

# Relative migration & coalescent rate hists

png(paste0(new.directory, "/rel.params.hist.png"), width = 2000, height = 1500)
layout(matrix(1:n_deme^2, n_deme, n_deme, byrow = TRUE))
x <- seq(0,2, 0.01)
y1 <- dgamma(x, cr_prior_shape, cr_prior_rate)
y2 <- dgamma(x, mm_prior_shape, mm_prior_rate)
for (i in 1:n_deme){
  for (j in 1:n_deme){
    if (i == j){
      hist(mig.eff.pop.sample[i, j, ], freq = FALSE, main = bquote(theta[.(i)]), xlab = bquote(theta[.(i)]))
      lines(x,y1, lty = 2, col = "red")
      legend("topright", c("Prior"), col = c("red"), lty = 2)
    } else{
      hist(mig.eff.pop.sample[i, j, ], freq = FALSE, main = bquote(lambda[paste("(", .(i), ",", .(j), ")")]), xlab = bquote(lambda[paste("(", .(i), ",", .(j), ")")]))
      lines(x,y2, lty = 2, col = "red")
      legend("topright", c("Prior"), col = c("red"), lty = 2)
    }
  }
}
dev.off()

# Relative migration and coalescent rate traces
png(paste0(new.directory, "/rel.params.traces.png"), width = 2000, height = 1500)
layout(matrix(1:n_deme^2, n_deme, n_deme, byrow = TRUE))
for (i in 1:n_deme){
  for (j in 1:n_deme){
    if (i == j){
      plot(mig.eff.pop.sample[i,j,], type = 'l', main = bquote(paste(theta[.(i)], ", mean =", .(round(mean(mig.eff.pop.sample[i,j,]), 4)))), ylab = bquote(theta[.(i)]) )
    } else{
      plot(mig.eff.pop.sample[i,j,], type = 'l', main = bquote(paste(lambda[paste("(", .(i), ",", .(j), ")")], ", mean =", .(round(mean(mig.eff.pop.sample[i,j,]), 4)))), ylab = bquote(lambda[paste("(", .(i), ",", .(j), ")")]))
    }
  }
}
dev.off()

# Number of migrations/coalescences per deme traces

n <- length(ED_sample)
m <- array(dim = c(n_deme, n_deme, n))

for (i in 1:n){
  ED <- ED_sample[[i]]
  max_label <- max(ED[,1])
  node_indices <- rep(0, max_label)
  for (j in 1 : dim(ED)[1]){
    node_indices[ED[j,1]] <- j
  }

  NC <- NodeCountC(ED, n_deme, node_indices)
  m[,,i] <- NC$m
  diag(m[,,i]) <- NC$c
}

png(paste0(new.directory, "/migrations.per.deme.png"), width = 2000, height = 1500)
layout(matrix(1:n_deme^2, n_deme, n_deme, byrow = TRUE))
for (i in 1:n_deme){
  for (j in 1:n_deme){
    if (i == j){
      plot(m[i,j,], type = 'l', ylab = bquote(theta[.(i)]) )
    } else{
      plot(m[i,j,], type = 'l', ylab = bquote(lambda[paste("(", .(i), ",", .(j), ")")]))
    }
  }
}
dev.off()

# Max posterior sampled tree
png(paste0(new.directory, "/max.post.sample.png"), width = 2000, height = 1500)
structured.plot(max.posterior.sample$ED, n_deme)
color.palette <- rainbow(n_deme)
coal.node.rows <- numeric(n.coal)
coalescence.nodes <- ED[!is.na(ED[,4]), 1]
for (i in 1 : n.coal){
  coal.node.rows[i] <- which(ED[,1] == coalescence.nodes[i])
}
nodelabels(node = coal.node.rows,pie = coal.node.deme.freq/N, piecol = color.palette, cex = 0.75)
dev.off()


# Timescale trace & Hist
png(paste0(new.directory, "/time_scale.plots.png"), width = 2000, height = 1500)
layout(matrix(1:2, 1))
plot(time_scale_sample, type = 'l')
hist(time_scale_sample)
dev.off()

# True migration & coalescent rate hists
png(paste0(new.directory, "/real.params.hists.png"), width = 2000, height = 1500)
layout(matrix(1:n_deme^2, n_deme, n_deme, byrow = TRUE))

for (i in 1:n_deme){
  for (j in 1:n_deme){
    if (i == j){
      hist(mig.eff.pop.sample[i, j, ]*time_scale_sample, freq = FALSE, main = bquote(theta[.(i)]), xlab = bquote(theta[.(i)]))
      abline(v=true_cr[i], col = "red", lty = 2, lwd = 2)

      #Prior convolved with time_scale_sample
      bin.width <- 0.005
      x <- seq(0, 1, bin.width)
      y <- numeric(length(x))
      for (k in 1 : length(x)){
        y[k] <- sum(dgamma(x[k]/time_scale_sample, cr_prior_shape, cr_prior_rate)) / length(time_scale_sample)
      }
      y <- y / (sum(y) * bin.width)
      lines(x,y, col = "blue", lty = 2)

      legend("topright", c("Simulation Parameter", "Prior"), col = c("red", "blue"), lty = 2)

    } else{
      hist(mig.eff.pop.sample[i, j, ]*time_scale_sample, freq = FALSE, main = bquote(lambda[paste("(", .(i), ",", .(j), ")")]), xlab = bquote(lambda[paste("(", .(i), ",", .(j), ")")]))
      abline(v=true_mm[i,j], col = "red", lty = 2, lwd = 2)

      #Prior convolved with time_scale_sample
      bin.width <- 0.005
      x <- seq(0, 1, bin.width)
      y <- numeric(length(x))
      for (k in 1 : length(x)){
        y[k] <- sum(dgamma(x[k]/time_scale_sample, cr_prior_shape, cr_prior_rate)) / length(time_scale_sample)
      }
      y <- y / (sum(y) * bin.width)
      lines(x,y, col = "blue", lty = 2)

      legend("topright", c("Simulation Parameter", "Prior"), col = c("red", "blue"), lty = 2)
    }
  }
}
dev.off()

################ Posterior Estimates on relative params
###### Bootstrap sample posterior density estimate
### Red = posterior using mean time_scale
### Blue = posterior using sampled time_scales

n.dens.samps <- 1000
density.estimate.samples <- sample(1:n.stored.samples, n.dens.samps, TRUE)

x <- seq(0, max(mig.eff.pop.sample), length.out = 101)
x2 <- seq(0, max(time_scale_sample), length.out = 101)

post.dens <- array(0, dim = c(n_deme, n_deme, length(x)))
ts.dens <- numeric(length(x))
mean_ts_post <- array(0, dim = c(n_deme, n_deme, length(x)))
mean_ts <- mean(time_scale_sample)

for (iter in density.estimate.samples){
  ED <- ED_sample[[iter]]
  time_scale <- time_scale_sample[iter]
  mig_coal_rate <- mig.eff.pop.sample[,,iter]


  max_label <- max(ED[,1])
  node_indices <- rep(0, max_label)
  for (l in 1 : dim(ED)[1]){
    node_indices[ED[l,1]] <- l
  }

  deme.decomp <- DemeDecompC(ED, n_deme, node_indices)
  node.count <- NodeCountC(ED, n_deme, node_indices)
  c <- node.count$c
  m <- node.count$m
  k <- deme.decomp$k

  time.increments <- deme.decomp$time.increments

  rate_consts <- colSums((k * (k-1) / 2) * time.increments)
  deme_length <- colSums(k*time.increments)

  coal_rate <- diag(mig_coal_rate)
  mig_mat <- mig_coal_rate
  diag(mig_mat) <- 0

  ts.dens <- ts.dens + dgamma(x2, shape = 0.001 + sum(m) + sum(c), rate = 0.001 + sum(rate_consts * coal_rate) + sum(deme_length * rowSums(mig_mat))) / n.dens.samps

  for (i in 1 : n_deme){
    for (j in 1 : n_deme){
      if (i == j){
        #Coal rate density
        post.dens[i,j,] <- post.dens[i,j,] + dgamma(x, cr_prior_shape + c[i], cr_prior_rate + time_scale * rate_consts[i]) / n.dens.samps

        mean_ts_post[i,j,] <- mean_ts_post[i,j,] + dgamma(x, cr_prior_shape + c[i], cr_prior_rate + mean_ts * rate_consts[i]) / n.dens.samps
      } else{
        #Mig mat density
        post.dens[i,j,] <- post.dens[i,j,] + dgamma(x, mm_prior_shape + m[i,j], mm_prior_rate + deme_length[i] * time_scale) / n.dens.samps

        mean_ts_post[i,j,] <- mean_ts_post[i,j,] + dgamma(x, mm_prior_shape + m[i,j], mm_prior_rate + deme_length[i] * mean_ts) / n.dens.samps
      }
    }
  }



}


png(paste0(new.directory, "/rel.params.post.density.png"), width = 2000, height = 1500)
  layout(matrix(1:n_deme^2, n_deme, n_deme, byrow = TRUE))
  for (i in 1:n_deme){
    for (j in 1:n_deme){
      if (i == j){
        hist(mig.eff.pop.sample[i, j, ], freq = FALSE, main = bquote(theta[.(i)]), xlab = bquote(theta[.(i)]))
        lines(x, post.dens[i,j,], col = "blue", lty = 2, lwd = 2)
        lines(x, mean_ts_post[i,j,], col = "red", lty = 2, lwd = 2)

        legend("topright", c("Mean timescale posterior", "Bootstrap Posterior"), col = c("red", "blue"), lty = 2, lwd = 2)
      } else{
        hist(mig.eff.pop.sample[i, j, ], freq = FALSE, main = bquote(lambda[paste("(", .(i), ",", .(j), ")")]), xlab = bquote(lambda[paste("(", .(i), ",", .(j), ")")]))
        lines(x, post.dens[i,j,], col = "blue", lty = 2, lwd = 2)
        lines(x, mean_ts_post[i,j,], col = "red", lty = 2, lwd = 2)

        legend("topright", c("Mean timescale posterior", "Bootstrap Posterior"), col = c("red", "blue"), lty = 2, lwd = 2)
      }
    }
  }
dev.off()

png(paste0(new.directory, "/time_scale.post.density.png"), width = 2000, height = 1500)
  layout(1)
  hist(time_scale_sample, freq = FALSE)
  lines(x2, ts.dens, col = "blue", lty = 2, lwd = 2)

  legend("topright", "Bootstrap Posterior", col = "blue", lty = 2)
dev.off()
