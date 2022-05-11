#' Runs and MCMC on Migration Histories
#'
#' Runs an MCMC to infer migration histories under the structured coalescent model
#'
#'
#' @param N0 Number of iterations to use as a burn-in
#' @param leaf.deme vector of demes for each leaf in \code{phylo}
#'
#' @return An object of class \code{phylo} augmented with the deme of each node in the tree.
#'
#' @export

migration.history.mcmc <- function(N0 = 1e5, N = 1e6,
                                   ED, effective.pop, gen.length, migration.matrix, n.deme = NA,
                                   proposal.rates = c(rep(5, 4), 1, 1),
                                   eff.pop.prior.mean = 1, eff.pop.prior.var = 1, mig.prior.mean = 1/10, mig.prior.var = 1/20,
                                   likelihood = "Structured",
                                   output.plots = TRUE,
                                   output.directory = "./MCMC_Results", create.new.directory = TRUE){

  if (is.na(n.deme)){
    n.deme <- length(effective.pop)
  }

  #Prior Parameters
  eff.pop.prior.shape <- eff.pop.prior.mean^2/eff.pop.prior.var + 2
  eff.pop.prior.rate <- eff.pop.prior.mean * (eff.pop.prior.mean^2/eff.pop.prior.var + 1)
  mig.prior.shape <- mig.prior.mean^2/mig.prior.var
  mig.prior.rate <- mig.prior.mean/mig.prior.var

  max.label <- max(ED[,1])
  node.indices <- rep(0, max.label)
  for (j in 1 : dim(ED)[1]){
    node.indices[ED[j,1]] <- j
  }

  if (likelihood == "Structured"){
    likelihood.func <- StructuredLikelihoodC
  } else if (likelihood == "DTA"){
    likelihood.func <- dta.likelihood
  } else if (likelihood == "Synthetic"){
    likelihood.func <- synth.likelihood
  }
  else {
    stop("Input likelihood selected from Structured, DTA or Synthetic")
  }

  ED.like <- likelihood.func(ED, effective.pop, gen.length, migration.matrix, node.indices)$log.likelihood

  freq <- matrix(0, 2, 9,
                 dimnames = list(NULL, c("MB", "MD", "MPB", "MPD", "CNS", "CNM", "BR", "EP", "MM")))
  proposal.probs <- cumsum(proposal.rates/sum(proposal.rates)) #Cumulative proposal probabilities for each reversible move (single birth/death : pair birth/death : merge/split : block recolour)

  # Output setup
  n.stored.samples <- min(1e4, N)
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

  max.posterior.sample <- list(ED = ED,
                               effective.population = effective.pop,
                               migration.matrix = migration.matrix,
                               iteration = -N0,
                               log.likelihood = ED.like,
                               log.posterior = ED.like + log.prior)

  #Create folders to store images
  start.time <- Sys.time()
  if (create.new.directory == TRUE){
    storage.directory <- file.path(output.directory, format(start.time, "%F_%H-%M"))
    dir.create(storage.directory)  #Create directory to store plots; directory name gives date and time
  } else{
    storage.directory <- output.directory
  }

  if (output.plots == TRUE){
    dir.create(paste0(storage.directory,"/GifOut"))
    k <- 100  #Number of trees to store throughout MCMC run
    png(paste0(storage.directory,"/GifOut/Frame_", sprintf("%d", 0),".png"))
      structured.plot(ed.to.phylo(ED), n.deme)
    dev.off()
    video.count <- 1
  }

  # Progress bar
  pb <- txtProgressBar(min = 0, max = N0 + N, initial = 0, style = 3)

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
      ED.like <- likelihood.func(ED, effective.pop, gen.length, migration.matrix, node.indices)$log.likelihood
    } else if (U < proposal.probs[6]){
      which.move <- 9
      migration.matrix <- mig.rate.update(ED, migration.matrix, n.deme, node.indices, shape = mig.prior.shape, rate = mig.prior.rate)
      mig.mat.prior <- dgamma(migration.matrix, shape = mig.prior.shape, rate = mig.prior.rate, log = TRUE)
      diag(mig.mat.prior) <- 0
      ED.like <- likelihood.func(ED, effective.pop, gen.length, migration.matrix, node.indices)$log.likelihood
    }

    freq[2, which.move] <- freq[2, which.move] + 1

    if ((which.move <= 7) && (proposal$prop.ratio > 0)){
      proposal.like <- likelihood.func(proposal$ED, effective.pop, gen.length, migration.matrix, proposal$node.indices)$log.likelihood
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

    if ((output.plots == TRUE) && (i %in% floor((1:(k+1)) * N/(k+1)))){ #Frames for gif output
      png(paste0(storage.directory,"/GifOut/Frame_", sprintf("%d", video.count), ".png"))
        structured.plot(ed.to.phylo(ED), n.deme)
      dev.off()
      video.count <- video.count + 1
    }

    #Progress bar
    setTxtProgressBar(pb, i + N0)
  }

  end.time <- Sys.time()
  close(pb)

  file.create(paste0(storage.directory, "/details.txt"))
  writeLines(c(paste("N0 =", N0),
               paste("N =", N),
               paste("likelihood =", likelihood),
               paste("Start time", format(start.time, "%H-%M")),
               paste("End time", format(end.time, "%H-%M")),
               paste("Time elapsed =", round(difftime(end.time, start.time, units = "mins"), digits = 2), "mins")),
             con = paste0(storage.directory, "/details.txt"))

  file.create(paste0(storage.directory, "/freq.txt"))
  write.table(freq, file=paste0(storage.directory, "/freq.txt"), row.names=FALSE, col.names=TRUE, sep = "\t")
  saveRDS(mig.eff.pop.sample, paste0(storage.directory, "/mig.eff.pop.sample.RDS"))
  saveRDS(ED.sample, paste0(storage.directory, "/ED.sample.RDS"))
  saveRDS(max.posterior.sample, paste0(storage.directory, "/max.posterior.sample.RDS"))

  file.create(paste0(storage.directory, "/coal.node.deme.freq.txt"))
  write.table(coal.node.deme.freq, file=paste0(storage.directory, "/coal.node.deme.freq.txt"), row.names=TRUE, col.names=FALSE, sep = "\t")


  if (output.plots == TRUE){
    # MCMC gif
    img_list <- sapply(paste0(storage.directory,"/GifOut/Frame_", sprintf("%d", 0:(k+1)), ".png"), image_read)
    img_joined <- image_join(img_list)
    img_animated <- image_animate(img_joined, fps = 2)
    image_write(image = img_animated, path = paste0(storage.directory,"/Vid.gif"))

    # Mig rates & eff pop histograms
    png(paste0(storage.directory, "/mig.eff.pop.sample.hist.png"), width = 2000, height = 1500)
    layout(matrix(1:n.deme^2, n.deme, n.deme, byrow = TRUE))
    for (i in 1:n.deme){
      for (j in 1:n.deme){
        if (i == j){
          upper <- ceiling(5*max(mig.eff.pop.sample[i, j, ]))/5
          hist(mig.eff.pop.sample[i, j, ], freq = FALSE,
               main = bquote(theta[.(i)]),
               xlab = bquote(theta[.(i)]),
               breaks = seq(0, upper, 1/5))
        } else{
          upper <- ceiling(2*max(mig.eff.pop.sample[i, j, ]))/2
          hist(mig.eff.pop.sample[i, j, ], freq = FALSE,
               main = bquote(lambda[paste("(", .(i), ",", .(j), ")")]),
               xlab = bquote(lambda[paste("(", .(i), ",", .(j), ")")]),
               breaks = seq(0, upper, 1/10), xlim = c(0,2))
        }
      }
    }
    dev.off()

    # Mig rates & eff pop trace plots
    png(paste0(storage.directory, "/mig.eff.pop.sample.trace.png"), width = 2000, height = 1500)
      layout(matrix(1:n.deme^2, n.deme, n.deme, byrow = TRUE))
      for (i in 1:n.deme){
        for (j in 1:n.deme){
          if (i == j){
            plot(mig.eff.pop.sample[i,j,], type = 'l',
                 main = bquote(paste(theta[.(i)], ", mean =", .(round(mean(mig.eff.pop.sample[i,j,]), 4)))),
                 ylab = bquote(theta[.(i)]) )
          } else{
            plot(mig.eff.pop.sample[i,j,], type = 'l',
                 main = bquote(paste(lambda[paste("(", .(i), ",", .(j), ")")], ", mean =", .(round(mean(mig.eff.pop.sample[i,j,]), 4)))),
                 ylab = bquote(lambda[paste("(", .(i), ",", .(j), ")")]))
          }
        }
      }
    dev.off()

    # No of migrations histogram & trace plot
    n.mig.sample <- numeric(n.stored.samples)
    for (i in 1 : n.stored.samples){
      n.mig.sample[i] <- dim(ED.sample[[i]])[1] - 2*n + 1
    }
    png(paste0(storage.directory, "/N.mig.plots.png"), width = 2000, height = 1500)
      layout(matrix(1:2, 1, 2))
      plot(n.mig.sample, type = 'l',
           main = "Trace Plot",
           xlab = "N.mig")
      hist(n.mig.sample, freq = FALSE,
           main = "Observed number of migrations",
           xlab = "N.mig",
           breaks = 0:max(n.mig.sample + 1) - 0.5)
    dev.off()


    #Maximum posterior sampled tree with coalescence node deme pie charts
    png(paste0(storage.directory, "/max.post.sample.png"), width = 2000, height = 1500)
      structured.plot(max.posterior.sample$ED, n.deme)
      color.palette <- rainbow(n.deme)
      coal.node.rows <- numeric(n.coal)
      coalescence.nodes <- ED[!is.na(ED[,4]), 1]
      for (i in 1 : n.coal){
        coal.node.rows[i] <- which(ED[,1] == coalescence.nodes[i])
      }
      nodelabels(node = coal.node.rows,pie = coal.node.deme.freq/N, piecol = color.palette, cex = 0.75)
    dev.off()
  }
}
