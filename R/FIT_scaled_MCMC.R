#' Runs an MCMC on Migration Histories
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

fit_scaled_MCMC <- function(N0 = 1e5, N = 1e6,
                        ED, coal_rate, time_scale, mig_mat, n_deme = NA,
                        proposal.rates = c(rep(5, 4), 1, 1, 1),
                        cr_prior_shape = 1, cr_prior_rate = 1, mm_prior_shape = 1, mm_prior_rate = 1,
                        likelihood = "Structured",
                        output.directory = "./MCMC_Results", create.new.directory = TRUE){

  node_indices <- NodeIndicesC(ED)

  ED_like <- ScaledLikelihoodC(ED, coal_rate, time_scale, mig_mat, node_indices)$log.likelihood


  freq <- matrix(0, 2, 10)  #Row 1 no. of accepted proposals, row 2 no. of proposals
  proposal.probs <- cumsum(proposal.rates/sum(proposal.rates)) #Cumulative proposal probabilities for each reversible move (single birth/death : pair birth/death : merge/split : block recolour)

  n.stored.samples <- min(N, 1e5)
  mig.eff.pop.sample <- array(0, c(n_deme, n_deme, n.stored.samples))
  time_scale_sample <- numeric(n.stored.samples)
  ED_sample <- list()
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

  if (create.new.directory){
    output.directory <- file.path(output.directory, format(Sys.time(), "%F_%H_%M"))
    dir.create(output.directory)  #Create directory to store plots; directory name gives date and time
  }



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
      coal_rate <- fit_scaled_cr(ED, coal_rate, mig_mat, time_scale, n_deme, node_indices, prior_shape = cr_prior_shape, prior_rate = cr_prior_rate, proposal_sd = time_scale)
      coal_rate_prior <- dgamma(coal_rate, shape = cr_prior_shape, rate = cr_prior_rate, log = TRUE)
      ED_like <- ScaledLikelihoodC(ED, coal_rate, time_scale, mig_mat, node_indices)$log.likelihood
      freq[1, which.move] <- freq[1, which.move] + 1
    } else if (U < proposal.probs[6]){
      which.move <- 9
      mig_mat <- fit_scaled_mm_gibbs(D, coal_rate, time_scale, n_deme, node_indices, mm_prior_shape, mm_prior_rate)
      mig_mat_prior <- dgamma(mig_mat, shape = mm_prior_shape, rate = mm_prior_rate, log = TRUE)
      diag(mig_mat_prior) <- 0
      ED_like <- ScaledLikelihoodC(ED, coal_rate, time_scale, mig_mat, node_indices)$log.likelihood
      freq[1, which.move] <- freq[1, which.move] + 1
    } else if (U < proposal.probs[7]){
      which.move <- 10
      time_scale <- fit_ts_gibbs(ED, coal_rate, mig_mat, n_deme, node_indices, prior_shape = 0.001, prior_rate = 0.001)
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
      ED_sample[[sample.count]] <- ED

      mig.eff.pop.sample[,,sample.count] <- mig_mat
      diag(mig.eff.pop.sample[,,sample.count]) <- coal_rate
      time_scale_sample[sample.count] <- time_scale
      sample.count <- sample.count + 1
    }

    #Progress bar
    setTxtProgressBar(pb, i + N0)
  }

  close(pb)

  file.create(paste0(output.directory, "/freq.txt"))
  write.table(freq, file=paste0(output.directory, "/freq.txt"), row.names=FALSE, col.names=TRUE, sep = "\t")
  saveRDS(mig.eff.pop.sample, paste0(output.directory, "/mm_cr_sample.RDS"))
  saveRDS(time_scale_sample, paste0(output.directory, "/ts_sample.RDS"))
  saveRDS(ED_sample, paste0(output.directory, "/ED_sample.RDS"))
  saveRDS(max.posterior.sample, paste0(output.directory, "/max_post_sample.RDS"))

  file.create(paste0(output.directory, "/coal_deme_freq.txt"))
  write.table(coal.node.deme.freq, file=paste0(output.directory, "/coal_deme_freq.txt"), row.names=TRUE, col.names=FALSE, sep = "\t")


  # Relative migration & coalescent rate hists
  png(paste0(output.directory, "/rel.params.hist.png"), width = 2000, height = 1500)
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
  png(paste0(output.directory, "/rel.params.traces.png"), width = 2000, height = 1500)
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

  png(paste0(output.directory, "/migrations.per.deme.png"), width = 2000, height = 1500)
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
  png(paste0(output.directory, "/max.post.sample.png"), width = 2000, height = 1500)
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
  png(paste0(output.directory, "/time_scale.plots.png"), width = 2000, height = 1500)
  layout(matrix(1:2, 1))
  plot(time_scale_sample, type = 'l')
  hist(time_scale_sample)
  dev.off()

  # True migration & coalescent rate hists
  png(paste0(output.directory, "/real.params.hists.png"), width = 2000, height = 1500)
  layout(matrix(1:n_deme^2, n_deme, n_deme, byrow = TRUE))

  for (i in 1:n_deme){
    for (j in 1:n_deme){
      if (i == j){
        hist(mig.eff.pop.sample[i, j, ]*time_scale_sample, freq = FALSE, main = bquote(theta[.(i)]), xlab = bquote(theta[.(i)]))
        #abline(v=true_cr[i], col = "red", lty = 2, lwd = 2)

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
        #abline(v=true_mm[i,j], col = "red", lty = 2, lwd = 2)

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

}
