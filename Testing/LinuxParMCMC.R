set.seed(100)
require(parallel)
require(doParallel)
require(ape)
require(magick)
devtools::load_all()

system.time({

#parallel::detectCores()
no_cores <- 8
no_top_props <- 8
no_cont_props <- 8
no_ts_props <- 1

N0 <- 0 #1e3  #Burn in
N <- 1e3 #1e6  #Main MCMC run

proposal.rates <- c(rep(5, 4), 1, 1, 1)

#Tree generation
n <- 100
n_deme <- 3
time_scale <- 1
coal_rate <- rep(0.1, n_deme)
mig_mat <-matrix(0.01, n_deme, n_deme)
diag(mig_mat) <- 0

leaf_data <- matrix(c(1:n,
                      sample(2018:2022, n, TRUE),
                      sample.int(n_deme, n, TRUE)), n)


true_phylo <- scaled_sim(leaf_data, coal_rate, time_scale, n_deme, mig_mat, FALSE)
true_ED <- phylo.to.ed(true_phylo)
true_cr <- coal_rate * time_scale
true_mm <- mig_mat * time_scale


#Prior Parameters
cr_prior_shape <- cr_prior_rate <- 1
mm_prior_shape <- mm_prior_rate <- 1

#Initial conditions
ED <- initial_ED <- initial.ed(strip.history(true_ED), leaf_data[,3])
coal_rate <- rep(1, n_deme)
mig_mat <- matrix(1, n_deme, n_deme); diag(mig_mat) <- 0
time_scale <- 1

node_indices <- NodeIndicesC(ED)
ED_like <- ScaledLikelihoodC(ED, coal_rate, time_scale, mig_mat, node_indices)$log.likelihood

freq <- matrix(0, 2, 10)  #Row 1 no. of accepted proposals, row 2 no. of proposals
top.probs <- cumsum(proposal.rates[1:4]/sum(proposal.rates[1:4]))
cont.probs <- cumsum(proposal.rates[5:6]/sum(proposal.rates[5:6]))
ts.probs <- 1

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

new.directory <- file.path("./MCMC_Results", format(Sys.time(), "%F_%H_%M"))
dir.create(new.directory)  #Create directory to store plots; directory name gives date and time

#Progress bar
pb <- txtProgressBar(min = 0, max = N0 + N, initial = 0, style = 3)

topology_proposals <- function(iter){
  if (U[iter, 1] < top.probs[1]){
    if (U[iter, 2] < 0.5){
      which.move <- 1
      proposal <- ed.mig.birth(ED, n_deme, TRUE, node_indices)
    } else{
      which.move <- 2
      proposal <- ed.mig.death(ED, n_deme, TRUE, node_indices)
    }
  } else if (U[iter, 1] < top.probs[2]){
    if (U[iter, 2] < 0.5){
      which.move <- 3
      proposal <- ed.mig.pair.birth(ED, n_deme, node_indices)
    } else{
      which.move <- 4
      proposal <- ed.mig.pair.death(ED, n_deme, node_indices)
    }
  } else if (U[iter, 1] < top.probs[3]){
    if (U[iter, 2] < 0.5){
      which.move <- 5
      proposal <- ed.coal.split(ED, n_deme, node_indices)
    } else{
      which.move <- 6
      proposal <- ed.coal.merge(ED, n_deme, node_indices)
    }
  } else if (U[iter, 1] < top.probs[4]){
    which.move <- 7
    proposal <- ed.block.recolour(ED, n_deme, TRUE, node_indices)
  }
  if (proposal$prop.ratio > 0){
    proposal$like <- ScaledLikelihoodC(proposal$ED, coal_rate, time_scale, mig_mat, proposal$node.indices)$log.likelihood
    log.accept.prob <- min(0, proposal$like - ED_like + proposal$log.prop.ratio)
  } else{
    proposal$like <- NA
    log.accept.prob <- -Inf
  }
  out <- list(accept = FALSE,
              which.move = which.move,
              proposal = proposal)
  if (log(U[iter, 3]) <= log.accept.prob){
    out$accept <- TRUE
  }
  return(out)
}

continuous_proposals <- function(iter){
  if (iter == 1){
    coal_rate <- scaled_coal_rate_gibbs_update(ED, time_scale, n_deme, node_indices, cr_prior_shape, cr_prior_rate)
    coal_rate_prior <- dgamma(coal_rate, shape = cr_prior_shape, rate = cr_prior_rate, log = TRUE)
    return(list(coal_rate = coal_rate,
                coal_rate_prior = coal_rate_prior))
  } else {
    which.move <- 9
    mig_mat <- scaled_mig_mat_gibbs_update(ED, time_scale, n_deme, node_indices, mm_prior_shape, mm_prior_rate)
    mig_mat_prior <- dgamma(mig_mat, shape = mm_prior_shape, rate = mm_prior_rate, log = TRUE)
    diag(mig_mat_prior) <- 0
    return(list(mig_mat = mig_mat,
                mig_mat_prior = mig_mat_prior))
  }
  # out <- list(which.move = which.move,
  #             coal_rate = coal_rate,
  #             coal_rate_prior = coal_rate_prior,
  #             mig_mat = mig_mat,
  #             mig_mat_prior = mig_mat_prior,
  #             ED_like = ED_like)
}

total_props <- no_top_props + no_cont_props + no_ts_props

for (iteration in -N0 : N){
  no_proposals <- 0
  while (no_proposals < no_top_props){
    U <- matrix(runif((no_top_props - no_proposals) * 3), nrow = no_top_props - no_proposals)
    top_proposals <- mclapply(1:(no_top_props - no_proposals), topology_proposals, mc.cores = no_cores)
    accepted <- sapply(1:(no_top_props - no_proposals), function(x){top_proposals[[x]]$accept})

    if (sum(accepted) > 0){
      first_accepted <- min(which(accepted == TRUE))
      no_proposals <- no_proposals + first_accepted

      freq[1, top_proposals[[first_accepted]]$which.move] <- freq[1, top_proposals[[first_accepted]]$which.move] + 1
      proposal <- top_proposals[[first_accepted]]$proposal
      ED <- proposal$ED
      ED_like <- proposal$like
      node_indices <- proposal$node.indices

    } else {
      first_accepted <- no_top_props
    }
  }

  cont_updates <- mclapply(1:2, continuous_proposals, mc.cores = no_cores)

  coal_rate <- cont_updates[[1]]$coal_rate
  coal_rate_prior <- cont_updates[[1]]$coal_rate_prior
  freq[1, 8] <- freq[1, 8] + 1

  mig_mat <- cont_updates[[2]]$mig_mat
  mig_mat_prior <- cont_updates[[2]]$mig_mat_prior
  freq[1, 9] <- freq[1,9] + 1


  time_scale <- scaled_time_scale_gibbs_update(coal_rate, mig_mat, ED, n_deme, node_indices, 0.001, 0.001)
  ED_like <- ScaledLikelihoodC(ED, coal_rate, time_scale, mig_mat, node_indices)$log.likelihood
  freq[1, 10] <- freq[1, 10] + 1

  #Progress bar
  setTxtProgressBar(pb, iteration)

}

close(pb)
#
# png(paste0(new.directory, "/mig.eff.pop.sample.hist.png"), width = 2000, height = 1500)
# layout(matrix(1:n_deme^2, n_deme, n_deme, byrow = TRUE))
# for (i in 1:n_deme){
#   for (j in 1:n_deme){
#     if (i == j){
#       upper <- ceiling(5*max(mig.eff.pop.sample[i, j, ]))/5
#       hist(mig.eff.pop.sample[i, j, ], freq = FALSE, main = bquote(theta[.(i)]), xlab = bquote(theta[.(i)]), breaks = seq(0, upper, 1/5))
#     } else{
#       upper <- ceiling(2*max(mig.eff.pop.sample[i, j, ]))/2
#       hist(mig.eff.pop.sample[i, j, ], freq = FALSE, main = bquote(lambda[paste("(", .(i), ",", .(j), ")")]), xlab = bquote(lambda[paste("(", .(i), ",", .(j), ")")]), breaks = seq(0, upper, 1/10), xlim = c(0,2))
#     }
#   }
# }
# dev.off()
#
# # Matrix plot for trace plots of migration rates and effective population sizes
# png(paste0(new.directory, "/mig.eff.pop.sample.trace.png"), width = 2000, height = 1500)
# layout(matrix(1:n_deme^2, n_deme, n_deme, byrow = TRUE))
# for (i in 1:n_deme){
#   for (j in 1:n_deme){
#     if (i == j){
#       plot(mig.eff.pop.sample[i,j,], type = 'l', main = bquote(paste(theta[.(i)], ", mean =", .(round(mean(mig.eff.pop.sample[i,j,]), 4)))), ylab = bquote(theta[.(i)]) )
#     } else{
#       plot(mig.eff.pop.sample[i,j,], type = 'l', main = bquote(paste(lambda[paste("(", .(i), ",", .(j), ")")], ", mean =", .(round(mean(mig.eff.pop.sample[i,j,]), 4)))), ylab = bquote(lambda[paste("(", .(i), ",", .(j), ")")]))
#     }
#   }
# }
# dev.off()
#
# # Matrix plot for histogram and trace plot of number of migrations
#
# n.mig.sample <- numeric(n.stored.samples)
# for (i in 1 : n.stored.samples){
#   n.mig.sample[i] <- dim(ED_sample[[i]])[1] - 2*n + 1
# }
# png(paste0(new.directory, "/N.mig.plots.png"), width = 2000, height = 1500)
# layout(matrix(1:2, 1, 2))
# plot(n.mig.sample, type = 'l', main = "Trace Plot", xlab = "N.mig")
# hist(n.mig.sample, freq = FALSE, main = "Observed number of migrations", xlab = "N.mig", breaks = 0:max(n.mig.sample + 1) - 0.5)
# dev.off()
#
#
# #Maximum posterior sampled tree with coalescence node deme pie charts
# png(paste0(new.directory, "/max.post.sample.png"), width = 2000, height = 1500)
# structured.plot(max.posterior.sample$ED, n_deme)
# color.palette <- rainbow(n_deme)
# coal.node.rows <- numeric(n.coal)
# coalescence.nodes <- ED[!is.na(ED[,4]), 1]
# for (i in 1 : n.coal){
#   coal.node.rows[i] <- which(ED[,1] == coalescence.nodes[i])
# }
# nodelabels(node = coal.node.rows,pie = coal.node.deme.freq/N, piecol = color.palette, cex = 0.75)
# dev.off()
#
#
# #TimeScale trace & Hist
# png(paste0(new.directory, "/time_scale.plots.png"), width = 2000, height = 1500)
# layout(matrix(1:2, 1))
# plot(time_scale_sample, type = 'l')
# hist(time_scale_sample)
# dev.off()
#
# #
# png(paste0(new.directory, "/mig.eff.pop.sample.hist.2.png"), width = 2000, height = 1500)
# layout(matrix(1:n_deme^2, n_deme, n_deme, byrow = TRUE))
# for (i in 1:n_deme){
#   for (j in 1:n_deme){
#     if (i == j){
#       upper <- ceiling(5*max(mig.eff.pop.sample[i, j, ]))/5
#       hist(mig.eff.pop.sample[i, j, ]*time_scale_sample, freq = FALSE, main = bquote(theta[.(i)]), xlab = bquote(theta[.(i)]))
#     } else{
#       upper <- ceiling(2*max(mig.eff.pop.sample[i, j, ]))/2
#       hist(mig.eff.pop.sample[i, j, ]*time_scale_sample, freq = FALSE, main = bquote(lambda[paste("(", .(i), ",", .(j), ")")]), xlab = bquote(lambda[paste("(", .(i), ",", .(j), ")")]))
#     }
#   }
# }
# dev.off()
#
#
# png(paste0(new.directory, "/real.params.hists.png"), width = 2000, height = 1500)
# layout(matrix(1:n_deme^2, n_deme, n_deme, byrow = TRUE))
# for (i in 1:n_deme){
#   for (j in 1:n_deme){
#     if (i == j){
#       bin.width <- 1/100
#       upper <- ceiling(max(mig.eff.pop.sample[i,j,] * time_scale_sample /bin.width)) * bin.width
#       hist(mig.eff.pop.sample[i, j, ]*time_scale_sample, freq = FALSE, main = bquote(theta[.(i)]), xlab = bquote(theta[.(i)]), breaks = seq(0, upper, bin.width), xlim = c(0, quantile(mig.eff.pop.sample[i,j,], 0.8)))
#       abline(v=true_cr[i], col = "red", lty = 2, lwd = 2)
#     } else{
#       bin.width <- 1/100
#       upper <- ceiling(max(mig.eff.pop.sample[i,j,] * time_scale_sample /bin.width)) * bin.width
#       hist(mig.eff.pop.sample[i, j, ]*time_scale_sample, freq = FALSE, main = bquote(lambda[paste("(", .(i), ",", .(j), ")")]), xlab = bquote(lambda[paste("(", .(i), ",", .(j), ")")]), breaks = seq(0, upper, bin.width), xlim = c(0, quantile(mig.eff.pop.sample[i,j,], 0.8)))
#       abline(v=true_mm[i,j], col = "red", lty = 2, lwd = 2)
#     }
#   }
# }
# dev.off()

})
