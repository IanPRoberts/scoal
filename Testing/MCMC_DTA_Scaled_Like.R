devtools::load_all()
require(ape)

output_dir <- paste0("/home/ian/R/scoal/MCMC_Results/",format(Sys.time(), "%F_%H_%M"))
dir.create(output_dir)

N0 <- 1e3
N <- 1e5

ED <- readRDS("data/MERS_ED.RDS")
ED <- initial.ed(strip.history(ED, NodeIndicesC(ED)), ED[is.na(ED[,3]), 5])
node_indices <- NodeIndicesC(ED)

n_deme <- max(ED[,5])

mig_mat <- matrix(1, n_deme, n_deme); diag(mig_mat) <- 0
coal_rate <- rep(1, n_deme)
time_scale <- 1

log_like <- ScaledDTALikelihoodC(ED, coal_rate, time_scale, mig_mat, NodeIndicesC(ED))$log.likelihood

proposal_rates <- c(rep(5, 4), 1, 1, 1) # Mig b/d : Mig pair b/d : Coal node s/m : Block recolour : Mig_mat : Coal_rate : Time_scale
cr_prior_shape <- 1
cr_prior_rate <- 1
mm_prior_shape <- 1
mm_prior_rate <- 1

freq <- matrix(0, 2, 10, dimnames = list(c("Accepted", "Proposed"), c("Mig birth", "Mig death", "Pair birth", "Pair death", "Coal split", "Coal merge", "Block recolour", "Mig mat", "Coal rate", "Time scale")))

n_samples <- min(1e5, N)
sample_indices <- round(seq.int(0,N, length.out = n_samples))
sample_count <- 1

ED_sample <- list()
mm_cr_sample <- array(NA, dim = c(n_deme, n_deme, n_samples))
ts_sample <- numeric(n_samples)

pb <- txtProgressBar(min = -N0, max = N, initial = 0, style = 3)
prop_rates <- c(rep(proposal_rates[1], 2), rep(proposal_rates[2], 2), rep(proposal_rates[3], 2), proposal_rates[4:7])

for (iter in (-N0):N){
  prop_choice <- sample(1:10, 1, prob = prop_rates)
  freq[2, prop_choice] <- freq[2, prop_choice] + 1

  if (prop_choice == 1){
    proposal <- ed.mig.birth(ED, n_deme, TRUE, node_indices)
  } else if(prop_choice == 2){
    proposal <- ed.mig.death(ED, n_deme, TRUE, node_indices)
  } else if(prop_choice == 3){
    proposal <- ed.mig.pair.birth(ED, n_deme, node_indices)
  } else if(prop_choice == 4){
    proposal <- ed.mig.pair.death(ED, n_deme, node_indices)
  } else if(prop_choice == 5){
    proposal <- ed.coal.split(ED, n_deme, node_indices)
  } else if(prop_choice == 6){
    proposal <- ed.coal.merge(ED, n_deme, node_indices)
  } else if(prop_choice == 7){
    proposal <- ed.block.recolour(ED, n_deme, TRUE, node_indices)
  } else if(prop_choice == 8){ #Mig_mat update
    for (i in 1 : n_deme){
      for (j in (1 : n_deme)[-i]){
        prop <- mig_mat
        prop[i,j] <- abs(rnorm(1, mig_mat[i,j], 10)) #rnorm.reflect(1, mig_mat[i,j], 100, 0, Inf)
        prop_log_like <- ScaledDTALikelihoodC(ED, coal_rate, time_scale, prop, node_indices)$log.likelihood
        log_ar <- min(0, prop_log_like - log_like)

        if (log(runif(1)) < log_ar){
          mig_mat <- prop
          log_like <- prop_log_like
          freq[1, prop_choice] <- freq[1, prop_choice] + 1
        }
      }
    }
  } else if(prop_choice == 9){ #Coal_rate update
    for (i in 2 : n_deme){
      prop <- coal_rate
      prop[i] <- abs(rnorm(1, coal_rate[i], 100)) #rnorm.reflect(1, coal_rate[i], 1000, 0, Inf)
      prop_log_like <- ScaledDTALikelihoodC(ED, prop, time_scale, mig_mat, node_indices)$log.likelihood
      log_ar <- min(0, prop_log_like - log_like)

      if (log(runif(1)) < log_ar){
        coal_rate <- prop
        log_like <- prop_log_like
        freq[1, prop_choice] <- freq[1, prop_choice] + 1
      }
    }
  } else if(prop_choice == 10){ #Time_scale update
    prop <- abs(rnorm(1, time_scale, 1)) #rnorm.reflect(1, time_scale, 0.1, 0, Inf)
    prop_log_like <- ScaledDTALikelihoodC(ED, coal_rate, prop, mig_mat, node_indices)$log.likelihood
    log_ar <- min(0, prop_log_like - log_like)
    if (log(runif(1)) < log_ar){
      time_scale <- prop
      log_like <- prop_log_like
      freq[1, prop_choice] <- freq[1, prop_choice] + 1
    }
  }

  if ((prop_choice <= 7) && (proposal$prop.ratio > 0)){
    prop_log_like <- ScaledDTALikelihoodC(proposal$ED, coal_rate, time_scale, mig_mat, proposal$node.indices)$log.likelihood
    log_ar <- min(0, prop_log_like - log_like + proposal$log.prop.ratio)
    if (log(runif(1)) < log_ar){
      freq[1, prop_choice] <- freq[1, prop_choice] + 1
      ED <- proposal$ED
      log_like <- prop_log_like
      node_indices <- proposal$node.indices
    }
  }

  if (iter == sample_indices[sample_count]){
    ED_sample[[sample_count]] <- ED
    mm_cr_sample[,,sample_count] <- mig_mat
    diag(mm_cr_sample[,,sample_count]) <- coal_rate
    ts_sample[sample_count] <- time_scale
    sample_count <- sample_count + 1
  }

  setTxtProgressBar(pb, iter)
}
close(pb)

file.create(paste0(output_dir, "/freq.txt"))
write.table(freq, file=paste0(output_dir, "/freq.txt"), row.names=FALSE, col.names=TRUE, sep = "\t")
saveRDS(mm_cr_sample, paste0(output_dir, "/mm_cr_sample.RDS"))
saveRDS(ts_sample, paste0(output_dir, "/ts_sample.RDS"))
saveRDS(ED_sample, paste0(output_dir, "/ED_sample.RDS"))
#saveRDS(max.posterior.sample, paste0(output_dir, "/max_post_sample.RDS"))

# file.create(paste0(output_dir, "/coal_deme_freq.txt"))
# write.table(coal.node.deme.freq, file=paste0(output_dir, "/coal_deme_freq.txt"), row.names=TRUE, col.names=FALSE, sep = "\t")


# Parameter traces (relative)
png(paste0(output_dir, "/traces_relative.png"), width = 2000, height = 1500)
layout(matrix(1:n_deme^2, n_deme, n_deme, byrow = TRUE))
for (i in 1:n_deme){
  for (j in 1:n_deme){
    if (i == j){
      plot(mm_cr_sample[i,j,], type = 'l', main = bquote(paste(s[.(i)], ", mean =", .(round(mean(mm_cr_sample[i,j,]), 4)))), ylab = bquote(s[.(i)]) )
    } else{
      plot(mm_cr_sample[i,j,], type = 'l', main = bquote(paste(r[paste("(", .(i), ",", .(j), ")")], ", mean =", .(round(mean(mm_cr_sample[i,j,]), 4)))), ylab = bquote(r[paste("(", .(i), ",", .(j), ")")]))
    }
  }
}
dev.off()

# Parameter traces (product)
png(paste0(output_dir, "/traces_product.png"), width = 2000, height = 1500)
layout(matrix(1:n_deme^2, n_deme, n_deme, byrow = TRUE))
for (i in 1:n_deme){
  for (j in 1:n_deme){
    if (i == j){
      plot(mm_cr_sample[i,j,] * ts_sample, type = 'l', main = bquote(paste(theta[.(i)], ", mean =", .(round(mean(mm_cr_sample[i,j,]), 4)))), ylab = bquote(theta[.(i)]) )
    } else{
      plot(mm_cr_sample[i,j,] * ts_sample, type = 'l', main = bquote(paste(lambda[paste("(", .(i), ",", .(j), ")")], ", mean =", .(round(mean(mm_cr_sample[i,j,]), 4)))), ylab = bquote(lambda[paste("(", .(i), ",", .(j), ")")]))
    }
  }
}
dev.off()

# Timescale trace & Hist
png(paste0(output_dir, "/time_scale.plots.png"), width = 2000, height = 1500)
layout(matrix(1:2, 1))
plot(ts_sample, type = 'l')
hist(ts_sample)
dev.off()

# True migration & coalescent rate hists
png(paste0(output_dir, "/hists_product.png"), width = 2000, height = 1500)
layout(matrix(1:n_deme^2, n_deme, n_deme, byrow = TRUE))
for (i in 1 : n_deme){
  for (j in 1 : n_deme){
    histogram <- hist(mm_cr_sample[i,j,] * ts_sample, plot = FALSE)

    bin.width <- 0.005
    x <- seq(min(histogram$breaks), max(histogram$breaks), bin.width)
    y <- numeric(length(x))
    if (i==j){
      #y <- sapply(x, function(a) sum(dgamma(a/ts_sample, cr_prior_shape, cr_prior_rate)) / length(ts_sample))
      y <- sapply(x, function(a) sum(dgamma(a/ts_sample, cr_prior_shape, cr_prior_rate)) / length(ts_sample))
      y <- y / (sum(y) * bin.width)
      plot(histogram, freq = FALSE, main = bquote(theta[.(i)]), xlab = bquote(theta[.(i)]))
      lines(x, y, col = "blue", lty = 2)
    } else {
      y <- sapply(x, function(a) sum(dgamma(a/ts_sample, mm_prior_shape, mm_prior_rate))) / length(ts_sample)
      y <- y / (sum(y) * bin.width)
      plot(histogram, freq = FALSE, main = bquote(lambda[paste("(", .(i), ",", .(j), ")")]), xlab = bquote(lambda[paste("(", .(i), ",", .(j), ")")]))
      lines(x, y, col = "blue", lty = 2)
    }
  }
}
dev.off()
