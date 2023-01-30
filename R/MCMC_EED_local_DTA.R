eed_mcmc_rates <- function(EED, N, n_deme, coal_rate, bit_mig_mat = NA, fit_mig_mat = NA, time_scale,
                           mm_shape = 1, mm_rate = 1, cr_shape = 1, cr_rate = 1, ts_shape = 1e-3, ts_rate = 1e-3,
                           update_rates = TRUE, output_dir = "./MCMC_Results", create.new.directory = TRUE){


  # Check migration matrix is input
  if ((is.null(dim(fit_mig_mat))) & (is.null(dim(bit_mig_mat)))) stop("Input forward-in-time or backward-in-time migration matrix")
  if (is.null(dim(fit_mig_mat))) fit_mig_mat <- FitMigMatC(bit_mig_mat, coal_rate)
  if (is.null(dim(bit_mig_mat))) bit_mig_mat <- BitMigMatC(fit_mig_mat, coal_rate)

  if (create.new.directory){
    output_dir <- file.path(output_dir, format(Sys.time(), "%F_%H_%M"))
    dir.create(output_dir)  #Create directory to store plots; directory name gives date and time
  }

  EED_sample <- list()
  accepted <- 0

  pb <- txtProgressBar(min = 0, max = N, initial = 0, style = 3)

  # First cycle through tree
  coal_nodes <- EED[!is.na(EED[,4]), 1]
  n_coal <- length(coal_nodes)
  jump_times <- numeric(N * n_coal)

  EED_NI <- NodeIndicesC(EED)

  # Initialise storage for rates & time scale estimates
  if (update_rates){
    mm_cr_sample <- array(NA, dim = c(n_deme, n_deme, N))
    ts_sample <- numeric(N)
  }

  current_sc_log_like <- ScaledLikelihoodC(EED, coal_rate, time_scale, bit_mig_mat, EED_NI)$log.likelihood
  current_dta_log_like <- ScaledDTALikelihoodC(EED, coal_rate, time_scale, bit_mig_mat, EED_NI)$log.likelihood

  for (x in 1 : N){
    setTxtProgressBar(pb, x)

    #Cycle through all coalescent nodes of tree
    for (y in 1 : n_coal){
      proposal <- EED_local_DTA(EED, fit_mig_mat, time_scale, coal_nodes[y])
      coal_node_dist <- proposal$node_dist
      proposal <- proposal$proposal
      prop_NI <- NodeIndicesC(proposal)
      prop_sc_log_like <- ScaledLikelihoodC(proposal, coal_rate, time_scale, bit_mig_mat, prop_NI)$log.likelihood
      prop_dta_log_like <- ScaledDTALikelihoodC(proposal, coal_rate, time_scale, bit_mig_mat, prop_NI)$log.likelihood

      # accept_ratio <- min(0, prop_sc_log_like - current_sc_log_like + current_dta_log_like - prop_dta_log_like)
      EED_deme <- EED[EED_NI[coal_nodes[y]], 5]
      prop_deme <- EED[prop_NI[coal_nodes[y]], 5]

      accept_ratio <- min(0, prop_sc_log_like - current_sc_log_like + current_dta_log_like - prop_dta_log_like + log(coal_node_dist[EED_deme]) - log(coal_node_dist[prop_deme]))

      if (log(runif(1)) < accept_ratio){ #ACCEPT
        current_sc_log_like <- prop_sc_log_like
        current_dta_log_like <- prop_dta_log_like
        EED <- proposal
        EED_NI <- prop_NI
        accepted <- accepted + 1
        jump_times[(x-1) * n_coal + y] <- 1
      }

      EED_sample[[(x-1) * n_deme + y]] <- EED
    }

    if (update_rates){
      bit_mig_mat <- scaled_mig_mat_gibbs_update(EED, time_scale, n_deme, EED_NI, prior_shape = mm_shape, prior_rate = mm_rate)
      coal_rate <- scaled_coal_rate_gibbs_update(EED, time_scale, n_deme, EED_NI, prior_shape = cr_shape, prior_rate = cr_rate)
      time_scale <- scaled_time_scale_gibbs_update(coal_rate, bit_mig_mat, EED, n_deme, EED_NI, prior_shape = ts_shape, prior_rate = ts_rate)
      fit_mig_mat <- FitMigMatC(bit_mig_mat, coal_rate)

      current_sc_log_like <- ScaledLikelihoodC(EED, coal_rate, time_scale, bit_mig_mat, EED_NI)$log.likelihood
      current_dta_log_like <- ScaledDTALikelihoodC(EED, coal_rate, time_scale, bit_mig_mat, EED_NI)$log.likelihood

      mm_cr_sample[,,x] <- bit_mig_mat
      diag(mm_cr_sample[,,x]) <- coal_rate
      ts_sample[x] <- time_scale
    }
  }

  close(pb)

  message(paste("Acceptance rate:", round(100 * accepted / (N * n_coal), 2), "%"))

  saveRDS(mm_cr_sample, paste0(output_dir, "/mm_cr_sample.RDS"))
  saveRDS(ts_sample, paste0(output_dir, "/ts_sample.RDS"))
  saveRDS(EED_sample, paste0(output_dir, "/ED_sample.RDS"))


  png(file.path(output_dir, "log_like_trace.png"), width = 2000, height = 1500)
    jump_chain_indices <- c(which(jump_times == 1), N * n_coal)
    sc_log_like <- sapply(jump_chain_indices, function(x){
      prop <- EED_sample[[x]]
      return(ScaledLikelihoodC(prop, coal_rate, time_scale, bit_mig_mat, NodeIndicesC(prop))$log.likelihood)
    })
    plot(jump_chain_indices,sc_log_like, type = 'l', xlab = "Index", ylab = "Log-likelihood", main = "EED MCMC")
  dev.off()

  png(file.path(output_dir, "relative_mm_cr_traces.png"), width = 2000, height = 1500)
  layout(matrix(1:n_deme^2, n_deme, n_deme, byrow = TRUE))
  for (i in 1:n_deme){
    for (j in 1:n_deme){
      if (i == j){
        plot(mm_cr_sample[i, j, ], type = 'l', main = bquote(theta[.(i)]), ylab = bquote(theta[.(i)]))
      } else{
        plot(mm_cr_sample[i, j, ], type = 'l', main = bquote(lambda[paste("(", .(i), ",", .(j), ")")]), ylab = bquote(lambda[paste("(", .(i), ",", .(j), ")")]))
      }
    }
  }
  dev.off()

  png(file.path(output_dir, "relative_mm_cr_hists.png"), width = 2000, height = 1500)
  layout(matrix(1:n_deme^2, n_deme, n_deme, byrow = TRUE))
  for (i in 1:n_deme){
    for (j in 1:n_deme){
      if (i == j){
        hist(mm_cr_sample[i, j, ], freq = FALSE, main = bquote(theta[.(i)]), xlab = bquote(theta[.(i)]))
        x <- seq(0, max(mm_cr_sample[i,j,]), length.out = 1e3)
        lines(x, dgamma(x, cr_shape, cr_rate), lty = 2, col = "red")
        legend("topright", c("Prior"), col = c("red"), lty = 2)
      } else{
        hist(mm_cr_sample[i, j, ], freq = FALSE, main = bquote(lambda[paste("(", .(i), ",", .(j), ")")]), xlab = bquote(lambda[paste("(", .(i), ",", .(j), ")")]))
        x <- seq(0, max(mm_cr_sample[i,j,]), length.out = 1e3)
        lines(x, dgamma(x, mm_shape, mm_rate), lty = 2, col = "red")
        legend("topright", c("Prior"), col = c("red"), lty = 2)
      }
    }
  }
  dev.off()

  png(file.path(output_dir, "ts_trace.png"), width = 2000, height = 1500)
  plot(ts_sample, type = 'l', main = "Time scale", ylab = bquote(gamma))
  x <- seq(0, max(ts_sample), length.out = 1e3)
  lines(x, dgamma(x, ts_shape, ts_rate))
  dev.off()

  png(file.path(output_dir, "true_mm_cr_traces.png"), width = 2000, height = 1500)
  layout(matrix(1:n_deme^2, n_deme, n_deme, byrow = TRUE))
  for (i in 1:n_deme){
    for (j in 1:n_deme){
      if (i == j){
        plot(mm_cr_sample[i, j, ] * ts_sample, type = 'l', main = bquote(theta[.(i)]), ylab = bquote(theta[.(i)]))
      } else{
        plot(mm_cr_sample[i, j, ] * ts_sample, type = 'l', main = bquote(lambda[paste("(", .(i), ",", .(j), ")")]), ylab = bquote(lambda[paste("(", .(i), ",", .(j), ")")]))
      }
    }
  }
  dev.off()

  png(file.path(output_dir, "true_mm_cr_hists.png"), width = 2000, height = 1500)
  #Could add convolved priors -> look at MCMC_Scaled_structured_coalescent.R
  layout(matrix(1:n_deme^2, n_deme, n_deme, byrow = TRUE))
  for (i in 1:n_deme){
    for (j in 1:n_deme){
      if (i == j){
        hist(mm_cr_sample[i, j, ] * ts_sample, freq = FALSE, main = bquote(theta[.(i)]), xlab = bquote(theta[.(i)]))
      } else{
        hist(mm_cr_sample[i, j, ] * ts_sample, freq = FALSE, main = bquote(lambda[paste("(", .(i), ",", .(j), ")")]), xlab = bquote(lambda[paste("(", .(i), ",", .(j), ")")]))
      }
    }
  }
  dev.off()
}
