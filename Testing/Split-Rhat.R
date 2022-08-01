#Split-Rhat from Gelman et al. (2013)
N <- 100000
M <- 5
n_deme <- 3

##### Continuous parameters (migration rates & coalescent rates)
cr_mms <- list()
for (run_id in 1 : M){
  cr_mms[[run_id]] <- readRDS(paste0("./MCMC_results/SameTopology/Run_", run_id, "/mm_cr_sample.RDS"))
}

within_means <- array(0, c(n_deme, n_deme, M))
for(run_id in 1 : M){
  for (i in 1 : n_deme){
    for (j in 1 : n_deme){
      within_means[i,j,run_id] <- mean(cr_mms[[run_id]][i,j,])
    }
  }
}

pooled_means <- matrix(0, n_deme, n_deme)
for (i in 1 : n_deme){
  for (j in 1 : n_deme){
    pooled_means[i,j] <- mean(within_means[i,j,])
  }
}

within_vars <- array(0, c(n_deme, n_deme, M))
for(k in 1 : M){
  for (i in 1 : n_deme){
    for (j in 1 : n_deme){
      within_vars[i,j,k] <- sum((cr_mms[[k]][i,j,] - within_means[i,j,k])^2)/(N-1)
    }
  }
}

R_hat_mm_cr <- matrix(0, n_deme, n_deme)

B <- matrix(0,n_deme,n_deme)
for (i in 1 : n_deme){
  for (j in 1 : n_deme){
    B[i,j] <- N / (M-1) * sum((within_means[i,j,] - pooled_means[i,j])^2)
  }
}

W <- matrix(0,n_deme,n_deme)
for (i in 1 : n_deme){
  for (j in 1 : n_deme){
    W[i,j] <- mean(within_vars[i,j,])
  }
}

R_hat_mm_cr <- sqrt(((N-1)/N * W + B/N)/W)
R_hat_mm_cr


#### Number of migrations per deme
m <- array(dim = c(n_deme, n_deme, N, M))

for (run_id in 1 : M){
  ED_sample <- readRDS(paste0("./MCMC_results/SameTopology/Run_", run_id, "/ED_sample.RDS"))
  for (i in 1:N){
    ED <- ED_sample[[i]]
    max_label <- max(ED[,1])
    node_indices <- rep(0, max_label)
    for (j in 1 : dim(ED)[1]){
      node_indices[ED[j,1]] <- j
    }

    NC <- NodeCountC(ED, n_deme, node_indices)
    m[,,i, run_id] <- NC$m
    diag(m[,,i, run_id]) <- NC$c
  }
  rm(ED_sample)
}

within_means <- array(0, c(n_deme, n_deme, M))
for(k in 1 : M){
  for (i in 1 : n_deme){
    for (j in 1 : n_deme){
      within_means[i,j,k] <- mean(m[i,j,,k])
    }
  }
}

pooled_means <- matrix(0, n_deme, n_deme)
for (i in 1 : n_deme){
  for (j in 1 : n_deme){
    pooled_means[i,j] <- mean(within_means[i,j,])
  }
}

within_vars <- array(0, c(n_deme, n_deme, M))
for(k in 1 : M){
  for (i in 1 : n_deme){
    for (j in 1 : n_deme){
      within_vars[i,j,k] <- sum((m[i,j,,k] - within_means[i,j,k])^2)/(N-1)
    }
  }
}

R_hat_mig_number <- matrix(0, n_deme, n_deme)

B <- matrix(0,n_deme,n_deme)
for (i in 1 : n_deme){
  for (j in 1 : n_deme){
    B[i,j] <- N / (M-1) * sum((within_means[i,j,] - pooled_means[i,j])^2)
  }
}

W <- matrix(0,n_deme,n_deme)
for (i in 1 : n_deme){
  for (j in 1 : n_deme){
    W[i,j] <- mean(within_vars[i,j,])
  }
}

R_hat_mig_number <- sqrt(((N-1)/N * W + B/N)/W)

R_hat_mig_number
