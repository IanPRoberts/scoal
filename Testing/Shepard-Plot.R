big_sample <- readRDS(paste0("./MCMC_results/SameTopology/Run_1/ED_sample.RDS"))
k <- 100
subsample <- c(1, 1:k * length(big_sample)/k)

ED_sample <- list()
for (i in 1 : (k+1)){
  ED_sample[[i]] <- big_sample[[subsample[i]]]
}


n_samples <- length(ED_sample)

dists <- matrix(0, n_samples, n_samples)

for (i in 1 : (n_samples - 1)){
  ED1 <- ED_sample[[i]]
  node_indices_1 <- rep(0, max(ED1[,1]))
  for (j in 1 : dim(ED1)[1]){
    node_indices_1[ED1[j,1]] <- j
  }
  for (j in ((i+1) : n_samples)){
    ED2 <- ED_sample[[j]]
    dists[i,j] <- dists[j,i] <- ED.dist(ED1, ED2, node_indices_1)
  }
}

MDS <- vegan::metaMDS(dists)
MDS_dists <- matrix(0, n_samples, n_samples)
for (i in 1 : (n_samples - 1)){
  for (j in ((i+1) : n_samples)){
    MDS_dists[i,j] <- sqrt(sum((MDS$points[i,] - MDS$points[j,])^2))
  }
}

plot(as.vector(dists), as.vector(MDS_dists))
