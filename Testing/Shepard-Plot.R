big_sample <- readRDS(paste0("C:/Users/ian_p/Documents/R/ScaledScoal/MCMC_results/SameTopology/Run_1/ED_sample.RDS"))
k <- 200
subsample <- c(1, 1:k * length(big_sample)/k)

ED_sample <- list()
for (i in 1 : (k+1)){
  ED_sample[[i]] <- big_sample[[subsample[i]]]
}


n_samples <- length(ED_sample)
n_deme <- 3

dists <- matrix(0, n_samples, n_samples)
root_deme <- numeric(n_samples)
majority_deme <- numeric(n_samples)
majority_coal <- numeric(n_samples)
majority_coal2 <- numeric(n_samples)


for (i in 1 : (n_samples - 1)){
  ED1 <- ED_sample[[i]]
  node_indices_1 <- NodeIndicesC(ED1)
  for (j in ((i+1) : n_samples)){
    ED2 <- ED_sample[[j]]
    node_indices_2 <- NodeIndicesC(ED2)
    dists[i,j] <- dists[j,i] <- ED_dist_C(ED1, ED2, n_deme, node_indices_1, node_indices_2)
  }

  root_node <- which(is.na(ED1[,2]))
  root_deme[i] <- ED1[root_node, 5]

  deme_decomp <- DemeDecompC(ED1, n_deme, node_indices_1)
  deme_lengths <- colSums(deme_decomp$k * deme_decomp$time.increments)
  majority_deme[i] <- which(deme_lengths == max(deme_lengths))

  node_count <- NodeCountC(ED1, n_deme, node_indices_1)
  majority_coal[i] <- sample.vector(which(node_count$c == max(node_count$c)), 1) #Chosen min in case of multiples...
  majority_coal2[i] <- min(which(node_count$c == max(node_count$c))) #Chosen min in case of multiples...
}

dists

ED1 <- ED_sample[[n_samples]]
root_node <- which(is.na(ED1[,2]))
root_deme[n_samples] <- ED1[root_node, 5]

MDS <- vegan::metaMDS(dists)
MDS_dists <- matrix(0, n_samples, n_samples)
for (i in 1 : (n_samples - 1)){
  for (j in ((i+1) : n_samples)){
    MDS_dists[i,j] <- MDS_dists[j,i] <-  sqrt(sum((MDS$points[i,] - MDS$points[j,])^2))
  }
}

plot(as.vector(dists), as.vector(MDS_dists))


plot(MDS$points, col = rainbow(n_deme)[root_deme], pch = 16, main = "MDS plot (root deme)")
plot(MDS$points, col = rainbow(n_deme)[majority_deme], pch = 16, main = "MDS plot (majority deme)")
plot(MDS$points, col = rainbow(n_deme)[majority_coal], pch = 16, main = "MDS plot (majority coalescent node deme)")

plot(MDS$points, col = rainbow(n_deme)[majority_coal2], pch = 16, main = "MDS plot (majority coalescent node deme)")

which(min(rowMeans(dists)) == rowMeans(dists))
MDS$points
points(MDS$points[190,], col = "black")
abline(v= MDS$points[190,1])
abline(h= MDS$points[190,2])
