storage.dir <- "./MCMC_Results"
n_task <- 15

n.deme <- 3
effective.pop <- 100 * 1:3

ESS <- array(dim = c(n.deme, n.deme, n_task))
quantiles <- array(dim = c(n.deme, n.deme, n_task, 3))
true_mig_mat <- array(dim = c(n.deme, n.deme, n_task))


mig.prior.mean = 1/10; mig.prior.var = 1/20
mig.prior.shape <- mig.prior.mean^2/mig.prior.var; mig.prior.rate <- mig.prior.mean/mig.prior.var
prior_quantiles <- qgamma(c(0.05, 0.5, 0.95), rate = mig.prior.rate, shape = mig.prior.shape)

for (task_id in 1 : n_task){
  mig.eff.pop.sample <- readRDS(paste0(storage.dir, "/Run_", task_id, "/mig.eff.pop.sample.RDS"))
  for (i in 1 : n.deme){
    for (j in 1 : n.deme){
      ESS[i,j, task_id] <- coda::effectiveSize(mig.eff.pop.sample[i, j, ])
      quantiles[i,j,task_id,] <- quantile(mig.eff.pop.sample[i,j,], probs = c(0.05, 0.5, 0.95))
    }
  }
  rm(mig.eff.pop.sample)

  true_mig_mat[,, task_id] <- readRDS(paste0(storage.dir, "/Run_", task_id, "/true.mig.mat.RDS"))
}

layout(matrix(1:9, 3, 3, byrow = TRUE))
for (i in 1 : n.deme){
  for (j in 1 : n.deme){
    if (i == j){
      plot(quantiles[i,j,,2], ylim = c(min(quantiles[i,j,, ]), max(quantiles[i,j,, ])), xlab = "Task_id", ylab = "Estimated eff_pop")
      arrows(1:n_task, quantiles[i,j,,1], 1:n_task, quantiles[i,j,,3], code = 3, angle = 90, length = 0.05)
      abline(h = effective.pop[i], lty = 2)
    } else{
      plot(true_mig_mat[i,j,], quantiles[i,j,,2], ylim = c(0, max(quantiles[i,j,, ])), xlab = "True mig_rate", ylab = "Estimated mig_rate")
      arrows(true_mig_mat[i,j,], quantiles[i,j,,1], true_mig_mat[i,j,], quantiles[i,j,,3], code = 3, angle = 90, length = 0.05)
      abline(0, 1, lty = 2)
      abline(h = prior_quantiles, col = "red", lty = 2)
    }
  }
}

layout(matrix(1:9, 3, 3, byrow = TRUE))
for (i in 1: n.deme){
  for (j in 1 : n.deme){
    colours <- rep("black", n_task)
    colours[ESS[i,j,] < 100] <- "red"
    plot(ESS[i,j,], type= 'h', ylim = c(0, 1.04 * max(ESS[i,j,])), col = colours)
    points(ESS[i,j,], pch = 4, col = colours)
    abline(h = 100, lty = 2)
  }
}
