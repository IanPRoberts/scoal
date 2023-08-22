#' DTA Sampling
#'
#' Samples N migration histories under the DTA model deme fixed at all leaves using belief propagation
#'
#'
#' @param N Number of MCMC iterations (including burn-in)
#' @param ED Extended data representation of a phylogeny including migration history
#' @param coal_rate Vector of coalescent rates in each deme
#' @param bit_mig_mat Matrix of backwards-in-time migration rates
#' @param fit_mig_mat (optional) Matrix of forwards-in-time migration rates
#' @param cr_mode Mode of prior on coalescent rates
#' @param cr_var Variance of prior on coalescent rates
#' @param mm_mode Mode of prior on backwards-in-time migration rates
#' @param mm_var Variance of prior on backwards-in-time migration rates
#' @param output_dir Directory to output log, tree and freq files to
#' @param run_name Prefix for log, tree and freq files. Outputs as 'run_name'.log etc
#' @param thin Thinning rate for MCMC samples
#' @param prop_rates Relative frequencies for each proposal type in the ratio Local_DTA:coal_rate update:bit_mig_mat update
#'
#' @export

local_DTA_mcmc <- function(N = 1e6, ED, coal_rate, bit_mig_mat, fit_mig_mat = FitMigMatC(bit_mig_mat, coal_rate),
                           st_radius = 5,
                           cr_mode = 1, cr_var = 1,
                           mm_mode = 0.05, mm_var = 0.5,
                           output_dir = '~', run_name = 'Local_DTA',
                           thin = 1e3,
                           prop_rates = c(10, 1, 1)){
  pb <- txtProgressBar(0, N, style = 3)

  # Prior parameters
  cr_rate <- (cr_mode + sqrt(cr_mode^2 + 4 * cr_var))/(2 * cr_var)
  cr_shape <- 1 + cr_mode * cr_rate
  mm_rate <- (mm_mode + sqrt(mm_mode^2 + 4 * mm_var))/(2 * mm_var)
  mm_shape <- 1 + mm_mode * mm_rate

  # Forward-in-time rates matrix
  fit_rates <- fit_mig_mat
  diag(fit_rates) <- - rowSums(fit_mig_mat)

  eigen_decomp <- eigen(fit_rates)
  inverse_vecs <- solve(eigen_decomp$vectors)

  # Initial likelihoods and priors
  ED_NI <- NodeIndicesC(ED)
  ED_SC <- SC_like_C(ED, coal_rate, bit_mig_mat, ED_NI)

  mm_prior <- sum(dgamma(bit_mig_mat, mm_shape, mm_rate, log = TRUE)[-(1 + 0:(n_deme - 1) * (n_deme + 1))])
  cr_prior <- sum(dgamma(coal_rate, cr_shape, cr_rate, log = TRUE))

  # Initialise output files
  freq_file <- file.path(output_dir,
                        paste0(run_name, ".freq"))
  freq <- matrix(0, 2, 3,
                 dimnames = list(c("#accept", "#reject"), c("DTA", "CR", "MM")))

  # log_file stores mig rates and coal rates etc
  log_file <- file.path(output_dir,
                        paste0(run_name, ".log"))
  # Column names
  cat("sample",
      "likelihood",
      "posterior",
      paste("coal_rate_", 1:n_deme, sep = "", collapse =","),
      paste("backward_migration_rate_", as.vector(outer(1:n_deme, 1:n_deme, paste, sep = "_")[-(1 + 0:(n_deme - 1) * (n_deme + 1))]), sep = ""),
      file = log_file, sep =",")

  # First row
  cat(paste0("\n", 0), #sample
      ED_SC, # likelihood
      ED_SC + mm_prior + cr_prior, # posterior
      coal_rate, #coal_rate
      as.vector(bit_mig_mat)[-(1 + 0:(n_deme - 1) * (n_deme + 1))], #Mig mat
      file = log_file, sep =",", append = TRUE)

  # tree_file stores posterior sampled trees
  tree_file <- file.path(output_dir,
                         paste0(run_name, ".trees"))
  phylo <- ed.to.phylo(ED)
  treedata <- treeio::as.treedata(phylo)
  treedata@data <- tidytree::tibble(type = paste0("\"", phylo$node.deme, "\""), node = 1:length(phylo$node.deme))

  # .nexus file header
  header <- capture.output(write.nexus(phylo, file = stdout(), translate = TRUE))
  header[2] <- paste("[R-package scoal, ", date(), "]\n\n", sep = "")
  cat(header[-(length(header)-0:1)], file = tree_file, sep = "\n")

  #Initial tree
  cat("\tTREE STATE_0 = ",
      treeio::write.beast.newick(treedata), "\n",
      file = tree_file, append = TRUE, sep = "")

  for (x in 1 : N){
    move_id <- sample(1:3, 1, prob = prop_rates)

    if (move_id == 1){
      subtree <- st_centre_dist(ED, st_radius, ED_NI)
      proposal <- local_DTA_subtree_proposal(subtree$EED, subtree$st_labels, fit_rates, eigen_decomp = eigen_decomp, inverse_vecs = inverse_vecs)

      # Early acceptance if prop == EED
      if ((nrow(ED) == nrow(proposal$proposal)) && (all(na.omit(as.vector(ED == proposal$proposal))))){
        freq[1, move_id] <- freq[1, move_id] + 1
      } else {
        prop_NI <- NodeIndicesC(proposal$proposal)
        prop_SC <- SC_like_C(proposal$proposal, coal_rate, bit_mig_mat, prop_NI)
        log_AR <- min(0, prop_SC - ED_SC +
                        local_DTA_likelihood(subtree$st_labels, coal_rate, bit_mig_mat)$log.likelihood -
                        proposal$prop_prob)

        if (log(runif(1)) < log_AR){
          # Accept
          freq[1, move_id] <- freq[1, move_id] + 1

          ED <- proposal$proposal
          ED_NI <- prop_NI
          ED_SC <- prop_SC
        }
      }
    } else {
      if (move_id == 2){
        coal_rate <- scaled_coal_rate_gibbs_update(ED, 1, n_deme, ED_NI, cr_shape, cr_rate)
        cr_prior <- sum(dgamma(coal_rate, cr_shape, cr_rate, log = TRUE))
      } else if (move_id == 3){
        bit_mig_mat <- scaled_mig_mat_gibbs_update(ED, 1, n_deme, ED_NI, mm_shape, mm_rate)
        mm_prior <- sum(dgamma(bit_mig_mat, mm_shape, mm_rate, log = TRUE)[-(1 + 0:(n_deme - 1) * (n_deme + 1))])
      }

      ED_SC <- SC_like_C(ED, coal_rate, bit_mig_mat, ED_NI)
      freq[1, move_id] <- freq[1, move_id] + 1

      fit_mig_mat <- FitMigMatC(bit_mig_mat, coal_rate)
      fit_rates <- fit_mig_mat
      diag(fit_rates) <- - rowSums(fit_mig_mat)

      eigen_decomp <- eigen(fit_rates)
      inverse_vecs <- solve(eigen_decomp$vectors)
    }

    freq[2, move_id] <- freq[2, move_id] + 1


    if (x %% thin == 0){ #x (mod thin) = 0, i.e. thin-many iterations have passed
      setTxtProgressBar(pb, x)
      cat(paste0("\n", x), #sample
          ED_SC, # likelihood
          ED_SC + mm_prior + cr_prior, # posterior
          coal_rate, #coal_rate
          as.vector(bit_mig_mat)[-(1 + 0:(n_deme - 1) * (n_deme + 1))], #Mig mat
          file = log_file, sep =",", append = TRUE)

      phylo <- ed.to.phylo(ED)
      treedata <- treeio::as.treedata(phylo)
      treedata@data <- tidytree::tibble(type = paste0("\"", phylo$node.deme, "\""), node = 1:length(phylo$node.deme))
      cat("\tTREE STATE_", x, " = ",
          treeio::write.beast.newick(treedata), "\n",
          file = tree_file, append = TRUE, sep = "")

      write.table(freq,
                  file = freq_file,
                  row.names = c('#ACCEPT', '#TOTAL'),
                  col.names = TRUE)
    }

  }
  close(pb)
  cat("END;",
      file = tree_file, append = TRUE, sep = "")
}
