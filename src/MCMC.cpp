#include <Rcpp.h>
#include "DemeDecomp.h"
#include "NodeCount.h"
#include "StructuredLikelihood.h"
using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]

NumericVector mcmc_cpp(int N0, int N,
              NumericMatrix ED, NumericVector eff_pop, double gen_len, NumericMatrix mig_mat, int n_deme,
              NumericVector prop_rates,
              double eff_pop_prior_mean, double eff_pop_prior_var, double mig_prior_mean, double mig_prior_var,
              CharacterVector likelihood,
              bool output_plots,
              CharacterVector output_folder){

  //Prior Parameters
  double eff_pop_prior_shape = pow(eff_pop_prior_mean, 2)/eff_pop_prior_var + 2;
  double eff_pop_prior_rate = eff_pop_prior_mean * (pow(eff_pop_prior_mean, 2)/eff_pop_prior_var + 1);
  double mig_prior_shape = pow(mig_prior_mean,2)/mig_prior_var;
  double mig_prior_rate = mig_prior_mean/mig_prior_var;

  NumericVector node_indices(max(ED(_,0)), -1);
  for (int j = 0; j < ED.nrow(); ++j){
    node_indices(ED(j, 0) - 1) = j;
  }

  List ED_like = StructuredLikelihoodC(ED, eff_pop, gen_len, mig_mat, node_indices + 1);

  IntegerMatrix freq(2,9);

  NumericVector prop_probs = cumsum(prop_rates);
  prop_probs = prop_probs/prop_probs[prop_probs.size() - 1];

  // Output Setup

  return prop_probs;
}


//
// # Output setup
//       n_stored_samples <- min(1e4, N)
//         mig_eff_pop_sample <- array(0, c(n_deme, n_deme, n_stored_samples))
//         ED_sample <- list()
//         samples_to_store <- round(seq_int(1, N, length_out = n_stored_samples))
//         sample_count <- 1
//
//       coalescence_nodes <- ED[!is_na(ED[,4]), 1]
//       n_coal <- length(coalescence_nodes)
//         coal_node_deme_freq <- matrix(0, n_coal, n_deme, dimnames = list(coalescence_nodes, NULL)) #Matrix to store freq of deme at each coalescence node during iteration
//
//         mig_mat_prior <- dgamma(mig_mat, shape = mig_prior_shape, rate = mig_prior_rate, log = TRUE)
//         diag(mig_mat_prior) <- 0
//       eff_pop_prior <- dgamma(1/effective_pop, shape = eff_pop_prior_shape, rate = eff_pop_prior_rate, log = TRUE)
//         log_prior <- sum(eff_pop_prior) + sum(mig_mat_prior)
//
//         max_posterior_sample <- list(ED = ED,
//                                      effective_population = effective_pop,
//                                      mig_mat = mig_mat,
//                                      iteration = -N0,
//                                      log_likelihood = ED_like,
//                                      log_posterior = ED_like + log_prior)
//
// #Create folders to store images
//         start_time <- Sys_time()
//           new_directory <- file_path(output_folder, format(start_time, "%F_%H-%M"))
//           dir_create(new_directory)  #Create directory to store plots; directory name gives date and time
//
//             if (output_plots == TRUE){
//               dir_create(paste0(new_directory,"/GifOut"))
//               k <- 100  #Number of trees to store throughout MCMC run
//               png(paste0(new_directory,"/GifOut/Frame_", sprintf("%d", 0),"_png"))
//               structured_plot(ed_to_phylo(ED), n_deme)
//               dev_off()
//               video_count <- 1
//             }
//
// # Progress bar
//             pb <- txtProgressBar(min = 0, max = N0 + N, initial = 0, style = 3)
//
//               for (i in -N0 : N){
//                 U <- runif(1)
//                 V <- runif(1)
//                 W <- runif(1)
//
//                 if (U < proposal_probs[1]){
//                   if (V < 0_5){
//                     which_move <- 1
//                     proposal <- ed_mig_birth(ED, n_deme, TRUE, node_indices)
//                   } else{
//                     which_move <- 2
//                     proposal <- ed_mig_death(ED, n_deme, TRUE, node_indices)
//                   }
//                 } else if (U < proposal_probs[2]){
//                   if (V < 0_5){
//                     which_move <- 3
//                     proposal <- ed_mig_pair_birth(ED, n_deme, node_indices)
//                   } else{
//                     which_move <- 4
//                     proposal <- ed_mig_pair_death(ED, n_deme, node_indices)
//                   }
//                 } else if (U < proposal_probs[3]){
//                   if (V < 0_5){
//                     which_move <- 5
//                     proposal <- ed_coal_split(ED, n_deme, node_indices)
//                   } else{
//                     which_move <- 6
//                     proposal <- ed_coal_merge(ED, n_deme, node_indices)
//                   }
//                 } else if (U < proposal_probs[4]){
//                   which_move <- 7
//                   proposal <- ed_block_recolour(ED, n_deme, TRUE, node_indices)
//                 } else if (U < proposal_probs[5]){
//                   which_move <- 8
//                   effective_pop <- eff_pop_update(ED, effective_pop, n_deme, node_indices, shape = eff_pop_prior_shape, rate = eff_pop_prior_rate)
//                   eff_pop_prior <- dgamma(1/effective_pop, shape = eff_pop_prior_shape, rate = eff_pop_prior_rate, log = TRUE)
//                   ED_like <- likelihood_func(ED, effective_pop, gen_len, mig_mat, node_indices)$log_likelihood
//                 } else if (U < proposal_probs[6]){
//                   which_move <- 9
//                   mig_mat <- mig_rate_update(ED, mig_mat, n_deme, node_indices, shape = mig_prior_shape, rate = mig_prior_rate)
//                   mig_mat_prior <- dgamma(mig_mat, shape = mig_prior_shape, rate = mig_prior_rate, log = TRUE)
//                   diag(mig_mat_prior) <- 0
//                   ED_like <- likelihood_func(ED, effective_pop, gen_len, mig_mat, node_indices)$log_likelihood
//                 }
//
//                 freq[2, which_move] <- freq[2, which_move] + 1
//
//                 if ((which_move <= 7) && (proposal$prop_ratio > 0)){
//                   proposal_like <- likelihood_func(proposal$ED, effective_pop, gen_len, mig_mat, proposal$node_indices)$log_likelihood
//                   log_accept_prob <- min(0, proposal_like - ED_like + proposal$log_prop_ratio)
//                   if (log(W) <= log_accept_prob){
//                     freq[1, which_move] <- freq[1, which_move] + 1
//                     ED <- proposal$ED
//                     ED_like <- proposal_like
//                     node_indices <- proposal$node_indices
//                   }
//                 }
//
//                 log_posterior <- sum(eff_pop_prior) + sum(mig_mat_prior) + ED_like
//                   if (log_posterior >= max_posterior_sample$log_posterior){
//                     max_posterior_sample <- list(ED = ED, effective_population = effective_pop, mig_mat = mig_mat, iteration = i, log_likelihood = ED_like, log_posterior = log_posterior)
//                   }
//
//                   if (i > 0){
// #Record deme at coalescent nodes for sampled tree
//                     coalescence_nodes <- ED[!is_na(ED[,4]), 1]
//                     for (j in 1:n_coal){
//                       coal_row <- which(ED[,1] == coalescence_nodes[j])
//                       coal_node_deme_freq[j, ED[coal_row, 5]] <- coal_node_deme_freq[j, ED[coal_row, 5]] + 1
//                     }
//                   }
//
//                   if (i %in% samples_to_store){
//                     ED_sample[[sample_count]] <- ED
//                     mig_eff_pop_sample[,,sample_count] <- mig_mat
//                     diag(mig_eff_pop_sample[,,sample_count]) <- effective_pop
//                     sample_count <- sample_count + 1
//                   }
//
//                   if ((output_plots == TRUE) && (i %in% floor((1:(k+1)) * N/(k+1)))){ #Frames for gif output
//                     png(paste0(new_directory,"/GifOut/Frame_", sprintf("%d", video_count), "_png"))
//                     structured_plot(ed_to_phylo(ED), n_deme)
//                     dev_off()
//                     video_count <- video_count + 1
//                   }
//
// #Progress bar
//                   setTxtProgressBar(pb, i + N0)
//               }
//
//               end_time <- Sys_time()
//                 close(pb)
//
//                 file_create(paste0(new_directory, "/details_txt"))
//                 writeLines(c(paste("N0 =", N0),
//                              paste("N =", N),
//                              paste("likelihood =", likelihood),
//                              paste("Start time", format(start_time, "%H-%M")),
//                              paste("End time", format(end_time, "%H-%M")),
//                              paste("Time elapsed =", round(difftime(end_time, start_time, units = "mins"), digits = 2), "mins")),
//                              con = paste0(new_directory, "/details_txt"))
//
//                 file_create(paste0(new_directory, "/freq_txt"))
//                 write_table(freq, file=paste0(new_directory, "/freq_txt"), row_names=FALSE, col_names=TRUE, sep = "\t")
//                 saveRDS(mig_eff_pop_sample, paste0(new_directory, "/mig_eff_pop_sample_RDS"))
//                 saveRDS(ED_sample, paste0(new_directory, "/ED_sample_RDS"))
//                 saveRDS(max_posterior_sample, paste0(new_directory, "/max_posterior_sample_RDS"))
//
//                 file_create(paste0(new_directory, "/coal_node_deme_freq_txt"))
//                 write_table(coal_node_deme_freq, file=paste0(new_directory, "/coal_node_deme_freq_txt"), row_names=TRUE, col_names=FALSE, sep = "\t")
//
//
//                 if (output_plots == TRUE){
// # MCMC gif
//                   img_list <- sapply(paste0(new_directory,"/GifOut/Frame_", sprintf("%d", 0:(k+1)), "_png"), image_read)
//                   img_joined <- image_join(img_list)
//                   img_animated <- image_animate(img_joined, fps = 2)
//                   image_write(image = img_animated, path = paste0(new_directory,"/Vid_gif"))
//
// # Mig rates & eff pop histograms
//                   png(paste0(new_directory, "/mig_eff_pop_sample_hist_png"), width = 2000, height = 1500)
//                     layout(matrix(1:n_deme^2, n_deme, n_deme, byrow = TRUE))
//                     for (i in 1:n_deme){
//                       for (j in 1:n_deme){
//                         if (i == j){
//                           upper <- ceiling(5*max(mig_eff_pop_sample[i, j, ]))/5
//                           hist(mig_eff_pop_sample[i, j, ], freq = FALSE,
//                                main = bquote(theta[_(i)]),
//                                xlab = bquote(theta[_(i)]),
//                                breaks = seq(0, upper, 1/5))
//                         } else{
//                           upper <- ceiling(2*max(mig_eff_pop_sample[i, j, ]))/2
//                           hist(mig_eff_pop_sample[i, j, ], freq = FALSE,
//                                main = bquote(lambda[paste("(", _(i), ",", _(j), ")")]),
//                                xlab = bquote(lambda[paste("(", _(i), ",", _(j), ")")]),
//                                breaks = seq(0, upper, 1/10), xlim = c(0,2))
//                         }
//                       }
//                     }
//                     dev_off()
//
// # Mig rates & eff pop trace plots
//                       png(paste0(new_directory, "/mig_eff_pop_sample_trace_png"), width = 2000, height = 1500)
//                         layout(matrix(1:n_deme^2, n_deme, n_deme, byrow = TRUE))
//                         for (i in 1:n_deme){
//                           for (j in 1:n_deme){
//                             if (i == j){
//                               plot(mig_eff_pop_sample[i,j,], type = 'l',
//                                    main = bquote(paste(theta[_(i)], ", mean =", _(round(mean(mig_eff_pop_sample[i,j,]), 4)))),
//                                    ylab = bquote(theta[_(i)]) )
//                             } else{
//                               plot(mig_eff_pop_sample[i,j,], type = 'l',
//                                    main = bquote(paste(lambda[paste("(", _(i), ",", _(j), ")")], ", mean =", _(round(mean(mig_eff_pop_sample[i,j,]), 4)))),
//                                    ylab = bquote(lambda[paste("(", _(i), ",", _(j), ")")]))
//                             }
//                           }
//                         }
//                         dev_off()
//
// # No of migrations histogram & trace plot
//                           n_mig_sample <- numeric(n_stored_samples)
//                             for (i in 1 : n_stored_samples){
//                               n_mig_sample[i] <- dim(ED_sample[[i]])[1] - 2*n + 1
//                             }
//                             png(paste0(new_directory, "/N_mig_plots_png"), width = 2000, height = 1500)
//                               layout(matrix(1:2, 1, 2))
//                               plot(n_mig_sample, type = 'l',
//                                    main = "Trace Plot",
//                                    xlab = "N_mig")
//                               hist(n_mig_sample, freq = FALSE,
//                                    main = "Observed number of migrations",
//                                    xlab = "N_mig",
//                                    breaks = 0:max(n_mig_sample + 1) - 0_5)
//                               dev_off()
//
//
// #Maximum posterior sampled tree with coalescence node deme pie charts
//                               png(paste0(new_directory, "/max_post_sample_png"), width = 2000, height = 1500)
//                                 structured_plot(max_posterior_sample$ED, n_deme)
//                                 color_palette <- rainbow(n_deme)
//                                 coal_node_rows <- numeric(n_coal)
//                                 coalescence_nodes <- ED[!is_na(ED[,4]), 1]
//                               for (i in 1 : n_coal){
//                                 coal_node_rows[i] <- which(ED[,1] == coalescence_nodes[i])
//                               }
//                               nodelabels(node = coal_node_rows,pie = coal_node_deme_freq/N, piecol = color_palette, cex = 0_75)
//                                 dev_off()
//                 }
// }
