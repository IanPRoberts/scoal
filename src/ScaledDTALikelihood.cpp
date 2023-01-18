#include <Rcpp.h>
using namespace Rcpp;

//' @title ScaledLikelihoodC
//' @description Computes the likelihood of a structured coalescent genealogy
//' @param ED NumericMatrix Extended data object representing structured phylogeny
//' @param eff_pop NumericVector Effective population of each deme
//' @param gen_length double Generation length of each individual in the global population
//' @param mig_mat NumericMatrix Backwards-in-time migration rates between pairs of demes
//' @param node_indices NumericVector Vector giving row numbers in ED for node labels
//' @returns List containing log-likelihood and likelihood of the structured coalescent genealogy
//'
//' @export
//'
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]

List ScaledDTALikelihoodC(NumericMatrix ED, NumericVector coal_rate, double time_scale, NumericMatrix mig_mat, NumericVector node_indices) {
  double log_like = 0;
  int n_deme = coal_rate.size();

  NumericMatrix forward_mm(n_deme, n_deme);
  NumericVector forward_mm_rowsum(n_deme);

  for (int i = 0; i < n_deme; ++i){
    for (int j = 0; j < n_deme; ++j){
      forward_mm(i,j) = time_scale * mig_mat(j,i) * coal_rate[i] / coal_rate[j];
      forward_mm_rowsum[i] += forward_mm(i,j);
    }
  }

  int parent_row;
  double time_inc;
  int parent_deme;
  int node_deme;

  for (int i = 0; i < ED.nrow(); ++i){
    if (ED(i,5) > 0){
      parent_row = node_indices[ED(i,1) - 1] - 1;
      time_inc = ED(i,5) - ED(parent_row, 5);

      parent_deme = ED(parent_row, 4) - 1;
      node_deme = ED(i, 4) - 1;

      log_like -= forward_mm_rowsum[parent_deme] * time_inc;

      if (parent_deme != node_deme){
        log_like += log(forward_mm(parent_deme, node_deme));
      }
    }
  }

  List out = List::create(_["log.likelihood"] = log_like , _["likelihood"] = exp(log_like));
  return out;
}
