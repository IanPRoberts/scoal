// [[Rcpp::interfaces(r, cpp)]]
#include <Rcpp.h>
#include "DemeDecomp.h"
#include "NodeCount.h"
#include "StructuredLikelihood.h"

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
// [[Rcpp::export]]

List ScaledLikelihoodC(NumericMatrix ED, NumericVector coal_rate, double time_scale, NumericMatrix mig_mat, NumericVector node_indices) {
  int n_deme = coal_rate.size();
  List deme_decomp = DemeDecompC(ED, n_deme, node_indices);
  List node_count = NodeCountC(ED, n_deme, node_indices);

  NumericMatrix k =  deme_decomp["k"];
  NumericVector time_increments = deme_decomp["time.increments"];
  NumericVector c = node_count["c"];
  NumericMatrix m = node_count["m"];

  double like = 0;
  int n_time_incs = time_increments.size();

  double coal_term;
  double deme_length;
  double mig_row_sum;

  for (int i = 0; i < n_deme; ++i){
    coal_term = 0;
    deme_length = 0;
    mig_row_sum = 0;

    for (int r = 0; r < n_time_incs; ++r){
      coal_term += k(r, i) * (k(r, i) - 1) * time_increments[r] / 2;
      deme_length += k(r, i) * time_increments[r];
    }

    for (int j = 0; j < n_deme; ++j){
      if (i != j){
        like += m(i, j) * log(mig_mat(i, j));
        mig_row_sum += mig_mat(i, j);
      }
    }
    like += c[i] * log(coal_rate[i]) - (coal_term * coal_rate[i] + deme_length * mig_row_sum) * time_scale;
  }

  like += (sum(c) + sum(m)) * log(time_scale);

  List out = List::create(_["log.likelihood"] = like , _["likelihood"] = exp(like));

  return out;
}
