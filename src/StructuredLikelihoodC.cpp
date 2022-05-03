#include <Rcpp.h>
#include "DemeDecomp.h"
using namespace Rcpp;

//' @title NodeCountC
//' @description Counts the number of migrations between pairs of demes, and coalescences within each deme
//' @param ED NumericMatrix Extended data object representing structured phylogeny
//' @param n_deme int Number of demes modelled under ED
//' @param node_indices NumericVector Vector giving row numbers in ED for node labels
//' @returns List containing matrix m with element (i,j) giving the number of migrations i -> j backwards in time, and a vector c with element i giving the number of coalescences occurring in deme i
//'
//' @export
//'
// [[Rcpp::export]]

NumericMatrix StructuredLikelihoodC(NumericMatrix ED, NumericVector eff_pop, double gen_len, NumericMatrix mig_mat, NumericVector node_indices) {
  int n_deme = eff_pop.size();
  NumericVector lambda = eff_pop * gen_len;

  // Case if lambda.size() == 1
  List deme_decomp;
  deme_decomp = DemeDecompC(ED, n_deme, node_indices); // Import DemeDecompC to use within separate C++ file

  NumericMatrix k = deme_decomp["k"];

  return k;
}
