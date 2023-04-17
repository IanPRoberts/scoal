// [[Rcpp::interfaces(r, cpp)]]
#include <Rcpp.h>
#include "DemeDecomp.h"
#include "Node_indices.h"

using namespace Rcpp;

//' @title ED_Dist
//' @description Computes the distance between pairs of migration histories on a common topology via proportion of branch lengths different between the two
//' @param ED1 NumericMatrix Extended data object representing first migration history
//' @param ED2 NumericMatrix Extended data object representing second migration history
//' @param node_indices_1 NumericVector Vector giving row numbers in ED1 for node labels
//' @param node_indices_2 NumericVector Vector giving row numbers in ED1 for node labels
//' @returns double Proportion of tree in different demes between ED1 and ED2
//'
//' @export
//'
// [[Rcpp::export]]

double ED_dist_C(NumericMatrix ED1,
                 NumericMatrix ED2,
                 int n_deme,
                 Nullable<NumericVector> node_indices_1 = R_NilValue,
                 Nullable<NumericVector> node_indices_2 = R_NilValue) {

  // Compute node_indices if not given in function call
  NumericVector ni_1;
  NumericVector ni_2;
  if (node_indices_1.isNull()){
    ni_1 = NodeIndicesC(ED1);
  } else{
    ni_1 = node_indices_1.get();
  }

  if (node_indices_2.isNull()){
    ni_2 = NodeIndicesC(ED2);
  } else{
    ni_2 = node_indices_2.get();
  }

  // Pooled event times (sorted & unique)
  List dd_1 = DemeDecompC(ED1, n_deme, ni_1);
  List dd_2 = DemeDecompC(ED2, n_deme, ni_2);

  NumericVector et_1 = dd_1["event.times"];
  NumericVector et_2 = dd_2["event.times"];

  NumericVector pool_et = et_1;
  for (int i = 0; i < et_2.size(); ++i){
    pool_et.push_back(et_2[i]);
  }
  NumericVector pool_et_su = sort_unique(pool_et);

  NumericVector pool_incs = diff(pool_et_su);
  int n_event = pool_incs.size();

  NumericMatrix k_1 = dd_1["k"];
  NumericMatrix k_2 = dd_2["k"];
  NumericMatrix k_diff(n_event, n_deme);
  double score = 0;

  for (int i = 0, count_1 = 0, count_2 = 0; i < n_event; ++i){

    for (int j = 0; j < n_deme; ++j){
      k_diff(i,j) = abs(k_1(count_1,j) - k_2(count_2,j));
    }

    score += sum(k_diff(i,_) * pool_incs[i]) / 2;

    if (pool_et_su[i+1] == et_1[count_1+1]){
      ++count_1;
    }
    if (pool_et_su[i+1] == et_2[count_2+1]){
      ++count_2;
    }
  }

  double tree_len = 0;
  NumericVector ti_1 = dd_1["time.increments"];
  n_event = ti_1.size();

  for (int i = 0; i < n_event; ++i){
    tree_len += sum(k_1(i,_) * ti_1[i]);
  }

  score /= tree_len;

  return score;
}
