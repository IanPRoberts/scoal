#include <Rcpp.h>
using namespace Rcpp;

//' @title DemeDecomp
//' @description Testing :)
//' @param x numeric vector blah blah
//' @returns mean(x) blah blah
//'
//' @export
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

NumericMatrix DemeDecompC(NumericMatrix ED, int n_deme, NumericVector node_indices) {
  int nrow = ED.nrow();
  NumericVector event_times = ED(_, 5);
  event_times = unique(event_times);
  std::sort(event_times.begin(), event_times.end());

  int n_event = event_times.length();
  NumericVector time_increments(n_event - 1);
  time_increments = tail(event_times, n_event - 1) - head(event_times, n_event - 1);

  NumericMatrix k(n_event - 1, n_deme);
  int root_row;

  for (int i = 0; i < nrow; ++i) {
    if(NumericVector::is_na(ED(i,1))){
      root_row = i;
      break;
    }
  }
  k(0, ED(root_row, 4) - 1) = 2;

  std::list<int> active_nodes;
  std::list<int> current_nodes;
  std::list<int>::iterator node;
  double current_time;
  int current_deme;
  int current_row;
  int child_row;
  int child_deme;

  for (int i = 2; i < 4; ++i){
    active_nodes.emplace_back(ED(root_row,i));
  }

  for (int i = 1; i < k.nrow(); ++i){
    //Printing active_nodes to console
    // Rprintf("Active nodes: ");
    // for (node = active_nodes.begin(); node != active_nodes.end(); ++ node){
    //   Rprintf("%d, ", *node);
    // }
    // Rprintf("\n");
    k(i,_) = k(i-1,_);
    current_nodes.clear();
    current_time = event_times[i];

    for (node = active_nodes.begin(); node != active_nodes.end(); ++node){
      // *node dereferences iterator "node" to get value at current entry
      if (ED(node_indices[*node - 1] - 1,5) == current_time){
        current_nodes.emplace_back(*node);
        active_nodes.erase(node);
      }
    }
    if (current_nodes.size() > 1){ // Multiple Leaves
      for (node = current_nodes.begin(); node != current_nodes.end(); ++node){
        current_deme = ED(node_indices[*node - 1] - 1, 4);
        k(i, current_deme - 1) -= 1;
      }
    } else{
      current_row = node_indices[*current_nodes.begin() - 1] - 1; // *current_nodes.begin() = dereferenced iterator giving the first (only) element in list; -1 to account for 0-counting
      current_deme = ED(current_row, 4) - 1;
      if (NumericMatrix::is_na(ED(current_row, 2)) == 1){ // Single leaf
        k(i, current_deme) -= 1;
      } else if (NumericMatrix::is_na(ED(current_row, 3)) == 0){ //Coalescence
        k(i, current_deme) += 1;
        for (int j = 2; j < 4; ++j){
          active_nodes.emplace_back(ED(current_row, j));
        }
      } else{ //Migration
        child_row = node_indices[ED(current_row, 2) - 1] - 1;
        child_deme = ED(child_row, 4) - 1;

        k(i, current_deme) -= 1;
        k(i, child_deme) += 1;
        active_nodes.emplace_back(ED(child_row, 0));
      }
    }
  }
  return k;
}
