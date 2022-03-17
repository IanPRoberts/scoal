#include <Rcpp.h>
using namespace Rcpp;

//' @title DemeDecompC
//' @description Decomposes an extended data object to give number of lineages in each deme during each time increment
//' @param ED NumericMatrix Extended data object representing structured phylogeny
//' @param n_deme int Number of demes modelled under ED
//' @param node_indices NumericVector Vector giving row numbers in ED for node labels
//' @returns List containing matrix k giving deme decomposition with number of lineages in each deme during each time interval, vector event.times containing sorted event times from phylogeny and vector time.increments giving time increments between events during phylogeny
//'
//' @export
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

List DemeDecompC(NumericMatrix ED, int n_deme, NumericVector node_indices) {
  int nrow = ED.nrow();
  NumericVector event_times = sort_unique(ED(_, 5));
  NumericVector time_increments = diff(event_times);
  NumericMatrix k(time_increments.length(), n_deme);
  int root_row;

  for (int i = 0; i < nrow; ++i) { //Find root row
    if(NumericVector::is_na(ED(i,1))){
      root_row = i;
      break; //Exit loop once (unique) root found
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

  for (int i = 2; i < 4; ++i){  //Add children of root to active_nodes
    active_nodes.emplace_back(ED(root_row,i));
  }

  for (int i = 1; i < k.nrow(); ++i){
    k(i,_) = k(i-1,_);
    current_nodes.clear();
    current_time = event_times[i];

    for (node = active_nodes.begin(); node != active_nodes.end(); ++node){ //Construct current_nodes from active_nodes
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
        for (int j = 2; j < 4; ++j){ //Add children of current_nodes to active_nodes
          active_nodes.emplace_back(ED(current_row, j));
        }
      } else{ //Migration
        child_row = node_indices[ED(current_row, 2) - 1] - 1;
        child_deme = ED(child_row, 4) - 1;

        k(i, current_deme) -= 1;
        k(i, child_deme) += 1;
        active_nodes.emplace_back(ED(child_row, 0)); //Add child of current_nodes to active_nodes
      }
    }
  }
  List out = List::create(_["k"] = k , _["event.times"] = event_times, _["time.increments"] = time_increments);

  return out;
}
