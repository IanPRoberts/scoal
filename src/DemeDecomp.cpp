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

std::list<int> DemeDecompC(NumericMatrix ED, int n_deme, NumericVector node_indices) {
  int nrow = ED.nrow();

  NumericVector event_times = ED(_,5);
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
  std::list<int>::iterator node;  //Initialise iterator "node" to iterate over list objects
  double current_time;
  int current_deme;
  int current_row;

  for (int i = 2; i < 4; ++i){
    active_nodes.emplace_back(ED(root_row,i));
  }

  for (int i = 1; i < k.nrow(); ++i){
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
    if (current_nodes.size() > 1){
      for (node = current_nodes.begin(); node != current_nodes.end(); ++node){
        current_deme = ED(node_indices[*node - 1] - 1,4);
        k(i, current_deme) -= 1;
      }
    } else{
      // TO DO STILL....
    }
  }

  return active_nodes;
}


//   for (i in 2 : (length(event.times) - 1)){
//     k[i,] <- k[i-1,]
//     active.rows <- node.indices[active.nodes]
//     current.indices <- which(ED[active.rows, 6] == event.times[i])
//     current.rows <- active.rows[current.indices]
//
//     if (length(current.rows) > 1){ #Multiple leaves
//       for (j in current.rows){
//         current.deme <- ED[j, 5]
//         k[i, current.deme] <- k[i, current.deme] - 1
//       }
//     } else{
//       if (anyNA(ED[current.rows, 3])){ #Leaf
//         current.deme <- ED[current.rows, 5]
//         k[i, current.deme] <- k[i, current.deme] - 1
//       } else if (!anyNA(ED[current.rows, 4])){ #Coalescence
//         current.deme <- ED[current.rows, 5]
//         k[i, current.deme] <- k[i, current.deme] + 1
//         active.nodes <- c(active.nodes, ED[current.rows, 3:4])
//       } else{ #Migration
//         current.deme <- ED[current.rows, 5]
//         k[i, current.deme] <- k[i, current.deme] - 1
//         current.child <- node.indices[ED[current.rows, 3]]
//         child.deme <- ED[current.child, 5]
//         k[i, child.deme] <- k[i, child.deme] + 1
//         active.nodes <- c(active.nodes, ED[current.rows, 3])
//       }
//     }
//     active.nodes <- active.nodes[-current.indices]
//   }
//
//   return(k)
// }
