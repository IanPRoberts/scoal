#include <Rcpp.h>
using namespace Rcpp;
//'
//' @export
// [[Rcpp::export]]

int LoopTestC(int max) {
  int k=0;
  for (int i = 0; i < max; ++i) {
    k = i;
  }
  return k;
}
