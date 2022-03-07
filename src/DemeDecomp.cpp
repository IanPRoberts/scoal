#include <Rcpp.h>
using namespace Rcpp;

//' @title DemeDecomp
//' @description Testing :)
//' @param x numeric vector blah blah
//' @returns mean(x) blah blah
//'
//' @export
// [[Rcpp::export]]

NumericVector DemeDecompC(NumericMatrix ED, int n_deme, NumericVector node_indices) {
  int nrow = ED.nrow();
  NumericVector EventTimes(nrow);
  for (int i = 0; i < nrow; i++) {
    EventTimes[i] = ED(i,5);
  }
  EventTimes = unique(EventTimes);
  return EventTimes;  //Change function type to NumericVector before compiling
}
