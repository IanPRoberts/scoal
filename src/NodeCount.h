#ifndef NodeCountC_H
#define NodeCountC_H

#include <RcppCommon.h>
using namespace Rcpp;

List NodeCount(NumericMatrix ED, int n_deme, NumericVector node_indices);

#endif // RCPP_NodeCount_H_GEN_
