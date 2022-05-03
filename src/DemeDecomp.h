#ifndef DemeDecompC_H
#define DemeDecompC_H

#include <RcppCommon.h>
using namespace Rcpp;

List DemeDecompC(NumericMatrix ED, int n_deme, NumericVector node_indices);

#endif // RCPP_scoal_H_GEN_
