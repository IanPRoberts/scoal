#ifndef DemeDecompC_H
#define DemeDecompC_H

#include <RcppCommon.h>
using namespace Rcpp;

List DemeDecomp(NumericMatrix ED, int n_deme, NumericVector node_indices);

#endif // RCPP_DemeDecompC_H_GEN_
