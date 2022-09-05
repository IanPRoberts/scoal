#ifndef ScaledLikelihoodC_H
#define ScaledLikelihoodC_H

#include <RcppCommon.h>
using namespace Rcpp;

List ScaledLikelihoodC(NumericMatrix ED, NumericVector coal_rate, double time_scale, NumericMatrix mig_mat, NumericVector node_indices);

#endif // RCPP_ScaledLikelihoodC_H_GEN_
