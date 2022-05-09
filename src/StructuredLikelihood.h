#ifndef StructuredLikelihoodC_H
#define StructuredLikelihoodC_H

#include <RcppCommon.h>
using namespace Rcpp;

List StructuredLikelihoodC(NumericMatrix ED, NumericVector eff_pop, double gen_len, NumericMatrix mig_mat, NumericVector node_indices);

#endif // RCPP_StructuredLikelihoodC_H_GEN_
