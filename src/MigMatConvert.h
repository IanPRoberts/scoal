#ifndef MigMatConversion_H
#define MigMatConversion_H

#include <RcppCommon.h>
using namespace Rcpp;

NumericMatrix FitMigMatC(NumericMatrix bit_mm, NumericVector coal_rate);
NumericMatrix BitMigMatC(NumericMatrix fit_mm, NumericVector coal_rate);

#endif // RCPP_FITMMC_H_GEN_
