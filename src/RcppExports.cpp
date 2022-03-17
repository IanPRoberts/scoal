// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// DemeDecompC
List DemeDecompC(NumericMatrix ED, int n_deme, NumericVector node_indices);
RcppExport SEXP _scoal_DemeDecompC(SEXP EDSEXP, SEXP n_demeSEXP, SEXP node_indicesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ED(EDSEXP);
    Rcpp::traits::input_parameter< int >::type n_deme(n_demeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type node_indices(node_indicesSEXP);
    rcpp_result_gen = Rcpp::wrap(DemeDecompC(ED, n_deme, node_indices));
    return rcpp_result_gen;
END_RCPP
}
// Test
List Test(NumericMatrix ED, int n_deme, NumericVector node_indices);
RcppExport SEXP _scoal_Test(SEXP EDSEXP, SEXP n_demeSEXP, SEXP node_indicesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ED(EDSEXP);
    Rcpp::traits::input_parameter< int >::type n_deme(n_demeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type node_indices(node_indicesSEXP);
    rcpp_result_gen = Rcpp::wrap(Test(ED, n_deme, node_indices));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_scoal_DemeDecompC", (DL_FUNC) &_scoal_DemeDecompC, 3},
    {"_scoal_Test", (DL_FUNC) &_scoal_Test, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_scoal(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
