// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/scoal.h"
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// DemeDecompC
List DemeDecompC(NumericMatrix ED, int n_deme, NumericVector node_indices);
static SEXP _scoal_DemeDecompC_try(SEXP EDSEXP, SEXP n_demeSEXP, SEXP node_indicesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ED(EDSEXP);
    Rcpp::traits::input_parameter< int >::type n_deme(n_demeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type node_indices(node_indicesSEXP);
    rcpp_result_gen = Rcpp::wrap(DemeDecompC(ED, n_deme, node_indices));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _scoal_DemeDecompC(SEXP EDSEXP, SEXP n_demeSEXP, SEXP node_indicesSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_scoal_DemeDecompC_try(EDSEXP, n_demeSEXP, node_indicesSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// ED_dist_C
double ED_dist_C(NumericMatrix ED1, NumericMatrix ED2, int n_deme, Nullable<NumericVector> node_indices_1, Nullable<NumericVector> node_indices_2);
static SEXP _scoal_ED_dist_C_try(SEXP ED1SEXP, SEXP ED2SEXP, SEXP n_demeSEXP, SEXP node_indices_1SEXP, SEXP node_indices_2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ED1(ED1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ED2(ED2SEXP);
    Rcpp::traits::input_parameter< int >::type n_deme(n_demeSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type node_indices_1(node_indices_1SEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type node_indices_2(node_indices_2SEXP);
    rcpp_result_gen = Rcpp::wrap(ED_dist_C(ED1, ED2, n_deme, node_indices_1, node_indices_2));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _scoal_ED_dist_C(SEXP ED1SEXP, SEXP ED2SEXP, SEXP n_demeSEXP, SEXP node_indices_1SEXP, SEXP node_indices_2SEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_scoal_ED_dist_C_try(ED1SEXP, ED2SEXP, n_demeSEXP, node_indices_1SEXP, node_indices_2SEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// mcmc_cpp
void mcmc_cpp(int N0, int N, NumericMatrix ED, NumericVector coal_rate, double time_scale, NumericMatrix mig_mat, int n_deme, NumericVector prop_rates, double cr_prior_shape, double cr_prior_rate, double mm_prior_shape, double mm_prior_rate);
static SEXP _scoal_mcmc_cpp_try(SEXP N0SEXP, SEXP NSEXP, SEXP EDSEXP, SEXP coal_rateSEXP, SEXP time_scaleSEXP, SEXP mig_matSEXP, SEXP n_demeSEXP, SEXP prop_ratesSEXP, SEXP cr_prior_shapeSEXP, SEXP cr_prior_rateSEXP, SEXP mm_prior_shapeSEXP, SEXP mm_prior_rateSEXP) {
BEGIN_RCPP
    Rcpp::traits::input_parameter< int >::type N0(N0SEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ED(EDSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type coal_rate(coal_rateSEXP);
    Rcpp::traits::input_parameter< double >::type time_scale(time_scaleSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mig_mat(mig_matSEXP);
    Rcpp::traits::input_parameter< int >::type n_deme(n_demeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prop_rates(prop_ratesSEXP);
    Rcpp::traits::input_parameter< double >::type cr_prior_shape(cr_prior_shapeSEXP);
    Rcpp::traits::input_parameter< double >::type cr_prior_rate(cr_prior_rateSEXP);
    Rcpp::traits::input_parameter< double >::type mm_prior_shape(mm_prior_shapeSEXP);
    Rcpp::traits::input_parameter< double >::type mm_prior_rate(mm_prior_rateSEXP);
    mcmc_cpp(N0, N, ED, coal_rate, time_scale, mig_mat, n_deme, prop_rates, cr_prior_shape, cr_prior_rate, mm_prior_shape, mm_prior_rate);
    return R_NilValue;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _scoal_mcmc_cpp(SEXP N0SEXP, SEXP NSEXP, SEXP EDSEXP, SEXP coal_rateSEXP, SEXP time_scaleSEXP, SEXP mig_matSEXP, SEXP n_demeSEXP, SEXP prop_ratesSEXP, SEXP cr_prior_shapeSEXP, SEXP cr_prior_rateSEXP, SEXP mm_prior_shapeSEXP, SEXP mm_prior_rateSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_scoal_mcmc_cpp_try(N0SEXP, NSEXP, EDSEXP, coal_rateSEXP, time_scaleSEXP, mig_matSEXP, n_demeSEXP, prop_ratesSEXP, cr_prior_shapeSEXP, cr_prior_rateSEXP, mm_prior_shapeSEXP, mm_prior_rateSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// FitMigMatC
NumericMatrix FitMigMatC(NumericMatrix bit_mm, NumericVector coal_rate);
static SEXP _scoal_FitMigMatC_try(SEXP bit_mmSEXP, SEXP coal_rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type bit_mm(bit_mmSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type coal_rate(coal_rateSEXP);
    rcpp_result_gen = Rcpp::wrap(FitMigMatC(bit_mm, coal_rate));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _scoal_FitMigMatC(SEXP bit_mmSEXP, SEXP coal_rateSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_scoal_FitMigMatC_try(bit_mmSEXP, coal_rateSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// BitMigMatC
NumericMatrix BitMigMatC(NumericMatrix fit_mm, NumericVector coal_rate);
static SEXP _scoal_BitMigMatC_try(SEXP fit_mmSEXP, SEXP coal_rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type fit_mm(fit_mmSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type coal_rate(coal_rateSEXP);
    rcpp_result_gen = Rcpp::wrap(BitMigMatC(fit_mm, coal_rate));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _scoal_BitMigMatC(SEXP fit_mmSEXP, SEXP coal_rateSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_scoal_BitMigMatC_try(fit_mmSEXP, coal_rateSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// NodeCountC
List NodeCountC(NumericMatrix ED, int n_deme, NumericVector node_indices);
static SEXP _scoal_NodeCountC_try(SEXP EDSEXP, SEXP n_demeSEXP, SEXP node_indicesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ED(EDSEXP);
    Rcpp::traits::input_parameter< int >::type n_deme(n_demeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type node_indices(node_indicesSEXP);
    rcpp_result_gen = Rcpp::wrap(NodeCountC(ED, n_deme, node_indices));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _scoal_NodeCountC(SEXP EDSEXP, SEXP n_demeSEXP, SEXP node_indicesSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_scoal_NodeCountC_try(EDSEXP, n_demeSEXP, node_indicesSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// NodeIndicesC
NumericVector NodeIndicesC(NumericMatrix ED);
static SEXP _scoal_NodeIndicesC_try(SEXP EDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ED(EDSEXP);
    rcpp_result_gen = Rcpp::wrap(NodeIndicesC(ED));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _scoal_NodeIndicesC(SEXP EDSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_scoal_NodeIndicesC_try(EDSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// ScaledDTALikelihoodC
List ScaledDTALikelihoodC(NumericMatrix ED, NumericVector coal_rate, double time_scale, NumericMatrix bit_mig_mat, NumericVector node_indices);
static SEXP _scoal_ScaledDTALikelihoodC_try(SEXP EDSEXP, SEXP coal_rateSEXP, SEXP time_scaleSEXP, SEXP bit_mig_matSEXP, SEXP node_indicesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ED(EDSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type coal_rate(coal_rateSEXP);
    Rcpp::traits::input_parameter< double >::type time_scale(time_scaleSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bit_mig_mat(bit_mig_matSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type node_indices(node_indicesSEXP);
    rcpp_result_gen = Rcpp::wrap(ScaledDTALikelihoodC(ED, coal_rate, time_scale, bit_mig_mat, node_indices));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _scoal_ScaledDTALikelihoodC(SEXP EDSEXP, SEXP coal_rateSEXP, SEXP time_scaleSEXP, SEXP bit_mig_matSEXP, SEXP node_indicesSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_scoal_ScaledDTALikelihoodC_try(EDSEXP, coal_rateSEXP, time_scaleSEXP, bit_mig_matSEXP, node_indicesSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// ScaledLikelihoodC
List ScaledLikelihoodC(NumericMatrix ED, NumericVector coal_rate, double time_scale, NumericMatrix mig_mat, NumericVector node_indices);
static SEXP _scoal_ScaledLikelihoodC_try(SEXP EDSEXP, SEXP coal_rateSEXP, SEXP time_scaleSEXP, SEXP mig_matSEXP, SEXP node_indicesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ED(EDSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type coal_rate(coal_rateSEXP);
    Rcpp::traits::input_parameter< double >::type time_scale(time_scaleSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mig_mat(mig_matSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type node_indices(node_indicesSEXP);
    rcpp_result_gen = Rcpp::wrap(ScaledLikelihoodC(ED, coal_rate, time_scale, mig_mat, node_indices));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _scoal_ScaledLikelihoodC(SEXP EDSEXP, SEXP coal_rateSEXP, SEXP time_scaleSEXP, SEXP mig_matSEXP, SEXP node_indicesSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_scoal_ScaledLikelihoodC_try(EDSEXP, coal_rateSEXP, time_scaleSEXP, mig_matSEXP, node_indicesSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// StructuredLikelihoodC
List StructuredLikelihoodC(NumericMatrix ED, NumericVector eff_pop, double gen_len, NumericMatrix mig_mat, NumericVector node_indices);
static SEXP _scoal_StructuredLikelihoodC_try(SEXP EDSEXP, SEXP eff_popSEXP, SEXP gen_lenSEXP, SEXP mig_matSEXP, SEXP node_indicesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ED(EDSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eff_pop(eff_popSEXP);
    Rcpp::traits::input_parameter< double >::type gen_len(gen_lenSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mig_mat(mig_matSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type node_indices(node_indicesSEXP);
    rcpp_result_gen = Rcpp::wrap(StructuredLikelihoodC(ED, eff_pop, gen_len, mig_mat, node_indices));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _scoal_StructuredLikelihoodC(SEXP EDSEXP, SEXP eff_popSEXP, SEXP gen_lenSEXP, SEXP mig_matSEXP, SEXP node_indicesSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_scoal_StructuredLikelihoodC_try(EDSEXP, eff_popSEXP, gen_lenSEXP, mig_matSEXP, node_indicesSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}

// validate (ensure exported C++ functions exist before calling them)
static int _scoal_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("List(*DemeDecompC)(NumericMatrix,int,NumericVector)");
        signatures.insert("double(*ED_dist_C)(NumericMatrix,NumericMatrix,int,Nullable<NumericVector>,Nullable<NumericVector>)");
        signatures.insert("void(*mcmc_cpp)(int,int,NumericMatrix,NumericVector,double,NumericMatrix,int,NumericVector,double,double,double,double)");
        signatures.insert("NumericMatrix(*FitMigMatC)(NumericMatrix,NumericVector)");
        signatures.insert("NumericMatrix(*BitMigMatC)(NumericMatrix,NumericVector)");
        signatures.insert("List(*NodeCountC)(NumericMatrix,int,NumericVector)");
        signatures.insert("NumericVector(*NodeIndicesC)(NumericMatrix)");
        signatures.insert("List(*ScaledDTALikelihoodC)(NumericMatrix,NumericVector,double,NumericMatrix,NumericVector)");
        signatures.insert("List(*ScaledLikelihoodC)(NumericMatrix,NumericVector,double,NumericMatrix,NumericVector)");
        signatures.insert("List(*StructuredLikelihoodC)(NumericMatrix,NumericVector,double,NumericMatrix,NumericVector)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _scoal_RcppExport_registerCCallable() { 
    R_RegisterCCallable("scoal", "_scoal_DemeDecompC", (DL_FUNC)_scoal_DemeDecompC_try);
    R_RegisterCCallable("scoal", "_scoal_ED_dist_C", (DL_FUNC)_scoal_ED_dist_C_try);
    R_RegisterCCallable("scoal", "_scoal_mcmc_cpp", (DL_FUNC)_scoal_mcmc_cpp_try);
    R_RegisterCCallable("scoal", "_scoal_FitMigMatC", (DL_FUNC)_scoal_FitMigMatC_try);
    R_RegisterCCallable("scoal", "_scoal_BitMigMatC", (DL_FUNC)_scoal_BitMigMatC_try);
    R_RegisterCCallable("scoal", "_scoal_NodeCountC", (DL_FUNC)_scoal_NodeCountC_try);
    R_RegisterCCallable("scoal", "_scoal_NodeIndicesC", (DL_FUNC)_scoal_NodeIndicesC_try);
    R_RegisterCCallable("scoal", "_scoal_ScaledDTALikelihoodC", (DL_FUNC)_scoal_ScaledDTALikelihoodC_try);
    R_RegisterCCallable("scoal", "_scoal_ScaledLikelihoodC", (DL_FUNC)_scoal_ScaledLikelihoodC_try);
    R_RegisterCCallable("scoal", "_scoal_StructuredLikelihoodC", (DL_FUNC)_scoal_StructuredLikelihoodC_try);
    R_RegisterCCallable("scoal", "_scoal_RcppExport_validate", (DL_FUNC)_scoal_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_scoal_DemeDecompC", (DL_FUNC) &_scoal_DemeDecompC, 3},
    {"_scoal_ED_dist_C", (DL_FUNC) &_scoal_ED_dist_C, 5},
    {"_scoal_mcmc_cpp", (DL_FUNC) &_scoal_mcmc_cpp, 12},
    {"_scoal_FitMigMatC", (DL_FUNC) &_scoal_FitMigMatC, 2},
    {"_scoal_BitMigMatC", (DL_FUNC) &_scoal_BitMigMatC, 2},
    {"_scoal_NodeCountC", (DL_FUNC) &_scoal_NodeCountC, 3},
    {"_scoal_NodeIndicesC", (DL_FUNC) &_scoal_NodeIndicesC, 1},
    {"_scoal_ScaledDTALikelihoodC", (DL_FUNC) &_scoal_ScaledDTALikelihoodC, 5},
    {"_scoal_ScaledLikelihoodC", (DL_FUNC) &_scoal_ScaledLikelihoodC, 5},
    {"_scoal_StructuredLikelihoodC", (DL_FUNC) &_scoal_StructuredLikelihoodC, 5},
    {"_scoal_RcppExport_registerCCallable", (DL_FUNC) &_scoal_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_scoal(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
