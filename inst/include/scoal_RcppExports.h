// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_scoal_RCPPEXPORTS_H_GEN_
#define RCPP_scoal_RCPPEXPORTS_H_GEN_

#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace scoal {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("scoal", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("scoal", "_scoal_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in scoal");
            }
        }
    }

    inline List DemeDecompC(NumericMatrix ED, int n_deme, NumericVector node_indices) {
        typedef SEXP(*Ptr_DemeDecompC)(SEXP,SEXP,SEXP);
        static Ptr_DemeDecompC p_DemeDecompC = NULL;
        if (p_DemeDecompC == NULL) {
            validateSignature("List(*DemeDecompC)(NumericMatrix,int,NumericVector)");
            p_DemeDecompC = (Ptr_DemeDecompC)R_GetCCallable("scoal", "_scoal_DemeDecompC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_DemeDecompC(Shield<SEXP>(Rcpp::wrap(ED)), Shield<SEXP>(Rcpp::wrap(n_deme)), Shield<SEXP>(Rcpp::wrap(node_indices)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline double ED_dist_C(NumericMatrix ED1, NumericMatrix ED2, int n_deme, Nullable<NumericVector> node_indices_1 = R_NilValue, Nullable<NumericVector> node_indices_2 = R_NilValue) {
        typedef SEXP(*Ptr_ED_dist_C)(SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_ED_dist_C p_ED_dist_C = NULL;
        if (p_ED_dist_C == NULL) {
            validateSignature("double(*ED_dist_C)(NumericMatrix,NumericMatrix,int,Nullable<NumericVector>,Nullable<NumericVector>)");
            p_ED_dist_C = (Ptr_ED_dist_C)R_GetCCallable("scoal", "_scoal_ED_dist_C");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_ED_dist_C(Shield<SEXP>(Rcpp::wrap(ED1)), Shield<SEXP>(Rcpp::wrap(ED2)), Shield<SEXP>(Rcpp::wrap(n_deme)), Shield<SEXP>(Rcpp::wrap(node_indices_1)), Shield<SEXP>(Rcpp::wrap(node_indices_2)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline void mcmc_cpp(int N0, int N, NumericMatrix ED, NumericVector coal_rate, double time_scale, NumericMatrix mig_mat, int n_deme, NumericVector prop_rates, double cr_prior_shape, double cr_prior_rate, double mm_prior_shape, double mm_prior_rate) {
        typedef SEXP(*Ptr_mcmc_cpp)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_mcmc_cpp p_mcmc_cpp = NULL;
        if (p_mcmc_cpp == NULL) {
            validateSignature("void(*mcmc_cpp)(int,int,NumericMatrix,NumericVector,double,NumericMatrix,int,NumericVector,double,double,double,double)");
            p_mcmc_cpp = (Ptr_mcmc_cpp)R_GetCCallable("scoal", "_scoal_mcmc_cpp");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_mcmc_cpp(Shield<SEXP>(Rcpp::wrap(N0)), Shield<SEXP>(Rcpp::wrap(N)), Shield<SEXP>(Rcpp::wrap(ED)), Shield<SEXP>(Rcpp::wrap(coal_rate)), Shield<SEXP>(Rcpp::wrap(time_scale)), Shield<SEXP>(Rcpp::wrap(mig_mat)), Shield<SEXP>(Rcpp::wrap(n_deme)), Shield<SEXP>(Rcpp::wrap(prop_rates)), Shield<SEXP>(Rcpp::wrap(cr_prior_shape)), Shield<SEXP>(Rcpp::wrap(cr_prior_rate)), Shield<SEXP>(Rcpp::wrap(mm_prior_shape)), Shield<SEXP>(Rcpp::wrap(mm_prior_rate)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
    }

    inline NumericMatrix FitMigMatC(NumericMatrix bit_mm, NumericVector coal_rate) {
        typedef SEXP(*Ptr_FitMigMatC)(SEXP,SEXP);
        static Ptr_FitMigMatC p_FitMigMatC = NULL;
        if (p_FitMigMatC == NULL) {
            validateSignature("NumericMatrix(*FitMigMatC)(NumericMatrix,NumericVector)");
            p_FitMigMatC = (Ptr_FitMigMatC)R_GetCCallable("scoal", "_scoal_FitMigMatC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_FitMigMatC(Shield<SEXP>(Rcpp::wrap(bit_mm)), Shield<SEXP>(Rcpp::wrap(coal_rate)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericMatrix >(rcpp_result_gen);
    }

    inline NumericMatrix BitMigMatC(NumericMatrix fit_mm, NumericVector coal_rate) {
        typedef SEXP(*Ptr_BitMigMatC)(SEXP,SEXP);
        static Ptr_BitMigMatC p_BitMigMatC = NULL;
        if (p_BitMigMatC == NULL) {
            validateSignature("NumericMatrix(*BitMigMatC)(NumericMatrix,NumericVector)");
            p_BitMigMatC = (Ptr_BitMigMatC)R_GetCCallable("scoal", "_scoal_BitMigMatC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_BitMigMatC(Shield<SEXP>(Rcpp::wrap(fit_mm)), Shield<SEXP>(Rcpp::wrap(coal_rate)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericMatrix >(rcpp_result_gen);
    }

    inline List NodeCountC(NumericMatrix ED, int n_deme, NumericVector node_indices) {
        typedef SEXP(*Ptr_NodeCountC)(SEXP,SEXP,SEXP);
        static Ptr_NodeCountC p_NodeCountC = NULL;
        if (p_NodeCountC == NULL) {
            validateSignature("List(*NodeCountC)(NumericMatrix,int,NumericVector)");
            p_NodeCountC = (Ptr_NodeCountC)R_GetCCallable("scoal", "_scoal_NodeCountC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_NodeCountC(Shield<SEXP>(Rcpp::wrap(ED)), Shield<SEXP>(Rcpp::wrap(n_deme)), Shield<SEXP>(Rcpp::wrap(node_indices)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline NumericVector NodeIndicesC(NumericMatrix ED) {
        typedef SEXP(*Ptr_NodeIndicesC)(SEXP);
        static Ptr_NodeIndicesC p_NodeIndicesC = NULL;
        if (p_NodeIndicesC == NULL) {
            validateSignature("NumericVector(*NodeIndicesC)(NumericMatrix)");
            p_NodeIndicesC = (Ptr_NodeIndicesC)R_GetCCallable("scoal", "_scoal_NodeIndicesC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_NodeIndicesC(Shield<SEXP>(Rcpp::wrap(ED)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericVector >(rcpp_result_gen);
    }

    inline double SC_like_C(NumericMatrix ED, NumericVector coal_rate, NumericMatrix bit_mig_mat, NumericVector node_indices) {
        typedef SEXP(*Ptr_SC_like_C)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_SC_like_C p_SC_like_C = NULL;
        if (p_SC_like_C == NULL) {
            validateSignature("double(*SC_like_C)(NumericMatrix,NumericVector,NumericMatrix,NumericVector)");
            p_SC_like_C = (Ptr_SC_like_C)R_GetCCallable("scoal", "_scoal_SC_like_C");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_SC_like_C(Shield<SEXP>(Rcpp::wrap(ED)), Shield<SEXP>(Rcpp::wrap(coal_rate)), Shield<SEXP>(Rcpp::wrap(bit_mig_mat)), Shield<SEXP>(Rcpp::wrap(node_indices)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline List ScaledDTALikelihoodC(NumericMatrix ED, NumericVector coal_rate, double time_scale, NumericMatrix bit_mig_mat, NumericVector node_indices) {
        typedef SEXP(*Ptr_ScaledDTALikelihoodC)(SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_ScaledDTALikelihoodC p_ScaledDTALikelihoodC = NULL;
        if (p_ScaledDTALikelihoodC == NULL) {
            validateSignature("List(*ScaledDTALikelihoodC)(NumericMatrix,NumericVector,double,NumericMatrix,NumericVector)");
            p_ScaledDTALikelihoodC = (Ptr_ScaledDTALikelihoodC)R_GetCCallable("scoal", "_scoal_ScaledDTALikelihoodC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_ScaledDTALikelihoodC(Shield<SEXP>(Rcpp::wrap(ED)), Shield<SEXP>(Rcpp::wrap(coal_rate)), Shield<SEXP>(Rcpp::wrap(time_scale)), Shield<SEXP>(Rcpp::wrap(bit_mig_mat)), Shield<SEXP>(Rcpp::wrap(node_indices)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline List ScaledLikelihoodC(NumericMatrix ED, NumericVector coal_rate, double time_scale, NumericMatrix mig_mat, NumericVector node_indices) {
        typedef SEXP(*Ptr_ScaledLikelihoodC)(SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_ScaledLikelihoodC p_ScaledLikelihoodC = NULL;
        if (p_ScaledLikelihoodC == NULL) {
            validateSignature("List(*ScaledLikelihoodC)(NumericMatrix,NumericVector,double,NumericMatrix,NumericVector)");
            p_ScaledLikelihoodC = (Ptr_ScaledLikelihoodC)R_GetCCallable("scoal", "_scoal_ScaledLikelihoodC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_ScaledLikelihoodC(Shield<SEXP>(Rcpp::wrap(ED)), Shield<SEXP>(Rcpp::wrap(coal_rate)), Shield<SEXP>(Rcpp::wrap(time_scale)), Shield<SEXP>(Rcpp::wrap(mig_mat)), Shield<SEXP>(Rcpp::wrap(node_indices)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline List StructuredLikelihoodC(NumericMatrix ED, NumericVector eff_pop, double gen_len, NumericMatrix mig_mat, NumericVector node_indices) {
        typedef SEXP(*Ptr_StructuredLikelihoodC)(SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_StructuredLikelihoodC p_StructuredLikelihoodC = NULL;
        if (p_StructuredLikelihoodC == NULL) {
            validateSignature("List(*StructuredLikelihoodC)(NumericMatrix,NumericVector,double,NumericMatrix,NumericVector)");
            p_StructuredLikelihoodC = (Ptr_StructuredLikelihoodC)R_GetCCallable("scoal", "_scoal_StructuredLikelihoodC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_StructuredLikelihoodC(Shield<SEXP>(Rcpp::wrap(ED)), Shield<SEXP>(Rcpp::wrap(eff_pop)), Shield<SEXP>(Rcpp::wrap(gen_len)), Shield<SEXP>(Rcpp::wrap(mig_mat)), Shield<SEXP>(Rcpp::wrap(node_indices)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

}

#endif // RCPP_scoal_RCPPEXPORTS_H_GEN_
