// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_GKF_RCPPEXPORTS_H_GEN_
#define RCPP_GKF_RCPPEXPORTS_H_GEN_

#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace GKF {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("GKF", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("GKF", "_GKF_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in GKF");
            }
        }
    }

    inline Rcpp::List nearPD(arma::mat x_, arma::vec w_, bool corr = false, bool keepDiag = false, bool EnforcePosDef = true, bool doSym = false, bool ensureSymmetry = false, double eig_tol = 1e-6, double conv_tol = 1e-7, double posd_tol = 1e-3, unsigned maxit = 100, bool chol = true) {
        typedef SEXP(*Ptr_nearPD)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_nearPD p_nearPD = NULL;
        if (p_nearPD == NULL) {
            validateSignature("Rcpp::List(*nearPD)(arma::mat,arma::vec,bool,bool,bool,bool,bool,double,double,double,unsigned,bool)");
            p_nearPD = (Ptr_nearPD)R_GetCCallable("GKF", "_GKF_nearPD");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_nearPD(Shield<SEXP>(Rcpp::wrap(x_)), Shield<SEXP>(Rcpp::wrap(w_)), Shield<SEXP>(Rcpp::wrap(corr)), Shield<SEXP>(Rcpp::wrap(keepDiag)), Shield<SEXP>(Rcpp::wrap(EnforcePosDef)), Shield<SEXP>(Rcpp::wrap(doSym)), Shield<SEXP>(Rcpp::wrap(ensureSymmetry)), Shield<SEXP>(Rcpp::wrap(eig_tol)), Shield<SEXP>(Rcpp::wrap(conv_tol)), Shield<SEXP>(Rcpp::wrap(posd_tol)), Shield<SEXP>(Rcpp::wrap(maxit)), Shield<SEXP>(Rcpp::wrap(chol)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::List >(rcpp_result_gen);
    }

}

#endif // RCPP_GKF_RCPPEXPORTS_H_GEN_