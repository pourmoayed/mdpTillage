// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// SolveMDPModel
SEXP SolveMDPModel(const List paramModel);
RcppExport SEXP mdpTillage_SolveMDPModel(SEXP paramModelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type paramModel(paramModelSEXP);
    rcpp_result_gen = Rcpp::wrap(SolveMDPModel(paramModel));
    return rcpp_result_gen;
END_RCPP
}
