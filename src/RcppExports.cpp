// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// computeCindex
double computeCindex(NumericVector beta_6, NumericMatrix data_1);
RcppExport SEXP _SA23229012_computeCindex(SEXP beta_6SEXP, SEXP data_1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta_6(beta_6SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data_1(data_1SEXP);
    rcpp_result_gen = Rcpp::wrap(computeCindex(beta_6, data_1));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SA23229012_computeCindex", (DL_FUNC) &_SA23229012_computeCindex, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_SA23229012(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
