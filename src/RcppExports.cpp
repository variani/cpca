// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// eigenPower_Rcpp
List eigenPower_Rcpp(const arma::mat A, arma::vec v0, const double tol, const int maxit, const int verbose);
RcppExport SEXP cpca_eigenPower_Rcpp(SEXP ASEXP, SEXP v0SEXP, SEXP tolSEXP, SEXP maxitSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type v0(v0SEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< const int >::type verbose(verboseSEXP);
    __result = Rcpp::wrap(eigenPower_Rcpp(A, v0, tol, maxit, verbose));
    return __result;
END_RCPP
}
