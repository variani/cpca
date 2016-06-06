// @ http://stackoverflow.com/questions/26234055/cohabitation-of-rcpparmadillo-and-rcppparallel

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

using namespace Rcpp;

// Scenario 1
//#define vector NumericVector

// Scenario 2
#define vector arma::vec

// compute values i/i+1 for i = 0 to n-1
// [[Rcpp::export]]
vector f1(int n) {
  vector x(n);
  for(int i = 0; i < n; i++) x[i] = (double) i/ (i+1);
  return x;
}

struct mytry : public Worker {
  vector& output;

  mytry(vector& output) : 
    output(output) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(int i = begin; i < end; i++) {
      output[i] = (double) i/ (i+1);
    }
  }
};

// [[Rcpp::export]]
vector f2(int n) 
{
  vector x(n);

  mytry A(x);

  parallelFor(0, n, A);

  return x;
}  
