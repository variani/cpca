# include <Rcpp.h>

using namespace Rcpp;

NumericVector numericProdMatVec(const NumericMatrix & mat, const NumericVector & vec);
double numericNorm(const NumericVector & x);
NumericVector numericMultVec(const NumericVector & x, double a);

// [[Rcpp::export()]]
List eigenPower_Rcpp(const NumericMatrix & A, const NumericVector & v0,
  const double tol = 1e-6, const int maxit = 1e3,
  const int verbose = 0)
{
  // inputs
  int N = A.nrow();
  int K = A.ncol();
  
  // containers
  NumericVector v(K);
  NumericVector b(K);
  double bnorm = 0;
        
  // algorithm
  double lambda0 = 0;
  double lambda = lambda0;
  double delta = 0;
  bool converged = false;
  
  // initialize & normalize `v`
  double v0norm = numericNorm(v0);
  v = numericMultVec(v0, 1.0 / v0norm);
  
  // loop  
  int it = 1;
  for(; it <= maxit ; it++) { 
    if(verbose > 1) { 
      Rcout << " * it " << it << std::endl;
    }
    
    b = numericProdMatVec(A, v);
    bnorm = numericNorm(b);
    v = numericMultVec(b, 1.0 / bnorm);
    
    lambda = bnorm;
        
    /*    
    // step 1: 
    // - product of `A` & `v` (to be stored in `b`)
    // - norm of the resulted product (`bnorm`)
    double sum = 0;
    for(int i = 0; i < N; i++) {
      double elem = 0;
      for(int j = 0; j < K; j++) {
        double val = A(i, j) * v(j);
        elem += val;
      }
      b(i) = elem;
      sum += elem * elem;
    }
    bnorm = sqrt(sum);
    
    // step 2: computed a new vector `v`
    for(int i = 0; i < N; i++) {
      v(i) = b(i) / bnorm;
    }
    
    // step 3: assigned `bnorm` to `lambda`
    lambda = bnorm;
    */
    
    if(verbose > 2) { 
      Rcout << "  -- lambda " << lambda << std::endl;
    }
      
    delta = fabs((lambda - lambda0) / lambda);
    if(delta < tol) {
      break;
    }
      
    lambda0 = lambda;
    
    Rcpp::checkUserInterrupt();
  }
    
  converged = (it < maxit);

  // return
  List ret;
  ret["v0"] = v0;
  ret["tol"] = tol;
  ret["maxit"] = maxit;
  ret["it"] = it;
  ret["delta"] = delta;
  ret["converged"] = converged;
    
  ret["lambda"] = lambda;
  ret["v"] = v;
   
  return(ret) ;
}
