# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
  
// [[Rcpp::export()]]
List eigenPower_RcppArmadillo(const arma::mat & A, const arma::vec & v0,
  const double tol = 1e-6, const int maxit = 1e3,
  const int verbose = 0)
{
  // inputs
  const int N = A.n_rows;
  const int K = A.n_cols;

  // containers
  arma::vec v = v0;
  arma::vec b = v0;
  arma::vec v_lambda(1);
        
  // algorithm
  double lambda0 = 0;
  double lambda = lambda0;
  double delta = 0;
  bool converged = false;
    
  int it = 1;
  for ( ; it <= maxit ; it++) { 
    if(verbose > 1) { 
      Rcout << " * it " << it << std::endl;
    }
      
    b = A * v;
    v = normalise(b);
    v_lambda = trans(v) * b;
    lambda = v_lambda[0]; 
      
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
