#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
  
// [[Rcpp::export()]]
List eigenPower_RcppEigen(const Eigen::MatrixXd & A, const Eigen::VectorXd v0,
  const double tol = 1e-6, const int maxit = 1e3,
  const int verbose = 0)
{
  // containers
  Eigen::VectorXd v = v0;
  Eigen::VectorXd b = v0;
 
  // algorithm
  double lambda0 = 0;
  double lambda = lambda0;
  double delta = 0;
  bool converged = false;
  
  // loop
  int it = 1;
  for(; it <= maxit ; it++) { 
    if(verbose > 1) { 
      Rcout << " * it " << it << std::endl;
    }
      
    b.noalias() = A * v; // @ https://eigen.tuxfamily.org/dox/group__TopicAliasing.html
    lambda = b.norm();
    v = b / lambda;

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
