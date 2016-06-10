/* References
 - http://gallery.rcpp.org/articles/eigen-eigenvalues/
 - https://cran.r-project.org/web/packages/RcppEigen/vignettes/RcppEigen-Introduction.pdf
 - https://eigen.tuxfamily.org/dox/group__TopicAliasing.html
*/
 

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

using Eigen::Map;           // 'maps' rather than copies 
using Eigen::MatrixXd;      // variable size matrix, double precision
using Eigen::VectorXd;      // variable size vector, double precision


// [[Rcpp::export()]]
List eigenPower_RcppEigen(
  const Eigen::Map<Eigen::MatrixXd> & A, const Eigen::Map<Eigen::VectorXd> & v0,
  double tol = 1e-6, int maxit = 1e3,
  int verbose = 0)
{
  // inputs
  int N = A.rows();
  int K = A.cols();
  
  // containers
  VectorXd v(K);
  VectorXd b(K);
  double bnorm = 0;
        
  // algorithm
  double lambda0 = 0;
  double lambda = lambda0;
  double delta = 0;
  bool converged = false;

  // initialize & normalize `v`
  double v0norm = v0.norm();
  v.noalias() = v0 / v0norm;

  // loop  
  int it = 1;
  for(; it <= maxit ; it++) { 
    if(verbose > 1) { 
      Rcout << " * it " << it << std::endl;
    }
    
    b.noalias() = A * v;
    lambda = b.norm();
    v.noalias() = b / lambda;
  
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
