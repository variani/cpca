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
using Eigen::VectorXi;

// [[Rcpp::export()]]
List eigenPower_RcppEigen(
  const Eigen::Map<Eigen::MatrixXd> & A, const Eigen::Map<Eigen::VectorXd> & v0,
  double tol = 1e-6, unsigned int maxit = 1e3, 
  unsigned int ncomp = 1,
  bool symmetric = false,
  unsigned int verbose = 0)
{
  // inputs
  unsigned int N = A.rows();
  unsigned int K = A.cols();
  
  if(ncomp == 0) {
    ncomp = K;
  }
  
  // containers
  VectorXd v(K);
  VectorXd b(K);
  double bnorm = 0;
        
  // algorithm
  double lambda0 = 0;
  double lambda = lambda0;
  double delta = 0;

  bool converged = false;

  VectorXi itComp = VectorXi::Zero(ncomp);
  VectorXi convergedComp = VectorXi::Zero(ncomp);
  
  VectorXd D = VectorXd::Zero(ncomp); 
  MatrixXd CPC = MatrixXd::Zero(K, ncomp);
  MatrixXd Qw = MatrixXd::Identity(K, K);
  
  // initialize & normalize `v`
  for(unsigned int comp = 0; comp < ncomp; comp++) {
    if(verbose > 1) { 
      Rcout << " * * component: " << comp + 1 << std::endl;
    }
  
    double v0norm = v0.norm();
    v.noalias() = v0 / v0norm;

    // loop  
    unsigned int it = 1;
    for(; it <= maxit ; it++) { 
      Rcpp::checkUserInterrupt();
    
      if(verbose > 1) { 
        Rcout << " * it " << it << std::endl;
      }
    
      if(symmetric) b.noalias() = A * v;
      else          b.noalias() = A.transpose() * v;
      
      if(comp > 0) { 
        b = Qw.transpose() * b; // `Qw` is symmetric
      }
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
    }
    
    // post-process
    itComp(comp) = it;
    if(it < maxit) {
      convergedComp(comp) = 1;
    }
    
    D(comp) = lambda;
    CPC.col(comp) = v;
     
    Qw -= v * v.transpose();
  }
  
  // return
  List ret;
  ret["v0"] = v0;
  ret["tol"] = tol;
  ret["maxit"] = maxit;
  ret["itComp"] = itComp;
  ret["convergedComp"] = convergedComp;
  
  ret["ncomp"] = ncomp;
  ret["values"] = D;
  ret["vectors"] = CPC;
  
  ret["lambda"] = D(0);
  ret["v"] = CPC.col(0);
   
  return(ret);
}
