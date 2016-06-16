#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include <Eigen/Eigenvalues> 

using namespace Rcpp;

using Eigen::Map;           // 'maps' rather than copies 
using Eigen::MatrixXd;      // variable size matrix, double precision
using Eigen::VectorXd;      // variable size vector, double precision

// [[Rcpp::export()]]
List eigen_RcppEigen(const Eigen::Map<Eigen::MatrixXd> & mat)
{
  // call the eigen solver `EigenSolver` @ https://eigen.tuxfamily.org/dox/classEigen_1_1EigenSolver.html
  Eigen::EigenSolver<MatrixXd> es(mat);

  // return
  List ret;

  ret["solver"] = "EigenSolver";
  ret["values"] = es.eigenvalues();
  ret["vectors"] = es.eigenvectors();
   
  return(ret);
}

// [[Rcpp::export()]]
List eigenSelfAdjoint_RcppEigen(const Eigen::Map<Eigen::MatrixXd> & mat)
{
  // call the eigen solver `SelfAdjointEigenSolver` @ https://eigen.tuxfamily.org/dox/classEigen_1_1SelfAdjointEigenSolver.html
  Eigen::SelfAdjointEigenSolver<MatrixXd> es(mat);

  // return
  List ret;

  ret["solver"] = "SelfAdjointEigenSolver";
  ret["values"] = es.eigenvalues();
  ret["vectors"] = es.eigenvectors();
   
  return(ret);
}

