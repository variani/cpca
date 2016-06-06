// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;


struct ProdMatVecEigen : public Worker 
{
  // inputs: two vectors
  const RMatrix<double> mat;
  const RVector<double> vec;

  // output
  RVector<double> rvec;
  
  // constructors
  ProdMatVecEigen(const Eigen::MatrixXd mat, const Eigen::VectorXd vec, Eigen::VectorXd rvec) 
    : mat(mat), vec(vec), rvec(rvec) {}
   
  // process just the elements of the range Ive been asked to
  void operator()(std::size_t begin, std::size_t end) 
  {
    ; //
  }
   
};

// [[Rcpp::export]]
Eigen::VectorXd ProdMatVecEigenParallel(Eigen::MatrixXd mat, Eigen::VectorXd vec) 
{
  // allocate the matrix we will return
  Eigen::VectorXd rvec(mat.nrow());
    
  // create the worker
  ProdMatVecEigen prodMatVecEigen(mat, vec, rvec);
     
  // call it with parallelFor
  parallelFor(0, mat.nrow(), prodMatVecEigen);

  return rvec;
}
