// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

#define vector arma::vec
#define matrix arma::mat

struct ProdMatVecArma : public Worker 
{
  // inputs: two vectors
  const matrix& mat;
  const vector& vec;

  // output
  vector& rvec;
  
  // constructors
  ProdMatVecArma(const matrix& mat, const vector& vec, vector& rvec) 
    : mat(mat), vec(vec), rvec(rvec) {}
   
  // process just the elements of the range Ive been asked to
  void operator()(std::size_t begin, std::size_t end) 
  {
    cout << begin << ", " << end << endl;
    
    matrix mat_i = mat.rows(begin, end - 1);
    
    vector rvec_i = mat_i * vec;
    
    // assign
    for(std::size_t k = begin; k < end; k++) {
      std::size_t ki = k - begin;
      rvec[k] = rvec_i[ki];
    }
  }
};

// [[Rcpp::export]]
vector prodMatVecArmaParallel(matrix& mat, vector& vec, unsigned int size) 
{
  // allocate the matrix we will return
  vector rvec(mat.n_rows);
    
  // create the worker
  ProdMatVecArma prodMatVecArma(mat, vec, rvec);
     
  // call it with parallelFor
  parallelFor(0, mat.n_rows, prodMatVecArma, size);

  return rvec;
}
