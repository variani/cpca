// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

using namespace Rcpp;

//-----------------------------
// Inner Product in Parallel
//-----------------------------

struct InnerProduct : public Worker 
{
  // inputs: two vectors
  const RVector<double> x;
  const RVector<double> y;

  // output
  double product;
  
  // constructors
  InnerProduct(const NumericVector x, const NumericVector y) 
    : x(x), y(y), product(0) {}
  InnerProduct(const InnerProduct& innerProduct, Split) 
    : x(innerProduct.x), y(innerProduct.y), product(0) {}
   
  // process just the elements of the range Ive been asked to
  void operator()(std::size_t begin, std::size_t end) 
  {
    product += std::inner_product(x.begin() + begin, x.begin() + end, 
      y.begin() + begin, 0.0);
  }
   
  // join my value with that of another InnerProduct
  void join(const InnerProduct& rhs) 
  { 
    product += rhs.product; 
  }
};

// [[Rcpp::export]]
double innerProductParallel(NumericVector x, NumericVector y, unsigned int chunkSize = 1) 
{
  // declare the InnerProduct instance that takes a pointer to the vector data
  InnerProduct innerProduct(x, y);
   
  // call paralleReduce to start the work
  parallelReduce(0, x.length(), innerProduct, chunkSize);
   
  // return the computed product
  return innerProduct.product;
}

//-----------------------------
// Norm in Parallel
//-----------------------------

struct InnerNorm : public Worker 
{
  // inputs: two vectors
  const RVector<double> x;

  // output
  double norm;
  
  // constructors
  InnerNorm(const NumericVector x) 
    : x(x), norm(0) {}
  InnerNorm(const InnerNorm& innerNorm, Split) 
    : x(innerNorm.x), norm(0) {}
   
  // process just the elements of the range Ive been asked to
  void operator()(std::size_t begin, std::size_t end) 
  {
    norm += std::inner_product(x.begin() + begin, x.begin() + end, 
      x.begin() + begin, 0.0);
  }
   
  // join my value with that of another InnerProduct
  void join(const InnerNorm& rhs) 
  { 
    norm += rhs.norm; 
  }
};

// [[Rcpp::export]]
double innerNormParallel(NumericVector x, unsigned int chunkSize = 1) 
{
  // declare the InnerProduct instance that takes a pointer to the vector data
  InnerNorm innerNorm(x);
   
  // call paralleReduce to start the work
  parallelReduce(0, x.length(), innerNorm, chunkSize);
   
  // return the computed product
  return std::sqrt(innerNorm.norm);
}

//-----------------------------
// ProdMatVec in Parallel
//-----------------------------

struct ProdMatVec : public Worker 
{
  // input to read from
  const RMatrix<double> mat;
  const RVector<double> vec;
   
  // output to write to
  RVector<double> rvec;
   
  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  ProdMatVec(const NumericMatrix mat, const NumericVector vec, NumericVector rvec)
    : mat(mat), vec(vec), rvec(rvec) {}
    
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) 
  {
    for (std::size_t i = begin; i < end; i++) {
      // rows we will operate on
      RMatrix<double>::Row rowi = mat.row(i);
        
      // calculate inneer product
      double prod = std::inner_product(rowi.begin(), rowi.end(), vec.begin(), 0.0);
        
      // write to output matrix
      rvec[i] = prod;
    }
  }
};

// [[Rcpp::export]]
NumericVector ProdMatVecParallel(NumericMatrix mat, NumericVector vec, unsigned int chunkSize = 1) 
{
  // allocate the matrix we will return
  NumericVector rvec(mat.nrow());
    
  // create the worker
  ProdMatVec prodMatVec(mat, vec, rvec);
     
  // call it with parallelFor
  parallelFor(0, mat.nrow(), prodMatVec, chunkSize);

  return rvec;
}

// [[Rcpp::export()]]
List eigenPower_Rcpp_Parallel(const NumericMatrix A, const NumericVector v0,
  const double tol = 1e-6, const int maxit = 1e3,
  unsigned int chunkSize = 1,
  const int verbose = 0)
{
  // containers
  NumericVector v = v0;
  NumericVector b = v0;
        
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
    
    b = ProdMatVecParallel(A, v, chunkSize);
    v = b / innerNormParallel(b, chunkSize);
    lambda = innerProductParallel(v, b, chunkSize);
    
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
   
  return(ret);
}

//-----------------------------
// EigenPowerArma in Parallel
//-----------------------------

struct EigenPowerArma : public Worker 
{
  // input to read from
  const arma::mat& A;
  const arma::vec& v;
   
  // output to write to
  arma::vec& b;
  double b2_sum;
   
  unsigned int nchunks; 
   
  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  EigenPowerArma(const arma::mat& A, const arma::vec& v, arma::vec& b)
    : A(A), v(v), b(b), b2_sum(0), nchunks(0) {}
    
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) 
  {
    arma::mat A_i = A.rows(begin, end - 1);
        
    arma::vec b_i = A_i * v;
    
    // assign
    double b2_sum_i = 0;
    for(std::size_t k = begin; k < end; k++) {
      std::size_t ki = k - begin;

      double val = b_i[ki];
      
      b[k] = val;
      b2_sum_i += val*val;
    }
    b2_sum += b2_sum_i;
    
    ++nchunks;
  }
};

// [[Rcpp::export]]
List eigenPowerIt_Arma_Parallel(const arma::mat& A, const arma::vec& v, unsigned int chunkSize = 1) 
{
  // allocate the vector to be returned
  arma::vec b(A.n_rows);
    
  // create the worker
  EigenPowerArma eigenPowerArma(A, v, b);
     
  // call it with parallelFor
  parallelFor(0, A.n_rows, eigenPowerArma, chunkSize);

  double lambda = std::sqrt(eigenPowerArma.b2_sum);
  b = b / lambda;
  
  // return
  List ret;
    
  ret["lambda"] = lambda;
  ret["v"] = b;
   
  return(ret);
  
}

// [[Rcpp::export()]]
List eigenPower_Arma_Parallel(arma::mat& A, arma::vec& v0,
  const double tol = 1e-6, const int maxit = 1e3,
  unsigned int chunkSize = 1,
  const int verbose = 0)
{
  // containers
  arma::vec v = v0;
  arma::vec b = v0;
        
  // algorithm
  double lambda0 = 0;
  double lambda = lambda0;
  double delta = 0;
  bool converged = false;
  double nchunks = 0;  
  double csize = 0;
    
  int it = 1;
  int nchunks_sum = 0;
  for ( ; it <= maxit ; it++) { 
    if(verbose > 1) { 
      Rcout << " * it " << it << std::endl;
    }
    
    // compute in parallel
    EigenPowerArma eigenPowerArma(A, v, b);
    parallelFor(0, A.n_rows, eigenPowerArma, chunkSize);

    lambda = std::sqrt(eigenPowerArma.b2_sum);
    v = b * (1/lambda);
    
    nchunks_sum += eigenPowerArma.nchunks;
    
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
  nchunks = nchunks_sum / it;
  csize = (double)(A.n_rows) / nchunks;
  
  // return
  List ret;
  ret["v0"] = v0;
  ret["tol"] = tol;
  ret["maxit"] = maxit;
  ret["it"] = it;
  ret["delta"] = delta;
  ret["converged"] = converged;
  ret["nchunks"] = nchunks;
  ret["csize"] = csize;
    
  ret["lambda"] = lambda;
  ret["v"] = v;
   
  return(ret);
}
