#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

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
double innerProductParallel(NumericVector x, NumericVector y) 
{
  // declare the InnerProduct instance that takes a pointer to the vector data
  InnerProduct innerProduct(x, y);
   
  // call paralleReduce to start the work
  parallelReduce(0, x.length(), innerProduct);
   
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
double innerNormParallel(NumericVector x) 
{
  // declare the InnerProduct instance that takes a pointer to the vector data
  InnerNorm innerNorm(x);
   
  // call paralleReduce to start the work
  parallelReduce(0, x.length(), innerNorm);
   
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
NumericVector ProdMatVecParallel(NumericMatrix mat, NumericVector vec) 
{
  // allocate the matrix we will return
  NumericVector rvec(mat.nrow());
    
  // create the worker
  ProdMatVec prodMatVec(mat, vec, rvec);
     
  // call it with parallelFor
  parallelFor(0, mat.nrow(), prodMatVec);

  return rvec;
}
