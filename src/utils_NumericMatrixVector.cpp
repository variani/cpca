# include <Rcpp.h>

using namespace Rcpp;

/*
  Functions: 
   - (1) input arguments are NumericMatrix/NumericVector;
   - (2) make use of std::inner_product
*/

/*
// [[Rcpp::export()]]
inline double stdNumericInnerProduct(const NumericVector & x, const NumericVector & y)
{
  return std::inner_product(x.begin(), x.end(), y.begin(), 0);
}

// [[Rcpp::export()]]
inline double stdNumericNorm(const NumericVector & x)
{
  return std::sqrt(std::inner_product(x.begin(), x.end(), x.begin(), 0));
}
*/

/*
  Functions: 
   - (1) input arguments are NumericMatrix/NumericVector;
   - (2) make use of `for` loops
*/


// [[Rcpp::export()]]
double numericNorm(const NumericVector & x)
{
  // looping variables
  unsigned int n = x.length();
  
  // sum
  double res = 0;
  for(unsigned int i = 0; i < n; i++) {
    double val = x(i);

    res += val * val;
  }
  res = sqrt(res);
  
  return res;
}

/*
// [[Rcpp::export()]]
double numericInnerProduct(const NumericVector & x, const NumericVector & y)
{
  // looping variables
  unsigned int n = x.length();
  
  // sum
  double res = 0;
  for(unsigned int i = 0; i < n; i++) {
    res += x(i) * y(i);
  }
  res = sqrt(res);
  
  return res;
}
*/

// [[Rcpp::export()]]
NumericVector numericProdMatVec(const NumericMatrix & mat, const NumericVector & vec)
{
  // allocate the output
  NumericVector rvec(vec.length());
  
  // looping variables
  unsigned int n = mat.nrow(), k = mat.ncol();
  
  // fill `rvec`
  for(unsigned int i = 0; i < n; i++) {
    double res = 0;
    for(unsigned int j = 0; j < k; j++) {
      res += mat(i, j) * vec(j);
    }
    rvec[i] = res;
  }
  
  // return the output
  return rvec;
}


/*
  Other functions: 
   - (1) input arguments are NumericMatrix/NumericVector;
*/

// [[Rcpp::export()]]
NumericVector numericMultVec(const NumericVector & x, double a)
{
  // allocate the output
  NumericVector y(x.length());

  for(int i = 0; i < x.length(); ++i) {
    y(i) = x(i) * a;
  }
  return y;
}

void numericUpdateMultVec(NumericVector & x, double a)
{
  for(int i = 0; i < x.length(); ++i) {
    x(i) = x(i) * a;
  }
}
