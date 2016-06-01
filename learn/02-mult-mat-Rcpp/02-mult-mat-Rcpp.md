# Matrix multiplication using Rcpp
Andrey Ziyatdinov  
`r Sys.Date()`  







# About 

## References

* [Post](http://stackoverflow.com/questions/24933290/elementwise-matrix-multiplication-r-versus-rcpp-how-to-speed-this-code-up) that inspired this post.

# Include


```r
library(Rcpp)
library(RcppArmadillo)
```

# Functions


```r
sourceCpp(code ='
  #include <RcppArmadillo.h>
  // [[Rcpp::depends(RcppArmadillo)]]

  using namespace Rcpp;
  using namespace arma;

  // [[Rcpp::export]]
  arma::vec A_matvec_mult(const arma::mat & X, const arma::vec & y){
    arma::vec out = X * y;  
    return out;
  }'
)
```


```r
sourceCpp(code ='
  # include <RcppArmadillo.h>
  // [[Rcpp::depends(RcppArmadillo)]]

  using namespace Rcpp;
  
  // [[Rcpp::export()]]
  List eigenPower_Rcpp(const arma::mat A, arma::vec v0,
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
}')
```


```r
eigenPower_wrapper <- function(A, v0, ...)
{
  stopifnot(!missing(A))
  if(missing(v0)) {
    v0 <- runif(ncol(A))
  }
  
  eigenPower_Rcpp(A, v0, ...)
}
```


# A small example


```r
A <- matrix(c(-261, 209, -49, 
  -530, 422, -98,
  -800, 631, -144),
  ncol = 3, nrow = 3, byrow = TRUE)
  
v <- c(1, 2, 3)
```

## Matrix-vector multiplication



```r
A %*% v
```

```
     [,1]
[1,]   10
[2,]   20
[3,]   30
```

```r
A_matvec_mult(A, v)
```

```
     [,1]
[1,]   10
[2,]   20
[3,]   30
```

# Bigger examples


```r
n <- seq(250, 1500, by = 250)

out <- lapply(n, function(ni) {
  set.seed(1)
  M <-  matrix(runif(ni * ni), ni, ni)
  
  lt <- lower.tri(M)
  M[lt] <- t(M)[lt]
  
  list(n = ni,
    t.eigenPower = system.time(eigenPower(M))[["elapsed"]],
    t.eigenPower_Rcpp = system.time(eigenPower_wrapper(M))[["elapsed"]]
  )
})  
```


```r
df <- ldply(out, function(x) data.frame(n = x$n,
  t_eigenPower = x$t.eigenPower,
  t_eigenPower_Rcpp = x$t.eigenPower_Rcpp))
  
pf <- melt(df, id.vars = "n")

ggplot(pf, aes(n, value, color = variable)) + geom_point() + geom_line()
```

![](figures/b_ex_plot-1.png) 

## Direct call of eigenPower_Rcpp


```r
n <- seq(250, 1500, by = 250)

out <- lapply(n, function(ni) {
  set.seed(1)
  M <-  matrix(runif(ni * ni), ni, ni)
  
  lt <- lower.tri(M)
  M[lt] <- t(M)[lt]
  
  list(n = ni,
    t.eigenPower = system.time(eigenPower(M))[["elapsed"]],
    t.eigenPower_Rcpp_wrapper = system.time(eigenPower_wrapper(M))[["elapsed"]],    
    t.eigenPower_Rcpp = system.time(eigenPower_Rcpp(M, runif(ncol(M))))[["elapsed"]]
  )
})  
```


```r
df <- ldply(out, function(x) data.frame(n = x$n,
  t_eigenPower = x$t.eigenPower,
  t_eigenPower_Rcpp_wrapper = x$t.eigenPower_Rcpp_wrapper,
  t_eigenPower_Rcpp = x$t.eigenPower_Rcpp))
  
pf <- melt(df, id.vars = "n")

ggplot(pf, aes(n, value, color = variable)) + geom_point() + geom_line()
```

![](figures/b_ex_plot2-1.png) 
