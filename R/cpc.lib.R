#-------------------------
# Main function `cpc`
#-------------------------

#' Function cpc. 
#'
#' This function computes the CPCA from a given set of covariance matrices 
#' (of different groups). 
#'
#' Currently, the only the stepwise algorithm by Trendafilov is supported.
#'
#' @name cpc
#' @param X An array of three dimensions: the 3rd dimension encodes the groups
#'   and the first two dimension contain the covariance matrices.
#' @param method The name of the method for computing the CPCA.
#'   The default value is \code{"stepwise"}, which is the stepwise algorithm by Trendafilov.
#' @param k The number of components to be computed (all if it is \code{0}).
#'   This paramter is valid if the given method supports 
#'   built-in ordering of the eigvenvectors.
#'   The default value is \code{0}, that means computing of all the components.
#' @param iter The maximum number of iterations.
#'   The parameter is valid for the stepwise algorithm by Trendafilov,
#'   that is applied in the power algorithm for estimation a single component.
#'   The default value is 30.
#' @param threshold The threshold value of the captured variance,
#'   which is reserved for further extensions.
#' @param ... Other parameters.
#' @return A list several slots: \code{CPC} rotation matrix with eigenvectors in columns;
#'   \code{ncomp} the number of components evaluated (equal to the number of columns in \code{CPC}).
#' @note This function adpats the original code in matlab written by Dr N. T. Trendafilov.
#' @references Trendafilov (2010). Stepwise estimation of common principal components. 
#'   Computational Statistics & Data Analysis, 54(12), 3446-3457. 
#'   doi:10.1016/j.csda.2010.03.010
#' @example inst/examples/function-cpc.R
#' @export
cpc <- function(X, method = "stepwise", k = 0, iter = 30, threshold = 0, ...)
{
  ### processing input argumets
  stopifnot(length(dim(X)) == 3)
  
  method <- match.arg(method, c("stepwise"))
  
  switch(method,
    stepwise = cpc_stepwise(X, apply(X, 3, nrow), k, iter, ...),
    stop("Error in swotch."))
}

#-------------------------------------------------------------------------------
# `cpc_stepwise` function was the first version of the power algorithm for CPCA.
# The code was adopted from the Matlab version written by Dr. Trendafilov.
# 
# In the versions of `cpca` earlier than 0.2.0, this code was the only 
# implementation available. The missing parts of the algorithms iclude:
#  - the cost function is missing;
#  - as a consequence, the number of iterations is fixed by the user;
#  - `%*%` operator is used instead of a more efficient ones, i.e. `crossprod`.
#-------------------------------------------------------------------------------

cpc_stepwise <- function(X, n_g, k = 0, iter = 30, ...)
{
  p <- dim(X)[1]
  mcas <- dim(X)[3]

  # If k = 0 retrieve all components
  if(k == 0) {
    k <- p
  }
  
  # parameters: number of interations
  iter <- 15
  n <- n_g / sum(n_g)
  
  # output variables
  D <- array(0, dim=c(p, mcas))
  CPC <- array(0, dim=c(p, p))
  Qw <- diag(1, p)
  
  # components `s`
  s <- array(0, dim = c(p, p))
  for(m in 1:mcas) {
    s <- s + n[m]*X[, , m]
  }
  
  # variables
  res <- eigen(s)
  q0 <- res$vectors
  d0 <- diag(res$values, p)
  if(d0[1, 1] < d0[p, p]) {
    q0 <- q0[, ncol(q0):1]
  }
  
  # loop 'for ncomp=1:p'
  # Replaced by k so that only the first k components are retrieved
  for(ncomp in 1:k) {
   q <- q0[, ncomp]
   d <- array(0, dim=c(1, mcas))
   for(m in 1:mcas) {
    d[, m] <- t(q) %*% X[, , m] %*% q
   }
   
   # loop 'for i=1:iter'
   for(i in 1:iter) {
    s <- array(0, dim=c(p, p))
    for(m in 1:mcas) {
      s <- s + n_g[m] * X[, , m] / d[, m]
    }
           
    w <- s %*% q
    if( ncomp != 1) {
      w <- Qw %*% w
    }

    q <- w / as.numeric(sqrt((t(w) %*% w)))
    for(m in 1:mcas) {
      d[, m]  <- t(q) %*% X[, , m] %*% q
    }
    
   }
   # end of loop 'for i=1:iter'
   
   D[ncomp, ] <- d
   CPC[, ncomp] <- q
   Qw <- Qw - q %*% t(q)
  }
  # end of loop 'for ncomp=1:k'
  
  ### return
  out <- list(D = D[1:ncomp, ], CPC = CPC[, 1:ncomp], ncomp = ncomp)
  return(out)
}

#-------------------------------------------------------------------------------
# `cpca_stepwise_base` is an updated version of `cpc_stepwise`.
# - the cost function is introduced;
# - new arguments `maxit`, `tol` are added;
# - the input covariance matrices are passed in a list.
# 
# Note:
# - `%*%` operator is still used.
# 
# This function is mainly needed for comparison with more efficient
# implementations, e.g. `cpc_stepwise_eigen`.
#-------------------------------------------------------------------------------

cpca_stepwise_base <- function(cov, ng, ncomp = 0, 
  tol = 1e-6, maxit = 1e3,
  start = c("eigenPower", "eigen", "random"), symmetric = TRUE,
  useCrossprod = TRUE,
  verbose = 0, ...)
{
  ### par
  stopifnot(length(cov) == length(ng))
  stopifnot(all(laply(cov, class) == "matrix"))
  
  stopifnot(symmetric)
  
  start <- match.arg(start)
    
  ng <- as.numeric(ng) # needed for further multiplication with double precision like `ng / n` 
  
  ### var
  k <- length(cov) # number of groups `k`
  p <- nrow(cov[[1]]) # number of variables `p` (the covariance matrices p x p)
  n <- sum(ng) # the number of samples
  
  # If ncomp = 0 retrieve all components
  if(ncomp == 0) {
    ncomp <- p
  }

  ### output variables
  D <- matrix(0, nrow = ncomp, ncol = k)
  CPC <- matrix(0, nrow = p, ncol = ncomp)
  Qw <- diag(1, p)
  
  convergedComp <- rep(FALSE, ncomp)
  itComp <- rep(0, ncomp)
  
  ### step 1: compute the staring estimation
  if(start == "eigenPower") {
    S <- matrix(0, nrow = p, ncol = p)  
    for(i in 1:k) {
        S <- S + (ng[i] / n) * cov[[i]]
    }
    
    res <- eigenPower(S, ncomp = ncomp)
    all(order(res$values[1:ncomp], decreasing = TRUE) == seq(1, ncomp))
      
    q0 <- res$vectors[, 1:ncomp, drop = FALSE]
  }
  else if(start == "eigen") {
    S <- matrix(0, nrow = p, ncol = p)  
    for(i in 1:k) {
        S <- S + (ng[i] / n) * cov[[i]]
    }
    
    res <- eigen(S, symmetric = symmetric) # `?eigen`: a vector containing the p eigenvalues of ‘x’, sorted in _decreasing_ order
    all(order(res$values[1:ncomp], decreasing = TRUE) == seq(1, ncomp))
      
    q0 <- res$vectors[, 1:ncomp, drop = FALSE]
  } else {
    q0 <- matrix(runif(p * ncomp), nrow = p, ncol = ncomp)  
  }
  
  #### step 2: estimation of `p` components in a loop
  for(comp in 1:ncomp) {
    if(verbose > 1) {
      cat(" * component:", comp, "/", ncomp, "\n")
    }
    
    q <- q0[, comp]
  
    d <- rep(0, k) 
    for(i in 1:k) {
      if(useCrossprod) d[i] <- as.numeric(crossprod(q, crossprod(cov[[i]], q)))
      else             d[i] <- as.numeric(t(q) %*% cov[[i]] %*% q)
    }
   
    # loop along `it`
    cost0 <- 0
    for(it in 1:maxit) {
      if(verbose > 1) {
        cat(" * it:", it, "/", maxit, "\n")
      }
    
      S <- matrix(0, nrow = p, ncol = p)
      for(i in 1:k) {
        S <- S + (ng[i] / d[i]) * cov[[i]]
      }
      
      if(useCrossprod) w <- crossprod(S, q) # (1) S %*% q; (2) tcrossprod(S, t(q)), as S is symmetric
      else             w <- S %*% q

      if(comp != 1) { 
        if(useCrossprod) w <- Qw %*% w
        else             w <- crossprod(Qw, w) # (1) Qw %*% w; (2) tcrossprod(Qw, t(w)), as Qw is symmetric
      }

      if(useCrossprod) q <- w / sqrt(as.numeric(t(w) %*% w)) # normalize 
      else             q <- as.numeric(w) / sqrt(as.numeric(crossprod(w))) # normalize 
      
      # compute `cost` & `d`
      for(i in 1:k) {
        if(useCrossprod) d[i] <- as.numeric(crossprod(q, crossprod(cov[[i]], q)))
        else             d[i] <- as.numeric(t(q) %*% cov[[i]] %*% q)
      }
      
      cost <- sum(log(d) * ng)
      
      delta <- abs((cost - cost0) / cost)
      if(verbose > 1) {
        cat("  --  delta:", delta, "\n")
      }
      if(delta < tol) {
        break
      }
    
      cost0 <- cost
    }
    # end of loop along `it`
    itComp[comp] <- it
    convergedComp[comp] <- (it < maxit)
   
    D[comp, ] <- d
    CPC[, comp] <- q

    if(useCrossprod) Qw <- Qw - q %*% t(q)
    else             Qw <- Qw - tcrossprod(q)
    
  }
  # end of loop along `comp`
  
  ### return
  out <- list(D = D, CPC = CPC, ncomp = ncomp,
    convergedComp = convergedComp, converged = all(convergedComp),
    itComp = itComp, maxit = maxit)
  
  oldClass(out) <- c("CPCAPowerBase", "CPCAPower")
  
  return(out)
}

#-------------------------------------------------------------------------------
# `cpca_stepwise_eigen` is an updated version of `cpc_stepwise_base`.
# - the input covariance matrices are of one of `Matrix` types;
# - the property, that the covariance matrices are symmetric, is used for speed up.
# 
# Note:
# - `%*%` operator is replaced by more efficient ones, e.g. `crossprod`;
#-------------------------------------------------------------------------------

#' @importFrom Matrix Matrix isSymmetric crossprod tcrossprod
cpca_stepwise_eigen <- function(cov, ng, ncomp = 0, 
  tol = 1e-6, maxit = 1e3,
  start = c("eigen", "random"), symmetric = TRUE,
  verbose = 0, ...)
{
  ### par
  stopifnot(length(cov) == length(ng))
  
  cl <- laply(cov, function(x) as.character(class(x)))
  stopifnot(all(grepl("*Matrix$", cl)))

  stopifnot(symmetric)
  covSymm <- laply(cov, Matrix::isSymmetric)
  stopifnot(all(covSymm))
  
  start <- match.arg(start)
  
  ng <- as.numeric(ng) # needed for further multiplication with double precision like `ng / n` 
  
  ### var
  k <- length(cov) # number of groups `k`
  p <- nrow(cov[[1]]) # number of variables `p` (the covariance matrices p x p)
  n <- sum(ng) # the number of samples
  
  # If ncomp = 0 retrieve all components
  if(ncomp == 0) {
    ncomp <- p
  }

  ### output variables
  D <- matrix(0, nrow = ncomp, ncol = k)
  CPC <- matrix(0, nrow = p, ncol = ncomp)
  
  Qw <- Matrix::Matrix(diag(1, p))
  
  convergedComp <- rep(FALSE, ncomp)
  itComp <- rep(0, ncomp)
  
  ### step 1: compute the staring estimation
  if(start == "eigen") {
    S <- Matrix::Matrix(0, nrow = p, ncol = p)
    stopifnot(Matrix::isSymmetric(S))

    for(i in 1:k) {
      S <- S + (ng[i] / n) * cov[[i]]
    }
    
    res <- eigen(S, symmetric = symmetric) # `?eigen`: a vector containing the p eigenvalues of ‘x’, sorted in _decreasing_ order
    all(order(res$values[1:ncomp], decreasing = TRUE) == seq(1, ncomp))
      
    q0 <- Matrix::Matrix(res$vectors[, 1:ncomp, drop = FALSE])
  } else {
    q0 <- Matrix::Matrix(runif(p * ncomp), nrow = p, ncol = ncomp)
  }
  
  #### step 2: estimation of `p` components in a loop
  for(comp in 1:ncomp) {
    if(verbose > 1) {
      cat(" * component:", comp, "/", ncomp, "\n")
    }
    
    q <- q0[, comp]
  
    d <- rep(0, k) 
    for(i in 1:k) {
      d[i] <- as.numeric(Matrix::crossprod(q, Matrix::crossprod(cov[[i]], q)))
    }
   
    # loop along `it`
    cost0 <- 0
    for(it in 1:maxit) {
      if(verbose > 1) {
        cat(" * it:", it, "/", maxit, "\n")
      }
    
      S <- Matrix::Matrix(0, nrow = p, ncol = p)
      for(i in 1:k) {
        S <- S + (ng[i] / d[i]) * cov[[i]]
      }
      
      w <- Matrix::crossprod(S, q) # (1) S %*% q; (2) tcrossprod(S, t(q)), as S is symmetric
      if(comp != 1) { 
        w <- Matrix::crossprod(Qw, w) # (1) Qw %*% w; (2) tcrossprod(Qw, t(w)), as Qw is symmetric
      }

      q <- as.numeric(w) / sqrt(as.numeric(Matrix::crossprod(w))) # normalize 
      
      # compute `cost` & `d`
      for(i in 1:k) {
        d[i] <- as.numeric(Matrix::crossprod(q, Matrix::crossprod(cov[[i]], q)))
      }
      
      cost <- sum(log(d) * ng)
      
      delta <- abs((cost - cost0) / cost)
      if(verbose > 1) {
        cat("  --  delta:", delta, "\n")
      }
      if(delta < tol) {
        break
      }
    
      cost0 <- cost
    }
    # end of loop along `it`
    itComp[comp] <- it
    convergedComp[comp] <- (it < maxit)
    
    D[comp, ] <- d
    CPC[, comp] <- q
    
    Qw <- Qw - Matrix::tcrossprod(q)
  }
  # end of loop along `comp`
  
  ### return
  out <- list(D = D, CPC = CPC, ncomp = ncomp,
    convergedComp = convergedComp, converged = all(convergedComp),
    itComp = itComp, maxit = maxit)
  
  oldClass(out) <- c("CPCAPowerEigen", "CPCAPower")
  
  return(out)
}
