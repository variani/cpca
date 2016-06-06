#-------------------------
# EVD functions
#-------------------------

#' Function eigenPower. 
#'
#' The function implements the power algorithm for EVD.
#'
#' @name eigenPower
#' @param A A two-dimensional square matrix, either of \code{matrix} or \code{Matrix} class.
#' @param v0 A numeric vector; the initial guess for eignevector.
#'   If it is missing, a random vector is generated.
#' @param tol The tolerance threshold used to stop when reaching no improvement if estmiation of eigenvalue.
#'   The default value is \code{1e-6}.
#' @param maxit The maximum number of iterations.
#'   The default value is \code{1e4}.
#' @param sparse The boolean value, whether to convert the input \code{A} matrix to one of the \code{Matrix} classes.
#'   The default value is \code{FALSE}.
#' @param sparseSymm The boolean value, whether to convert the input \code{A} matrix to one of the \code{Matrix} classes,
#'   while trying to convert to a symmetric matrix type.
#'   The default value is \code{FALSE}.
#' @param verbose The integer value indicating the verbose level.
#'   The default value is \code{0}.
#' @return A list several slots: \code{v} the first eigenvector; 
#'   \code{lambda} the first eigenvalue; etc.
#' @note This function is inspired by the post \url{http://blogs.sas.com/content/iml/2012/05/09/the-power-method.html}.
#' @example inst/examples/function-eigenPower.R
#' @importFrom Matrix t crossprod tcrossprod
#' @export
eigenPower <- function(A, v0, tol = 1e-6, maxit = 1e3, 
  sparse = FALSE, sparseSymm = FALSE, 
  verbose = 0)
{
  timing <- list()
  timing$args <- proc.time()
  
  ### arguments
  stopifnot(!missing(A))
  
  if(missing(v0)) {
    v0 <- runif(ncol(A))
  }
  
  ### convert into sparse matrices
  if(any(sparse, sparseSymm)) {
    stopifnot(requireNamespace("Matrix"))
    
    if(sparse) {
      A <- Matrix::Matrix(A, sparse = T)
    }
    
    if(sparseSymm) {
      A <- Matrix::Matrix(A, sparse = T)
      A <- as(A, "symmetricMatrix")
    }
  }
  
  ### vars
  sparseMatrix <- FALSE
  cl <- attr(class(A), "package")
  if(!is.null(cl)) {
    if(cl == "Matrix") {
      sparseMatrix <- TRUE
      
      stopifnot(requireNamespace("Matrix"))
    }
  }
  
  ### preparation before looping
  timing$algo <- proc.time()
  
  v0 <- as.numeric(v0)
  
  v <- v0
  v <- v / sqrt(v %*% v)

  ### loop
  lambda0 <- 0
  for(it in 1:maxit) {
    if(verbose > 1) {
      cat(" * it:", it, "/", maxit, "\n")
    }
    
    b <- tcrossprod(A, t(v)) # A %*% v
    v <- b / sqrt(as.numeric(crossprod(b)))
    lambda <- as.numeric(crossprod(v, b)) # t(v) %*% b
    
    delta <- abs((lambda - lambda0) / lambda)
    if(delta < tol) {
      break
    }
    
    lambda0 <- lambda
  }
  
  ### post-process
  converged <- (it < maxit)
  
  ### output
  timing$return <- proc.time()
  
  timing$targs <- timing$algo[["elapsed"]] - timing$args[["elapsed"]]
  timing$talgo <- timing$return[["elapsed"]] - timing$algo[["elapsed"]]
  
  out <- list(v0 = v0, tol = tol, maxit = maxit,
    it = it, delta = delta, converged = converged, 
    sparseMatrix = sparseMatrix,
    timing = timing,
    lambda = lambda, v = v)
  
  return(out)
}

#' Function eigenPowerRcpp. 
#'
#' The function implements the power algorithm for EVD using Rcpp.
#'
#' @name eigenPowerRcpp
#' @param A A two-dimensional square matrix, either of \code{matrix} or \code{Matrix} class.
#' @param v0 A numeric vector; the initial guess for eignevector.
#'   If it is missing, a random vector is generated.
#' @param tol The tolerance threshold used to stop when reaching no improvement if estmiation of eigenvalue.
#'   The default value is \code{1e-6}.
#' @param maxit The maximum number of iterations.
#'   The default value is \code{1e4}.
#' @param verbose The integer value indicating the verbose level.
#'   The default value is \code{0}.
#' @return A list several slots: \code{v} the first eigenvector; 
#'   \code{lambda} the first eigenvalue; etc.
#' @export
eigenPowerRcpp <- function(A, v0, tol = 1e-6, maxit = 1e3, 
  verbose = 0)
{
  stopifnot(!missing(A))
 
  if(missing(v0)) {
    v0 <- runif(ncol(A))
  }
  
  eigenPower_Rcpp(A, v0, tol = tol, maxit = maxit, verbose = verbose)
}
