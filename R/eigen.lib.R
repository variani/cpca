#-------------------------
# EVD functions
#-------------------------

#' Function eigenPower. 
#'
#' The function implements the power algorithm for EVD.
#'
#' @name eigenPower
#' @rdname eigenPower
#'
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
#' @param ncomp The number of eigenvectors to be extracted.
#'   The default value is \code{1}.
#'   The value of \code{0} means extract all eigenvectors.
#' @param symmetric A logical which value says explicetly if the input matrix \code{A} is symmetric.
#'    The default value is \code{FALSE}.
#' @param cores The number of cores (for parallel versions).
#'   The default value is \code{-1}.
#'   This argument is passed next to \code{RcppParallel::setThreadOptions(numThreads = cores)}.
#' @param chunkSize The minimal size of a chunk (for parallel versions).
#'   The default value is \code{1}.
#'   This argument is passed next to a wrapper \code{eigenPower_Rcpp_Parallel}.
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
  ncomp = 1,
  verbose = 0)
{
  timing <- list()
  timing$args <- proc.time()
  
  ### arguments
  stopifnot(!missing(A))
  stopifnot(nrow(A) == ncol(A))
  
  p <- ncol(A)

  if(ncomp != 1) {
    stopifnot(missing(v0))
  }
  
  if(missing(v0)) {
    v0 <- runif(p)
  }
  v0 <- as.numeric(v0)
  
  if(ncomp == 0) {
    ncomp <- p
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

  # output variables
  convergedComp <- rep(FALSE, ncomp)
  itComp <- rep(0, ncomp)
  
  D <- rep(0, ncomp)
  CPC <- matrix(0, nrow = p, ncol = ncomp)
  Qw <- diag(1, p)
  if(sparseMatrix) {
    Qw <- Matrix::Matrix(Qw, sparse = sparseMatrix)
  }
   
  ### computation
  timing$algo <- proc.time()

  for(comp in 1:ncomp) {
    if(verbose > 1) {
      cat(" * component:", comp, "/", ncomp, "\n")
    }
    
    ### preparation before looping
    v <- v0
    v <- v / sqrt(as.numeric(crossprod(v)))

    ### loop
    lambda0 <- 0
    for(it in 1:maxit) {
      if(verbose > 1) {
        cat(" * it:", it, "/", maxit, "\n")
      }
    
      b <- tcrossprod(A, t(v)) # A %*% v
      if(comp != 1) { 
        b <- crossprod(Qw, b)
      }
      v <- b / sqrt(as.numeric(crossprod(b)))
      lambda <- as.numeric(crossprod(v, b)) # t(v) %*% b
    
      delta <- abs((lambda - lambda0) / lambda)
      if(delta < tol) {
        break
      }
    
      lambda0 <- lambda
    }
  
    ### post-process
    itComp[comp] <- it
    convergedComp[comp] <- (it < maxit)
    
    D[comp] <- lambda
    CPC[, comp] <- as.numeric(v)
     
    Qw <- Qw - tcrossprod(v)
  }
  
  ### output
  timing$return <- proc.time()
  
  timing$targs <- timing$algo[["elapsed"]] - timing$args[["elapsed"]]
  timing$talgo <- timing$return[["elapsed"]] - timing$algo[["elapsed"]]

  timing$cputime.sec <- (timing$return - timing$args)[["elapsed"]]
  
  out <- list(v0 = v0, tol = tol, maxit = maxit,
    itComp = itComp, convergedComp = convergedComp,
    it = mean(it), converged = all(convergedComp), 
    sparseMatrix = sparseMatrix,
    timing = timing,
    ncomp = ncomp,
    values = D, vectors = CPC,
    lambda = D[1], v = CPC[, 1])
    
  oldClass(out) <- c("EigenPower")
   
  return(out)
}

#' Function eigenPowerRcpp. 
#'
#' @rdname eigenPower
#'
#' @export
eigenPowerRcpp <- function(A, v0, tol = 1e-6, maxit = 1e3, mode = 1,
  verbose = 0)
{
  ### args
  timing <- list()
  timing$args <- proc.time()
  
  stopifnot(!missing(A))
 
  if(missing(v0)) {
    v0 <- runif(ncol(A))
  }
  
  ### run
  out <- eigenPower_Rcpp(A, v0, tol = tol, maxit = maxit, mode = mode, verbose = verbose)

  ### return
  timing$return <- proc.time()
  
  timing$cputime.sec <- (timing$return - timing$args)[["elapsed"]]
  
  
  out$timing <- timing
  oldClass(out) <- c("EigenPowerRcpp", "EigenPower")
  
  return(out)
}

#' Function eigenPowerRcppEigen. 
#'
#' @rdname eigenPower
#'
#' @export
eigenPowerRcppEigen <- function(A, v0, tol = 1e-6, maxit = 1e3, 
  ncomp = 1, symmetric = FALSE,
  verbose = 0)
{
  ### args
  timing <- list()
  timing$args <- proc.time()
  
  stopifnot(!missing(A))
 
  if(missing(v0)) {
    v0 <- runif(ncol(A))
  }
  
  ### run
  out <- eigenPower_RcppEigen(A, v0, tol = tol, maxit = maxit,   ncomp = ncomp, symmetric = symmetric, verbose = verbose)

  ### return
  timing$return <- proc.time()
  
  timing$cputime.sec <- (timing$return - timing$args)[["elapsed"]]
  
  out$timing <- timing
  oldClass(out) <- c("EigenPowerRcppEigen", "EigenPower")
  
  return(out)
}

#-------------------------
# EVD parallelized functions
#-------------------------

#' Function eigenPowerRcppParallel. 
#'
#' @rdname eigenPower
#'
#' @export
eigenPowerRcppParallel <- function(A, v0, tol = 1e-6, maxit = 1e3, 
  cores = -1, chunkSize = 1,
  verbose = 0)
{
  ### args
  timing <- list()
  timing$args <- proc.time()
 
  stopifnot(!missing(A))
 
  if(missing(v0)) {
    v0 <- runif(ncol(A))
  }
  
  if(!missing(cores)) { 
    RcppParallel::setThreadOptions(numThreads = cores)
  }
  
  ### run
  stop()
  #out <- eigenPower_Rcpp_Parallel(A, v0, tol = tol, maxit = maxit, verbose = verbose)
  
  out$v <- as.numeric(out$v)
  
  ### return
  timing$return <- proc.time()
  
  timing$cputime.sec <- (timing$return - timing$args)[["elapsed"]]
  
  
  out$timing <- timing
  oldClass(out) <- c("EigenPowerRcppParallel", "EigenPower")
  
  return(out)
}

#' Function eigenPowerEigenParallel. 
#'
#' @rdname eigenPower
#'
#' @export
eigenPowerEigenParallel <- function(A, v0, tol = 1e-6, maxit = 1e3, 
  cores = -1, chunkSize = 1,
  verbose = 0)
{
  ### args
  timing <- list()
  timing$args <- proc.time()
 
  stopifnot(!missing(A))
 
  if(missing(v0)) {
    v0 <- runif(ncol(A))
  }
  
  if(!missing(cores)) { 
    RcppParallel::setThreadOptions(numThreads = cores)
  }
  
  ### run
  stop()
  #out <- eigenPower_RcppEigen_Parallel(A, v0, tol = tol, maxit = maxit, chunkSize = chunkSize, verbose = verbose)
  
  out$v <- as.numeric(out$v)
  
  ### return
  timing$return <- proc.time()
  
  timing$cputime.sec <- (timing$return - timing$args)[["elapsed"]]
  
  
  out$timing <- timing
  oldClass(out) <- c("EigenPowerEigenParallel", "EigenPower")
  
  return(out)
}
