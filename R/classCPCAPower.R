#--------------------
# CPCAPower class
#--------------------

#' S3 class CPCAPower.
#'
#' @name CPCAPowerClass
#' @rdname CPCAPowerClass
#'
#' @param x 
#'    An object of class \code{CPCAPower}.
#' @param digits
#'   The number of digits to be printed.
#' @param ...
#'    Additional arguments.
#'
#' @exportClass CPCAPower

#' @rdname CPCAPowerClass
#' @export
print.CPCAPower <- function(x, digits = 2, ...)
{
  cat("\n Power CPCA of class:", class(x)[1], "\n")
  if(x$converged) {
    cat(" - method converged in (", paste(x$itComp, collapse = ", "), ") iterations ( maximum iterations", x$maxit, ")\n")
  } else {
    cat(" - (!) method not converged for some/all components ( maximum iterations", x$maxit, ")\n")
  }  

  cat("\n Input data:\n")
  cat(" - covariance matrices of size:", nrow(x$CPC), "\n")
  
  cat("\n Output parameters:\n")
  #cat(" - eigen value (lambda):", round(x$lambda, digits), "\n")
  #cat(" - head of eigen vector (v):", paste(head(round(x$v, digits)), collapse = ", "), "\n")
  
  #cat("\n CPU time:", modelParCPUtime(x, "auto"), "\n")
}
