#--------------------
# EigenPower class
#--------------------

#' S3 class EigenPower.
#'
#' @name EigenPowerClass
#' @rdname EigenPowerClass
#'
#' @param x 
#'    An object of class \code{EigenPower}.
#' @param digits
#'   The number of digits to be printed.
#' @param ...
#'    Additional arguments.
#'
#' @exportClass EigenPower

#' @rdname EigenPowerClass
#' @export
print.EigenPower <- function(x, digits = 2, ...)
{
  cat("\n Power EVD of class:", class(x)[1], "\n")
  if(x$converged) {
    cat(" - method converged in", x$it, "iterations ( maximum", x$maxit, ")\n")
  } else {
    cat(" - (!) method not converged in", x$it, "iterations ( maximum", x$maxit, ")\n")
  }  

  cat("\n Input data:\n")
  cat(" - matrix of size:", length(x$v), "rows\n")
  if(!is.null(x$sparseMatrix)) {
    if(x$sparseMatrix) {
      cat("  -- sparse matrix\n")
    }
  }
  
  cat("\n Output parameters:\n")
  cat(" - eigen value (lambda):", round(x$lambda, digits), "\n")
  cat(" - head of eigen vector (v):", paste(head(round(x$v, digits)), collapse = ", "), "\n")
  
  cat("\n CPU time:", modelParCPUtime(x, "auto"), "\n")
}

#------------------------
# EigenPowerRcpp class
#------------------------

#' S3 class EigenPowerRcpp.
#'
#' @name EigenPowerRcppClass
#' @rdname EigenPowerRcppClass
#'
#' @param x 
#'    An object of class \code{EigenPowerRcpp}.
#' @param ...
#'    Additional arguments.
#'
#' @exportClass EigenPowerRcpp

#' @rdname EigenPowerRcppClass
#' @export
print.EigenPowerRcpp <- function(x, ...)
{
  print.EigenPower(x, ...)
}

#------------------------
# modelPar* functions
#------------------------

modelParCPUtime <- function(x, format = "sec", ...) 
{ 
  t <- x$timing$cputime.sec
  
  if(format == "auto") {
    format <- ifelse(t < 60, "sec_str", "POSIX")
  }
  
  switch(format,
    "sec" = t,
    "sec_str" = paste(round(t, 3), "sec"),
    "POSIX" = format(.POSIXct(t, tz = "GMT"), "%H:%M:%S"),
    stop("swith error for `format`"))
}

