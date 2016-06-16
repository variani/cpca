#----------------------------------
# Data simulations
#----------------------------------

#' Function hilbert. 
#'
#' This function is used to produce a symmetric matrix.
#'
#' @param n The matrix size.
#' @return A symmetric matrix
#' @references \url{https://stat.ethz.ch/R-manual/R-devel/library/base/html/svd.html}
#' @export
hilbert <- function(n) 
{ 
  i <- 1:n
  1 / outer(i - 1, i, "+") 
}
