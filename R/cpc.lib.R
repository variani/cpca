#-------------------------
# Main function `cpc`
#-------------------------

#' @export
cpc <- function(X, method = "power", k = 0, threshold = 0, ...)
{
  ### processing input argumets
  stopifnot(length(dim(X)) == 3)
  
  method <- match.arg(method, c("power"))
  
  switch(method,
    power = cpc_power(X, apply(X, 3, nrow), k, ...),
    stop("Error in swotch."))
}

cpc_power <- function(X, n_g, k = 0, iter = 15, ...)
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
  out <- list(D = D, CPC = CPC, ncomp = ncomp)
  return(out)
}
