# References:
# @ https://github.com/road2stat/imgsvd/edit/master/server.R
# @ https://en.wikipedia.org/wiki/Principal_component_analysis
# @ (ggplot + png) https://dgkontopoulos.wordpress.com/2012/11/22/setting-a-background-image-in-ggplot2/
#
# Notes:
# - Center data using EVD? Try `compress_eigen` varying `center` argument.
#   On `lena_std.png` image and with k = 5, one can see the effect of centering.
#   In particular, the centered approach gives a more contrast image.

### inc 
library(plyr)

library(png)

library(ggplot2)
library(gridExtra)

compress_pca <- function(mat, k, center = TRUE, transpose = FALSE) 
{
  factorize <- function(mat, k, center)
  {
    if(center) {
      means <- colMeans(mat)
      mat <- sweep(mat, 2, means, "-") 
    }
      
    fac <- eigen(crossprod(mat), symmetric = TRUE)
    
    m <- mat %*% fac$vectors[, 1:k] %*% t(fac$vectors[, 1:k])
    
    if(center) {
      m <- sweep(m, 2, means, "+") 
    }
          
    m[m < 0] <- 0
    m[m > 1] <- 1
        
    return(m)        
  }

  dat <- matrix(as.numeric(NA), nrow(mat), ncol(mat))
  for(i in 1:nrow(mat)) {
    cat(" * ", i, "/", nrow(mat), "\n")
    
    mati <- matrix(mat[i, ], 64, 64)
    if(transpose) {
      mati <- t(mati)
    }

    m <- factorize(mati, k = k, center = center)

    if(transpose) {
      m <- t(m)
    }
    
    dat[i, ] <- as.vector(m)
  }
  
  return(dat)
}

compress_cpca <- function(mat, k, center = TRUE, transpose = FALSE) 
{
  ### prepare data
  means <- list()
  cov <- list()
  ng <- rep(as.integer(NA), nrow(mat))
  
  for(i in 1:nrow(mat)) {
    cat(" * ", i, "/", nrow(mat), "\n")
    
    mati <- matrix(mat[i, ], 64, 64)
    if(transpose) {
      mati <- t(mati)
    }
    
    meansi <- colMeans(mati)
    if(center) {
      mati <- sweep(mati, 2, meansi, "-") 
    }
    
    means[[i]] <- meansi
    cov[[i]] <- cov(mati)
    ng[i] <- nrow(mati)
  }
  
  ### run CPCA  
  out <- cpca(cov, ng, ncomp = k)
    
  ### reconstruct
  dat <- matrix(as.numeric(NA), nrow(mat), ncol(mat))
  for(i in 1:nrow(mat)) {
    cat(" * ", i, "/", nrow(mat), "\n")
    
    mati <- matrix(mat[i, ], 64, 64)
    
    if(transpose) {
      mati <- t(mati)
    }
    
    meansi <- means[[i]]
    
    if(center) {
      mati <- sweep(mati, 2, meansi, "-") 
    }
    
    m <- mati %*% out$CPC[, 1:k, drop = FALSE] %*% t(out$CPC[, 1:k, drop = FALSE])
    
    if(center) {
      m <- sweep(m, 2, means[[i]], "+") 
    }

    if(transpose) {
      m <- t(m)
    }
              
    m[m < 0] <- 0
    m[m > 1] <- 1
        
    dat[i, ] <- as.vector(m)
  }
  
  return(dat)

}

#plot(1:64, 1:64, type = 'n')
#rasterImage(img, 1, 1, 64, 64)

plotFaces <- function(dat, ncol = 10) 
{
  grobs <- llply(1:nrow(dat), function(i) {
    rasterGrob(matrix(dat[i, ], 64, 64))
  })
  p <- marrangeGrob(grobs, nrow = ceiling(nrow(dat) / ncol), ncol = ncol, top = NULL)
  
  return(p)
}

### data
data(faces, package = "RnavGraphImageData")

faces <- faces / 255 # 0-255 values, i.e. 8-bit encoding

faces <- t(faces)

viewFaces <- function(dat, k, center = TRUE, transpose = FALSE, ncol = 10) 
{
  plotFaces(rbind(dat,
    compress_pca(dat, k, center = center, transpose = transpose),
    compress_cpca(dat, k, center = center, transpose = transpose)), ncol = ncol)
}

