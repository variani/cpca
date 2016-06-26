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
library(png)

library(ggplot2)
library(gridExtra)

### par
imgfile <- "data/lena_std.png"

### fucntions copied from ref #1
compress_svd <- function(m, k) 
{
  factorize <- function(m, c, k)
  {
    cat(" SVD on", c, "...\n")
    fac <- svd(m[, , c])

    m <- fac$u[, 1:k] %*% diag(fac$d[1:k]) %*% t(fac$v[, 1:k])
    
    m[m < 0] <- 0
    m[m > 1] <- 1
        
    return(m)        
  }
    
  rimg <- array(c(factorize(m, 1, k), factorize(m, 2, k), factorize(m, 3, k)), 
    dim(img))
  return(rimg)
}

compress_eigen <- function(m, k, center = TRUE) 
{
  factorize <- function(m, c, k)
  {
    cat(" EVD on", c, "...\n")
    mat <- m[, , c]
    
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
    
  rimg <- array(c(factorize(m, 1, k), factorize(m, 2, k), factorize(m, 3, k)), 
    dim(img))
  return(rimg)
}


viewPNG <- function(img, tmpfile = "tmp.png") 
{
  writePNG(img, tmpfile)
  system("eog tmp.png")
}  


### read file
img <- readPNG(imgfile)

### create a figure
p <- marrangeGrob(rasterGrob(compress_eigen(img, 5)), 
  rasterGrob(compress_eigen(img, 25)), 
  rasterGrob(img), 
  nrow = 1, ncol = 3, 
  #top = "Compressed image from lenna.org, using EVD/SVD with k = {5, 25, all} components")
  top = NULL)
  
ggsave("lena.pdf", p, width = 12, height = 4)



