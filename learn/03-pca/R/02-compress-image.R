# @ https://github.com/road2stat/imgsvd/edit/master/server.R
# @ https://en.wikipedia.org/wiki/Principal_component_analysis

### inc 
library(png)

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

compress_eigen <- function(m, k) 
{
  factorize <- function(m, c, k)
  {
    cat(" EVD on", c, "...\n")
    mat <- m[, , c]
    
    fac <- eigen(crossprod(mat), symmetric = TRUE)
    
    m <- mat %*% fac$vectors[, 1:k]
    
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
