### inc 
library(png)

library(ggplot2)
library(grid)
library(gridExtra)

### par
imgfile <- "data/lena_std.png"
imgfile2 <- "data/img-rarpack.png"

### local fucntions
compress_pca <- function(m, k, center = TRUE) 
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

compress_cpca <- function(img, k, center = TRUE) 
{
  K <- dim(img)[3]
  
  ### prepare data
  means <- list()
  cov <- list()
  ng <- rep(as.integer(NA), K)
  
  for(i in 1:K) {
    cat(" * ", i, "/", K, "\n")
    
    m <- img[, , i]
    
    meansi <- colMeans(m)
    if(center) {
      m <- sweep(m, 2, meansi, "-") 
    }
    
    means[[i]] <- meansi
    cov[[i]] <- crossprod(m)
    ng[i] <- nrow(m)
  }
  
  ### run CPCA  
  out <- cpca(cov, ng, ncomp = k)
    
  ### reconstruct
  dat <- array(as.numeric(NA), dim(img))
  for(i in 1:K) {
    cat(" * ", i, "/", K, "\n")
    
    m <- img[, , i]
    
    meansi <- means[[i]]
    
    if(center) {
      m <- sweep(m, 2, meansi, "-") 
    }
    
    mat <- m %*% out$CPC[, 1:k, drop = FALSE] %*% t(out$CPC[, 1:k, drop = FALSE])
    
    if(center) {
      mat <- sweep(mat, 2, means[[i]], "+") 
    }

    mat[mat < 0] <- 0
    mat[mat > 1] <- 1
        
    dat[, , i] <- mat
  }
  
  return(dat)
}

viewImg <- function(img) 
{
  marrangeGrob(list(rasterGrob(img)), nrow = 1, ncol = 1, top = NULL)
}  

lineImg <- function(img, k, ncol = 10, center = TRUE) 
{
  imgs <- c(llply(k, function(x) compress_pca(img, x, center = center)), list(img),
    llply(k, function(x) compress_cpca(img, x, center = center)), list(img))
    
  marrangeGrob(llply(imgs, rasterGrob),
    nrow = 2, ncol = length(k) + 1, top = NULL)  
}

### read file
img <- readPNG(imgfile)
#img <- readPNG(imgfile2)

