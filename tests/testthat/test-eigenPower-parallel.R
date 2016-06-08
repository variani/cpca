context("eigePower in parallel")

test_that("eigenPowerRcppParallel: the same results on 1 it", {
  library(RcppParallel)
  
  ### par  
  n <- 1000
  digits <- 5
  
  cores <- defaultNumThreads()
  if(cores > 1) {
    cores <- 2
    
    ### data  
    hilbert <- function(n) { 
      i <- 1:n
      1 / outer(i - 1, i, "+") 
    }

    mat <- hilbert(n)
    v <- 1:n
  
    ### compute models
    out0 <- eigenPower(mat, v, maxit = 1)
    
    out1 <- eigenPowerRcppParallel(mat, v, maxit = 1, cores = cores)
    out2 <- eigenPowerRcppParallel(mat, v, maxit = 1, cores = cores)
    out3 <- eigenPowerRcppParallel(mat, v, maxit = 1, cores = cores)

    ### testing
    # head of eigen vector: 0.05991430 0.05952567 0.05919683 0.05889783 0.05861869 0.05835440
    expect_true(all(head(round(out1$v, digits)) == head(round(out0$v, digits))))
    expect_true(all(head(round(out2$v, digits)) == head(round(out0$v, digits))))
    expect_true(all(head(round(out3$v, digits)) == head(round(out0$v, digits))))
  }
})

test_that("eigenPowerEigenRcppParallel: the same results on 1 it", {
  skip("eigenPowerEigenParallel is not working")
  
  library(RcppParallel)

  ### par  
  n <- 1000
  digits <- 5
      
  cores <- defaultNumThreads()
  if(cores > 1) {
    cores <- 2
    
    ### data  
    hilbert <- function(n) { 
      i <- 1:n
      1 / outer(i - 1, i, "+") 
    }

    mat <- hilbert(n)
    v <- 1:n
  
    ### compute models
    out0 <- eigenPower(mat, v, maxit = 1)
    
    out1 <- eigenPowerEigenParallel(mat, v, maxit = 1, cores = cores)
    out2 <- eigenPowerEigenParallel(mat, v, maxit = 1, cores = cores)
    out3 <- eigenPowerEigenParallel(mat, v, maxit = 1, cores = cores)

    ### testing
    # head of eigen vector: 0.05991430 0.05952567 0.05919683 0.05889783 0.05861869 0.05835440
    expect_true(all(all(head(round(out1$v, digits)) == head(round(out0$v, digits)))),
      all(head(round(out2$v, digits)) == head(round(out0$v, digits))),
      all(head(round(out3$v, digits)) == head(round(out0$v, digits))))
  }
})
