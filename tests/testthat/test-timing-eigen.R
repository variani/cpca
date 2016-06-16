context("Timing for EVD")

test_that("Speed up on a large sparse matrix", {
  skip("better benchmarks are needed")
  
  ### inc
  library(Matrix)
    
  ### par
  N <- 1250
  prop <- 0.25
  
  ### simulate data 
  Mat <- rsparsematrix(N, N, nnz = ceiling(N * prop), symmetric = T)
  mat <- as.matrix(Mat)
  
  out1 <- eigenPower(mat)
  out2 <- eigenPower(Mat)
  
  expect_true(out1$timing$talgo > out2$timing$talgo)
})

test_that("Speed up on 2 cores", {
  skip("eigenPowerEigenParallel is not working")
  
  library(RcppParallel)
  
  cores <- defaultNumThreads()

  if(cores > 1) {
    cores <- 2
    
    ### par  
    n <- 1000
  
    hilbert <- function(n) { 
      i <- 1:n
      1 / outer(i - 1, i, "+") 
    }

    mat <- hilbert(n)
    v <- 1:n
  
    ### compute models
    out1 <- eigenPower(mat, v)
    out2 <- eigenPowerRcppParallel(mat, v, cores = cores)
    out3 <- eigenPowerEigenParallel(mat, v, cores = cores)

    ### testing
    # converged in `< 20` iterations
    #expect_true(all(out1$it < 20, out2$it < 20, out3$it < 20)) 

    # eigenvalue is 2.44
    #expect_true(all(round(out1$lambda, 2) == 2.44, round(out2$lambda, 2) == 2.44,
    #  round(out3$lambda, 2) == 2.44))
  }
})
