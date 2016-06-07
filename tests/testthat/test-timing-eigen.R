context("Timing for EVD")

test_that("Speed up on a large sparse matrix", {
  ### inc
  library(Matrix)
    
  ### par
  N <- 1250
  prop <- 0.25
  
  ### simulate data 
  Mat <- rsparsematrix(N, N, nnz = ceiling(N * prop), symmetric = T)
  mat <- as.matrix(M2)
  
  out1 <- eigenPower(mat)
  out2 <- eigenPower(Mat)
  
  expect_true(out1$timing$talgo > out2$timing$talgo)
})
