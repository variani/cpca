context("Timing for EVD")

test_that("CPU gain on a large sparse matrix", {
  ### inc
  library(Matrix)
    
  ### par
  N <- 1250
  prop <- 0.25
  
  ### simulate data 
  M2 <- rsparsematrix(N, N, nnz = ceiling(N * prop), symmetric = T)
  M1 <- as.matrix(M2)
  
  out1 <- eigenPower(M1)
  out2 <- eigenPower(M2)
  
  expect_true(out1$timing$talgo > out2$timing$talgo)
  expect_true(1 == 2)
})

