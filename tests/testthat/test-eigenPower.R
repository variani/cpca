context("Power method for EVD")

test_that("tiny exampe with 3x3 matrix", {
  ### inc
  library(Matrix)
  
  ### simulate data for testing
  A <- matrix(c(-261, 209, -49, 
    -530, 422, -98,
    -800, 631, -144),
    ncol = 3, nrow = 3, byrow = TRUE)
  v0 <- c(1, 2, 3)  

  ### compute models
  out1 <- eigenPower(A, v0, maxit = 40)
  out2 <- eigenPower(A, v0, maxit = 40, sparse = TRUE)

  ### testing
  # converged in 2 iterations
  expect_true(all(out1$it == 2, out2$it == 2)) 
  
  # eigenvalue is 10
  expect_true(all(round(out1$lambda, 2) == 10, round(out2$lambda, 2) == 10)) 

  # eigenvector is [0.27, 0.53, 0.80]
  expect_true(all(all(as.numeric(round(out1$v, 2)) == c(0.27, 0.53, 0.80)), 
    all(as.numeric(round(out2$v, 2)) == c(0.27, 0.53, 0.80))))
})

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
})



