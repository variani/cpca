context("Power method for EVD")

test_that("tiny exampe with 3x3 matrix", {
  ### inc
  library(Matrix)
  
  ### simulate data for testing
  A <- matrix(c(-261, 209, -49, 
    -530, 422, -98,
    -800, 631, -144),
    ncol = 3, nrow = 3, byrow = TRUE)
  v <- c(1, 2, 3)  

  ### compute models
  out1 <- eigenPower(A, v, maxit = 40)
  out2 <- eigenPower(A, v, maxit = 40, sparse = TRUE)
  out3 <- eigenPower(A, v, maxit = 40)
  
  ### testing
  # converged in `< 5` iterations
  expect_true(all(out1$it < 5, out2$it < 5, out3$it < 5)) 
  
  # eigenvalue is 10
  expect_true(all(round(out1$lambda, 2) == 10, round(out2$lambda, 2) == 10,
    round(out3$lambda, 2) == 10)) 

  # eigenvector is [0.27, 0.53, 0.80]
  expect_true(all(all(as.numeric(round(out1$v, 2)) == c(0.27, 0.53, 0.80)), 
    all(as.numeric(round(out2$v, 2)) == c(0.27, 0.53, 0.80)),
    all(as.numeric(round(out3$v, 2)) == c(0.27, 0.53, 0.80))))
})

test_that("Exampe with a matrix of a moderate size", {
  skip("eigenPowerEigen is to be compiled properly")
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
  out2 <- eigenPowerRcpp(mat, v)
  out3 <- eigenPowerRcppEigen(mat, v)

  ### testing
  # converged in `< 15` iterations
  expect_true(all(out1$it < 15, out2$it < 15, out3$it < 15)) 

  # eigenvalue is 2.44
  expect_true(all(round(out1$lambda, 2) == 2.44, round(out2$lambda, 2) == 2.44,
    round(out3$lambda, 2) == 2.44))
})
