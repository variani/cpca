context("Power method for CPCA")

test_that("cpc_stepwise_base: ncomp = 1", {
  library(plyr)
  
  # data
  data(iris)

  cov <- dlply(iris, "Species", function(x) cov(x[, -ncol(x)]))
  ng <- daply(iris, "Species", function(x) nrow(x))
  
  out <- cpc_stepwise_base(cov, ng, maxit = 20, ncomp = 1)

  ### tests  
  expect_true(out$converged)
})

test_that("cpc_stepwise_base: start = eigen", {
  library(plyr)
  
  # data
  data(iris)

  cov <- dlply(iris, "Species", function(x) cov(x[, -ncol(x)]))
  ng <- daply(iris, "Species", function(x) nrow(x))
  
  out <- cpc_stepwise_base(cov, ng, maxit = 20, start = "eigen")
  
  ### known results of EVD for iris
  EV <- matrix(c(0.75, 0.09, 0.63, 0.20,
    0.44, -0.79, -0.33, -0.26,
    0.47, 0.60, -0.54, -0.34,
    0.15, -0.02, -0.45, 0.88),
    nrow = 4, ncol = 4, byrow = TRUE)


  ### tests  
  expect_true(out$converged)

  expect_true(all(abs(round(out$CPC, 2)) == abs(EV)))
})

test_that("cpc_stepwise_base: start = random", {
  library(plyr)
  
  # data
  data(iris)

  cov <- dlply(iris, "Species", function(x) cov(x[, -ncol(x)]))
  ng <- daply(iris, "Species", function(x) nrow(x))
  
  out <- cpc_stepwise_base(cov, ng, maxit = 30, start = "random")
  
  ### known results of EVD for iris
  EV <- matrix(c(0.75, 0.09, 0.63, 0.20,
    0.44, -0.79, -0.33, -0.26,
    0.47, 0.60, -0.54, -0.34,
    0.15, -0.02, -0.45, 0.88),
    nrow = 4, ncol = 4, byrow = TRUE)


  ### tests  
  expect_true(out$converged)

  expect_true(all(abs(round(out$CPC, 2)) == abs(EV)))
})
