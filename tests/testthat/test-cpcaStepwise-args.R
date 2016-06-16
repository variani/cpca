#------------------------------------------------------------
# Testing arguments of `cpca_stepwise_base` function
#------------------------------------------------------------

context("`cpca_stepwise_base` and its arguments")

test_that("`start` argument", {
  library(plyr)
  
  # data
  data(iris)
  ncomp <- ncol(iris) - 1
  
  cov <- dlply(iris, "Species", function(x) cov(x[, -ncol(x)]))
  ng <- daply(iris, "Species", function(x) nrow(x))

  # run  
  out1 <- cpca_stepwise_base(cov, ng, ncomp = 1, start = "eigen")
  out2 <- cpca_stepwise_base(cov, ng, ncomp = 1, start = "eigenPower")
  out3 <- cpca_stepwise_base(cov, ng, ncomp = 1, start = "random")

  # test
  expect_true(all(round(abs(out1$CPC), 2) == round(abs(out2$CPC), 2)))
  expect_true(all(round(abs(out1$CPC), 2) == round(abs(out3$CPC), 2)))
})
