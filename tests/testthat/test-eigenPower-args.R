#------------------------------------------------------------
# Testing arguments of `eigenPower` function
#------------------------------------------------------------

context("`eigenPower` and its arguments")

test_that("`ncomp` argument", {
  # data
  data(iris)
  ncomp <- ncol(iris) - 1
  
  cov <- cov(iris[, -ncol(iris)])
  
  # run
  out1 <- eigen(cov)
  out2 <- eigenPower(cov, ncomp = 0)
  
  expect_true(out2$ncomp == ncomp)
  
  expect_true(all(round(out1$values, 2) == round(out2$values, 2)))
  expect_true(all(round(abs(out1$vectors), 2) == round(abs(out2$vectors), 2)))
})
