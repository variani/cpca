context("CPCA applied to iris dataset")

test_that("cpca on iris", {
  library(plyr)
  library(abind)
  
  data(iris)

  C <- daply(iris, "Species", function(x) cov(x[, -ncol(x)]))
  C <- aperm(C, c(2, 3, 1)) # put the 1st dimension to the end

  mod <- cpc(C)
  V.cpc <- mod$CPC # eigen vectors
  
  # See Trendafilov (2010). Stepwise estimation of common principal components. 
  # Computational Statistics & Data Analysis, 54(12), 3446-3457. 
  # doi:10.1016/j.csda.2010.03.010
  # p. 10, Example 2
  V.article <- matrix(c(0.75, -0.09, 0.63, 0.20,
    0.44, 0.79, -0.33, -0.26,
    0.47, -0.60, -0.54, -0.34,
    0.15, 0.02, -0.45, 0.88), 
    nrow = 4, ncol = 4, byrow = TRUE)

  expect_true(all(abs(round(V.cpc, 2)) == abs(V.article)))
})


