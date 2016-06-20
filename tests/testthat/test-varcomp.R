context("Variance explained by components")

test_that("`varcomp` method of `comprcomp` class", {
  # data
  data(iris)
  
  X <- as.matrix(iris[, -5])
  Y <- iris[, 5]
  Xc <- scale(X, center = TRUE, scale = FALSE)

  ### model
  mod <- comprcomp(X, Y) # expected: (1) 4 components; (2) centering & not scaling

  ### tests for CPCA model. i.e. `grouping = FALSE`
  vars <- varcomp(mod, X, Y, prop = TRUE)
  expect_true(vars[1] < 0.5) # For PCA, the variance captured by PC1 > 90%
  expect_true(vars[2] > 0.3) # For PCA, the variance captured by PC2 < 6%
  
  ### tests for CPCA model. i.e. `grouping = TRUE`
  vars <- varcomp(mod, X, Y, grouping = TRUE)
  expect_true(all(round(mod$cpca$D, 2) == round(vars, 2)))

  vars <- varcomp(mod, X, Y, grouping = TRUE, prop = TRUE)
  expect_true(all(round(colSums(vars), 2) == 1))
})

