### inc 
library(pls)

load_all()

### par
center <- TRUE
scale <- FALSE

### data
data(iris)
X <- as.matrix(iris[, -5])
Y <- iris[, 5]

### PCA model
Xc <- scale(X, center = center, scale = scale)
mod <- prcomp(X, center = center, scale = scale)

covX <- cov(X)

### eigen values
values <- mod$sdev^2
round(values, 2)

# variance of scores
round(cov(mod$x), 2)

# computing scores by hand
round(var(Xc %*% mod$rotation), 2)

### captured variance via eigenvalues
round(100 * values / sum(values), 1)
#[1] 92.5  5.3  1.7  0.5

### commond PCA model
m <- comprcomp(X, Y) 

rot <- m$cpca$CPC

# computing scores by hand
round(var(Xc %*% rot), 2)

sum(var(Xc %*% mod$rotation))
sum(diag(var(Xc %*% rot)))
sum(diag(var(Xc)))

values2 <- diag(var(Xc %*% rot))
round(values2, 2)
round(values, 2)




