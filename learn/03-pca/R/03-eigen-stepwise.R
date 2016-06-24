# @ https://en.wikipedia.org/wiki/Principal_component_analysis
#
# Iterative SVD: X(k) = X - X w w' = X (I - ww')
# Iterative EVD: C(k) = C - C w w' = C (I - ww')

### inc
library(cpca)

### data
mat <- hilbert(4)
cmat <- scale(mat, center = TRUE, scale = FALSE)

covmat <- cov(mat)
cpmat <- crossprod(cmat)
pmat <- crossprod(mat)

### var
p <- ncol(mat)

# not-centered matrices: compare values of EVD and SVD 
evd <- eigen(pmat)
svd <- svd(mat) 
stopifnot(all(round(evd$values, 4) == round(svd$d^2, 4)))

# centered matrices: compare values of EVD and SVD
evd <- eigen(cpmat)
svd <- svd(cmat) 
stopifnot(all(round(evd$values, 4) == round(svd$d^2, 4)))

## EV/SVD: all vectors
evd <- eigen(pmat)
svd <- svd(mat)

### EV1/SVD1
evd1 <- eigen(pmat)
w1 <- evd1$vectors[, 1]

svd1 <- svd(mat)
v1 <- svd1$v[, 1]

prod <- crossprod(w1, v1)
stopifnot(round(abs(prod), 2) == 1)

### EV2/SVD2
pmat2 <- pmat %*% (diag(1, p) - tcrossprod(w1))
evd2 <- eigen(pmat2)

w2 <- evd2$vectors[, 1]

prod <- crossprod(evd$vectors[, 2], w2)
stopifnot(round(abs(prod), 2) == 1)

mat2 <- mat %*% (diag(1, p) - tcrossprod(w1))
svd2 <- svd(mat2)

v2 <- svd2$v[, 1]

prod <- crossprod(svd$v[, 2], v2)
stopifnot(round(abs(prod), 2) == 1)

### EV3/SVD3
pmat3 <- pmat %*% (diag(1, p) - tcrossprod(w1) - tcrossprod(w2))
evd3 <- eigen(pmat3)

w3 <- evd3$vectors[, 1]

prod <- crossprod(evd$vectors[, 3], w3)
stopifnot(round(abs(prod), 2) == 1)

mat3 <- mat %*% (diag(1, p) - tcrossprod(w1) - tcrossprod(w2))
svd3 <- svd(mat3)

v3 <- svd3$v[, 1]

prod <- crossprod(svd$v[, 3], v3)
stopifnot(round(abs(prod), 2) == 1)

