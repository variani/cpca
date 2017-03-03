# @ http://stats.stackexchange.com/questions/134282/relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca

### inc
library(rsvd)

### data
data(iris)

dat <- iris[, -5]
dat_c <- scale(dat, center = TRUE, scale = FALSE)

### PCA
mod_pca <- prcomp(dat_c, center = FALSE, scale. = FALSE)

ev_pca <- mod_pca$rotation # eigen vectors or PCs
sc_pca <- mod_pca$x # scores
scx_pca <- dat_c %*% ev_pca

### SVD
mod_svd <- svd(dat_c)

ev_svd <- mod_svd$v
sc_svd <- mod_svd$u %*% diag(mod_svd$d)
scx_svd <- dat_c %*% ev_svd

### Random SVD
mod_rsvd <- rsvd(dat_c, k = 2)

ev_rsvd <- mod_rsvd$v

### test
stopifnot(all(ev_pca - ev_svd < 1e-10))
stopifnot(all(sc_pca - sc_svd < 1e-10))
stopifnot(all(scx_pca - scx_svd < 1e-10))

stopifnot(all(abs(ev_pca[, 1:2]) - abs(ev_rsvd) < 1e-10))

