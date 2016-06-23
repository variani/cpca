# PCA

## Algorithms

* [SVD or EVD algorithm for computing PCA?](http://stats.stackexchange.com/questions/79043/why-pca-of-data-by-means-of-svd-of-the-data)
    * Summary from Matlab help: The EIG algorithm is faster than SVD when the number of observations, n, 
      exceeds the number of variables, p, but is less accurate because the condition number of the covariance 
      is the square of the condition number of X.


```
computes on      extract all PCs at once       sequential extraction    
X                SVD                           NIPALS    
X'X              EVD                           POWER
```

## Application to image compression

* http://www.johnmyleswhite.com/notebook/2009/12/17/image-compression-with-the-svd-in-r/
* https://github.com/road2stat/imgsvd
* http://www.r-bloggers.com/large-scale-eigenvalue-decomposition-and-svd-with-rarpack/
