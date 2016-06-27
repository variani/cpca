# PCA

## Algorithms

Tutorials

* http://scikit-learn.org/stable/modules/decomposition.html

Issues

* [SVD or EVD algorithm for computing PCA?](http://stats.stackexchange.com/questions/79043/why-pca-of-data-by-means-of-svd-of-the-data)
    * Summary from Matlab help: The EIG algorithm is faster than SVD when the number of observations, n, 
      exceeds the number of variables, p, but is less accurate because the condition number of the covariance 
      is the square of the condition number of X.


| computes on | extract all PCs at once | sequential extraction |
|-------------|-------------------------|-----------------------|
| `X`   | SVD | NIPALS |
| `X'X` | EVD | POWER  |


## Application to image compression

* http://www.johnmyleswhite.com/notebook/2009/12/17/image-compression-with-the-svd-in-r/
* https://github.com/road2stat/imgsvd
* http://www.r-bloggers.com/large-scale-eigenvalue-decomposition-and-svd-with-rarpack/

Sample images

* www.lenna.org or http://www.cs.cmu.edu/~chuck/lennapg/lenna.shtml
* http://www.cs.nyu.edu/~roweis/data.html
      * all faces in one figure http://www.cs.nyu.edu/~roweis/data/olivettifaces.gif
      * the dataset in R package http://artax.karlin.mff.cuni.cz/r-help/library/RnavGraphImageData/html/faces.html
      * visualization via t-SNE https://lvdmaaten.github.io/tsne/
      * more description: http://scikit-learn.org/stable/datasets/olivetti_faces.html

Methods similar to CPCA

* Concurent Subspace Analysis (CSA) [link](https://scholar.google.es/scholar?cluster=15503588781614231812&hl=en&as_sdt=2005&sciodt=0,5)
* Multi-block PCA method for image change detection [link](https://scholar.google.es/scholar?cluster=13099656007719521232&hl=en&as_sdt=0,5)
