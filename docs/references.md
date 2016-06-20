# Papers

## About

This documents contains a temporal list of references, that were encountered ocasionally and saved for further consideration.

## R packages

* [multigroup](http://cran.r-project.org/web/packages/multigroup/) package: This package includes several methods 
  to study multigroup data, where the same set of variables are measured on different groups of individuals.
    * See function `FCPCA` (that means Flury's Common Principal Component Analysis), where FG algorithm 
      proposed by Flury is implemented.
* [MetaPCA](https://github.com/donkang34/MetaPCA) Simultaneous dimension reduction using PCA when multiple studies are combined. ; available on [CRAN](http://cran.r-project.org/web/packages/MetaPCA/index.html)
  * eigenvalue maximization approach and angle minimization approach
  * extension for Robust PCA and Sparse PCA in the meta-analysis realm
  * [Dongwan Don Kang Google Scholar](http://scholar.google.es/citations?sortby=pubdate&hl=en&user=gU6J7wkAAAAJ&view_op=list_works)

  
## Multivariate Methods

* [Zhang et al., Feature Transformation with Class Conditional Decorrelation, 2013 (pdf)](http://www.nlpr.ia.ac.cn/pal/xyz/Publication/XYZ2013-class_conditional_decorrelation-ICDM.pdf):
  Fisher DA = Whitenning + PCA. Implementation of Whitenning was proposed to be Class Conditional Correlation (CCD). 
 CCD takes a list of class-related covariance matrices and diagonalizes them simultaneously.

## Tests

* [Rublik, A TEST OF THE HYPOTHESIS OF PARTIAL COMMON PRINCIPAL COMPONENTS, 2009 (pdf)](http://www.um.sav.sk/en/images/stories/dep03/doc/FRHypothesisPCPC.pdf):
  See Introduction for literature review.
* [Boente et al., Robust tests for the common principal components model, 2009 (pdf)](http://www.researchgate.net/publication/222548516_Robust_tests_for_the_common_principal_components_model/file/9fcfd50df235394ee4.pdf):
  See Section 7 Final Comments for a good summary.
* Diego et al., (2014) GAW19 contribution.

PCA tests

* Duscussion: http://stats.stackexchange.com/questions/33917/how-to-determine-significant-principal-components-using-bootstrapping-or-monte-c
* R package for permutation tests: http://cran.r-project.org/web/packages/jackstraw/jackstraw.pdf
* CCA tests (two groups needed): http://cran.r-project.org/web/packages/CCP/CCP.pdf
  * Example of application: http://alstatr.blogspot.com.es/2015/01/canonical-correlation-analysis-on.html
* Articles
    * http://pbil.univ-lyon1.fr/JTHome/CD/articles/SD805.pdf
    * Lassoed PC: http://projecteuclid.org/euclid.aoas/1223908049
