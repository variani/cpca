# cpca

[![travis-ci build status](https://travis-ci.org/variani/cpca.svg?branch=master)](https://travis-ci.org/variani/cpca) [![cran version](http://www.r-pkg.org/badges/version/cpca)](https://cran.r-project.org/web/packages/cpca) [![downloads](http://cranlogs.r-pkg.org/badges/cpca)](http://cranlogs.r-pkg.org/badges/cpca) [![total downloads](http://cranlogs.r-pkg.org/badges/grand-total/cpca)](http://cranlogs.r-pkg.org/badges/grand-total/cpca) [![research software impact](http://depsy.org/api/package/cran/cpca/badge.svg)](http://depsy.org/package/r/cpca)

## About

The `cpca` package approaches Common Principal Component Analysis (CPCA) using [the stepwise method](http://www.sciencedirect.com/science/article/pii/S016794731000112X) proposed by Trendafilov. In contrast to others, this method orders the components by the explained variance intrincically and allows computing a few first components. The later feature is beneficial in practice for high-dimensional data.

![](https://raw.githubusercontent.com/variani/cpca/master/docs/images/lena.png)

The figure above shows an application of CPCA to the image compression problem. The original figure is the famous one from [http://lenna.org](lenna.org), and we estimated CPCA for three RBG data matrices *simultaneously* rather than doing SVD/PCA for each matrix separately. (Spoiler: the performance of PCA and CPCA is the same, and this example is sought for the demonstration purpose only.)

Left panel on the Figure shows the CPCA-based compression using 5 components, central panel - 25 components, and right panel - all components. Script: [learn/03-pca/R/02-compress-image.R](https://github.com/variani/cpca/blob/master/learn/03-pca/R/02-compress-image.R).

The main function in the official release is `cpc`, while a new function `cpca` is used in the development branch.

## CPCA on the iris data

```
library(cpca)
demo(iris, package = "cpca")
```

The [demo.html](http://htmlpreview.github.io/?https://raw.github.com/variani/cpca/master/inst/doc/demo.html) shows that 
the eigenvectors obtained by the `cpc` function are exactly the same as reported 
in [Trendafilov et al., 2010](http://www.sciencedirect.com/science/article/pii/S016794731000112X), Section 5, Example 2. 

## Installation

The following commands install the development (master branch) version from Github.

```
library(devtools)
install_github("cpca", user = "variani")
```

## Citation

Currently, we don't have a specific publication for the `cpca` package. Please see the current citation information by the following command in R.

```
library(cpca)
citation(package = "cpca")
```

## References

A list of publications where the `cpca` package was used:

- Kanaan-Izquierdo *et al.*, Multiview approach to spectral clustering, [EMBC, 2012](https://doi.org/10.1109/EMBC.2012.6346165)
- Fernandez-Albert *et al.*, Intensity drift removal in LC/MS metabolomics by common variance compensation, [Bioinformatics, 2014](https://doi.org/10.1093/bioinformatics/btu423)

## License

The cpca package is licensed under the GPLv3. See COPYING file in the `inst` directory for additional details.

-   COPYING - cpca package license (GPLv3)

