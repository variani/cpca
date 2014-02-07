# cpca

`cpca` is an R package with methods to perform Common Principal Component Analysis (CPCA).

The main function to perform CPCA is called `cpc`. See `?cpc` for the help.

For now, the `cpc` function implements only one method based on Trendafilov, 2010.
This method estimates the Common Principal Components (CPCs) by a stepwise procedure 
based on the well-known power method for a single covariance/correlation matrix.
The feature of this method is that it orders the CPCs by the explained variance (intrincically),
and the user can estimate the few first components, e.g. 2-3, rather than all the components.
It is beneficial in practice when a data set has many variables.


## Demo

The `iris` demo shows an application of the `cpc` function to Fisher's iris data. 

```
library(cpca)
demo(iris, package = "cpca")
```

[demo.html](http://htmlpreview.github.io/?https://raw.github.com/variani/cpca/master/inst/doc/demo.html) stored in the `inst/doc` directory presents both the code and the resulted output of the demo.

Note that the eigenvectors obtained by the `cpc` function are exactly the same as reported in Trendafilov, 2010, Section 5, Example 2. That means that Trendafilov's method (which is default in the `cpc` function) is implemnted accurately (at least for iris data).  

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

The citation information is stored in the `CITATION` file in the `inst` directory and can be updated in the future.

* CITATION - cpca package citation information

## References

List of publications, where the `cpca` package was used:

* Kanaan-Izquierdo, S., Ziyatdinov, A., Massanet, R., & Perera, A. (2012). Multiview approach to spectral clustering.  In 2012 Annual International Conference of the IEEE Engineering in Medicine and Biology Society (pp. 1254–1257). IEEE. doi:10.1109/EMBC.2012.6346165
* Fernandez-Albert, F. et al. (to be appeared). A Common Variance Compensation method for intensity drift removal in LC / MS metabolomics.

Mathematical algorithms implemented in the `cpca` package:

* Trendafilov, N. T. (2010). Stepwise estimation of common principal components. Computational Statistics & Data Analysis, 54(12), 3446–3457. doi:10.1016/j.csda.2010.03.010

## License

The cpca package is licensed under the GPLv3. See COPYING file in the `inst` directory for additional details.

* COPYING - cpca package license (GPLv3)
