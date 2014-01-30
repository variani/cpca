\name{cpc}
\alias{cpc}
\title{Function cpc.}
\usage{
  cpc(X, method = "stepwise", k = 0, threshold = 0, ...)
}
\arguments{
  \item{X}{An array of three dimensions: the 3rd dimension
  encodes the groups and the first two dimension contain
  the covariance matrices.}

  \item{method}{The name of the method for computing the
  CPCA.  The default value is \code{"stepwise"}, which is
  the stepwise algorithm by Trendafilov.}

  \item{k}{The number of components to be computed (all if
  it is \code{0}).  This paramter is valid if the given
  method supports built-in ordering of the eigvenvectors.
  The default value is \code{0}, that means computing of
  all the components.}

  \item{threshold}{The threshold value of the captured
  variance, which is reserved for further extensions.}
}
\value{
  A list several slots: \code{CPC} rotation matrix with
  eigenvectors in columns; \code{ncomp} the number of
  components evaluated (equal to the number of columns in
  \code{CPC}).
}
\description{
  This function computes the CPCA from a given set of
  covariance matrices (of different groups).
}
\details{
  Currently, the only the stepwise algorithm by Trendafilov
  is supported.
}
\examples{
library(plyr)
library(abind)

data(iris)

C <- daply(iris, "Species", function(x) cov(x[, -ncol(x)]))

C <- aperm(C, c(2, 3, 1)) # put the 1st dimension to the end
dim(C)
dimnames(C)

mod <- cpc(C)
str(mod)

round(mod$CPC, 2)
# See Trendafilov (2010). Stepwise estimation of common principal components. 
# Computational Statistics & Data Analysis, 54(12), 3446–3457. 
# doi:10.1016/j.csda.2010.03.010
# p. 10, Example 2
#
#     [,1]  [,2]  [,3]  [,4]
#[1,] 0.75 -0.09  0.63  0.20
#[2,] 0.44  0.79 -0.33 -0.26
#[3,] 0.47 -0.60 -0.54 -0.34
#[4,] 0.15  0.02 -0.45  0.88
#
# The eigenvectors must be the same, as the default method in `cpc` function
# is the power algorithm proposed by Trendafilov.
}
\references{
  Trendafilov (2010). Stepwise estimation of common
  principal components.  Computational Statistics & Data
  Analysis, 54(12), 3446–3457.
  doi:10.1016/j.csda.2010.03.010
}

