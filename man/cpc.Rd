\name{cpc}
\alias{cpc}
\title{Function cpc.}
\usage{
  cpc(X, method = "stepwise", k = 0, iter = 30,
    threshold = 0, ...)
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

  \item{iter}{The maximum number of iterations.  The
  parameter is valid for the stepwise algorithm by
  Trendafilov, that is applied in the power algorithm for
  estimation a single component.  The default value is 30.}

  \item{threshold}{The threshold value of the captured
  variance, which is reserved for further extensions.}

  \item{...}{Other parameters.}
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
\note{
  This function adpats the original code in matlab written
  by Dr N. T. Trendafilov.
}
\examples{
require(plyr)
require(abind)

data(iris)

C <- daply(iris, "Species", function(x) cov(x[, -ncol(x)]))
C <- aperm(C, c(2, 3, 1)) # put the 1st dimension to the end

# default call
mod1 <- cpc(C)
round(mod1$CPC, 2)

# compute only first two CPCs
mod2 <- cpc(C, k = 2)
round(mod2$CPC, 2)

}
\references{
  Trendafilov (2010). Stepwise estimation of common
  principal components.  Computational Statistics & Data
  Analysis, 54(12), 3446-3457.
  doi:10.1016/j.csda.2010.03.010
}

