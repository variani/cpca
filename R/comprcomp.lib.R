
comprcomp <- function(X, Y, center = TRUE, scale = FALSE, ncomp = 0)
{
  ### arg
  stopifnot(!missing(Y))
  stopifnot(!missing(X))
  stopifnot(class(X) == "matrix")
   
  p <- ncol(X)
  nobs <- nrow(X)
  
  if(ncomp == 0) { 
    ncomp <- p
  }
  
  ### pre-processing
  desc <- llply(levels(Y), function(lvl) {
    ind = which(Y == lvl)
    list(ind = ind, nobs = length(ind),
      mean = apply(X[ind, ], 2, mean),
      sd = apply(X[ind, ], 2, sd))
  })
  names(desc) <- levels(Y)
  
  # compute covariances: (1) update `X`; (2) compute a list of covariance matrices (per group/level of `Y`) `cov`
  cov <- list()
  for(lvl in levels(Y)) {
    ind <- desc[[lvl]]$ind
    
    if(center) { 
      X[ind, ] <- sweep(X[ind, ], 2, desc[[lvl]]$mean, "-") 
    }
    if(scale) { 
      X[ind, ] <- sweep(X[ind, ], 2, desc[[lvl]]$sd, "/") 
    }
    
    cov[[lvl]] <- (1.0 / (desc[[lvl]]$nobs - 1)) * crossprod(X[ind, ])
  }
  
  # compute group size: vector `ng`
  ng <- laply(desc, function(x) x$nobs)
  
  ### EVD
  cpca <- cpca(cov, ng, ncomp = ncomp)
  
  ### output
  out <- list(desc = desc, cov = cov, ng = ng,
    center = center, scale = scale,
    cpca = cpca)
  
  oldClass(out) <- "comprcomp"
  return(out)
}

#--------------------
# comprcomp class
#--------------------

#' S3 class comprcomp.
#'
#' @name comprcompClass
#' @rdname comprcompClass
#'
#' @param object 
#'    An object of class \code{comprcomp}.
#' @param ...
#'    Additional arguments.
#'
#' @exportClass comprcomp

#--------------------------------
# Methods for `comprcomp` class
#--------------------------------

## @rdname comprcompClass
## @export scoreplot
##scores <- function(object, ...) UseMethod("scores")

#' @rdname comprcompClass
#' @export
scores.comprcomp <- function(object, X, Y, comp = 1:2, 
  center, scale, grouping = FALSE, ...) 
{
  ### args
  stopifnot(!missing(Y))

  stopifnot(!missing(X))
  stopifnot(class(X) == "matrix")
   
  p <- ncol(X)
  
  stopifnot(all(comp <= p))
  stopifnot(all(comp <= object$cpca$ncomp))
  
  if(missing(center)) {
    center <- object$center
  }
  
  if(missing(scale)) {
    scale <- object$scale
  }
  
  ### pre-process data in `X`
  if(grouping) {
    for(lvl in levels(Y)) {
      ind <- object$desc[[lvl]]$ind
    
      if(center) { 
        X[ind, ] <- sweep(X[ind, ], 2, object$desc[[lvl]]$mean, "-") 
      }
      if(scale) { 
        X[ind, ] <- sweep(X[ind, ], 2, object$desc[[lvl]]$sd, "/") 
      }
    }
  } else {
    meanX <- apply(X, 2, mean)
    sdX <- apply(X, 2, sd)  
  
    if(center) { 
      X <- sweep(X, 2, meanX, "-") 
    }
    if(scale) { 
      X <- sweep(X, 2, sdX, "/") 
    }
  }
  
  # output matrix of scores `S` 
  S <- X %*% object$cpca$CPC[, comp, drop = F]
  
  return(S)
}

#-------------------------------------------
# Plotting methods for `comprcomp` class
#-------------------------------------------

#' @rdname comprcompClass
#' @export varplot
varplot <- function(object, ...) UseMethod("varplot")

#' @rdname comprcompClass
#' @export
varplot.comprcomp <- function(object, X, Y, comp = 1:2, facet = TRUE, ...)
{
  S <- scores(object, X, Y, comp = comp, grouping = TRUE)
  
  ### prepare data.frame `df` for plotting
  df <- as.data.frame(S)
  colnames(df) <- paste0("comp", 1:2)

  df$lab <- Y
  
  ### plot
  p <- ggplot(df, aes(comp1, comp2, color = lab)) + geom_point()
  
  if(facet) {
    p <- p + facet_wrap(~ lab)
  }
  
  return(p)
}

## @rdname comprcompClass
## @export scoreplot
##scoreplot <- function(object, ...) UseMethod("scoreplot")

#' @rdname comprcompClass
#' @export
scoreplot.comprcomp <- function(object, X, Y, comp = 1:2, ...)
{
  S <- scores(object, X, Y, comp = comp)
  
  ### prepare data.frame `df` for plotting
  df <- as.data.frame(S)
  colnames(df) <- paste0("comp", 1:2)

  df$lab <- Y
  
  ### plot
  p <- ggplot(df, aes(comp1, comp2, color = lab)) + geom_point()
  
  return(p)
}

