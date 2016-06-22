
#' Function comprcomp. 
#'
#' This function is inspired by the base function \code{prcomp} for PCA and 
#' shares similarities in arguments (\code{center} and \code{scale} arguments),
#' and methods (e.g. \code{scores}).
#'
#' The difference is that \code{comprcomp} computes CPCA by calling \code{cpca} function.
#'
#' The advantage of using \code{comprcomp} instead of \code{cpca} function 
#' is that the input data are given the format observaions vs. variables,
#' and construction of covariance matrices, as well as pre-procesing routines,
#' are implemented inside of \code{comprcomp} function.
#' Consequently, projection of new data, plotting, etc are easier
#' (as easy as using \code{prcomp} function).
#'
#' \code{comprcomp} function is a constructor for \code{comprcomp} class.
#' 
#' @param X A matrix of data points (raw data).
#' @param Y A vector of classes (factor). 
#' @param center A logical value indicating centering (per variable/per column).
#' @param scale A logical value indicating scaling to unit variance (per variable/per column).
#' @param ncomp The numberof components.
#'    The default value is \code{0}, that means computing of all components.
#'
#' @return A object of class \code{comprcomp}.
#' 
#' @export
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
  out <- list(nobs = nobs, p = p, ncomp = ncomp,
    desc = desc, cov = cov, ng = ng,
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
#' @param object An object of class \code{comprcomp}.
#' @param x An object of class \code{comprcomp}.
#' @param X A matrix of data points (raw data).
#' @param Y A vector of classes (factor). 
#' @param comp A vector of components (indcices).
#' @param center A logical value indicating centering (per variable/per column).
#' @param scale A logical value indicating scaling to unit variance (per variable/per column).
#' @param grouping A logical value indicating doing PCA or CPCA (PCA per group).
#'    The default value is \code{FALSE}.
#' @param prop A logical indicating whether the data are to be reported in proportions.
#'    The default value is \code{FALSE}.
#' @param perc A logical indicating whether the data are to be reported in percentage.
#'    The default value is \code{FALSE}.
#' @param sorted A logical indicating wether to sort results or not.
#'    The default value is \code{FALSE}.
#' @param facet A logical indicating whether faceting is to applied in \code{varplot} method.
#' @param ... Additional arguments.
#'
#' @exportClass comprcomp

#' @rdname comprcompClass
#' @export
print.comprcomp <- function(x, ...)
{
  cat("\n Common prcomp object of class:", class(x)[1], "\n")
}

#--------------------------------
# Methods for `comprcomp` class
#--------------------------------

#' @rdname comprcompClass
#' @export
scores.comprcomp <- function(object, X, Y, comp = 1:2, 
  center, scale, grouping = FALSE, ...) 
{
  ### args
  stopifnot(!missing(X))
  stopifnot(class(X) == "matrix")

  stopifnot(!missing(Y))
   
  stopifnot(ncol(X) == object$p)  
  
  stopifnot(all(comp <= object$p))
  
  stopifnot(all(comp <= object$ncomp))
  
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

#' @rdname comprcompClass
#' @export
scores.lda <- function(object, X, Y, comp = 1:2, ...) 
{
  stopifnot(requireNamespace("MASS")) # for MASS:::predict.lda
  
  pred <- predict(object, X)
  scores <- pred$x
  
  scores <- scores[, comp, drop = F]
  
  return(scores)
}
  
#--------------------------------
# Methods `comp*` for components
#--------------------------------

#' @rdname comprcompClass
#' @export
compscore.comprcomp <- function(object, X, 
  sorted = FALSE, ...)
{
  K <- length(object$ng)
  nprop <- as.numeric(object$ng) / sum(object$ng)

  ncomp <- object$ncomp
  
  S <- cov(X)
  
  compscore <- rep(as.numeric(NA), ncomp)
  for(comp in seq(1, ncomp)) {
    q <- object$cpca$CPC[, comp]
    
    sum <- 0
    for(k in 1:K) {
      Sk <- object$cov[[k]]
      
      prod <- as.numeric(crossprod(q, crossprod(S, q)))
      prodk <- as.numeric(crossprod(q, crossprod(Sk, q)))

      sum <- sum + nprop[k] * (log(prod) - log(prodk))
    }
    
    compscore[comp] <- sum
  }
  names(compscore) <- paste0("CPC", seq(1, ncomp))

  ### return
  if(sorted) {
    ord <- order(compscore, decreasing = TRUE)
    compscore <- compscore[ord]
  }
  
  return(compscore)
}

#' @rdname comprcompClass
#' @export
compvar.comprcomp <- function(object, X, Y, comp, 
  center, scale, grouping = FALSE, 
  prop = FALSE, perc = FALSE, sorted = FALSE, ...) 
{
  ### args
  stopifnot(!missing(X))
  stopifnot(class(X) == "matrix")

  stopifnot(!missing(Y))

  stopifnot(ncol(X) == object$p)  
  
  if(missing(comp)) {
    comp <- seq(1, object$ncomp)
  }
  
  stopifnot(all(comp <= object$p))
  stopifnot(all(comp <= object$ncomp))
  
  if(missing(center)) {
    center <- object$center
  }
  
  if(missing(scale)) {
    scale <- object$scale
  }

  if(perc) {
    prop <- TRUE
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
  
  ### compute variances per PC
  if(grouping) {
    compvar <- matrix(as.numeric(NA), nrow = length(comp), ncol = nlevels(Y))
    
    for(i in seq(1, nlevels(Y))) {
      lvl <- levels(Y)[i]
      ind <- object$desc[[lvl]]$ind
      
      Xi <- X[ind, ]
      
      for(k in seq(1, length(comp))) {
        j <- comp[k]

        vec <- object$cpca$CPC[, j]
         
        proj <- as.numeric(tcrossprod(vec, Xi))
        
        compvar[k, i] <- var(proj)
      }
      
      if(prop) {
        vartoti <- sum(diag(var(Xi)))
             
        compvar[, i] <- compvar[, i] / vartoti
      }
    }
    
    rownames(compvar) <- paste0("CPC", comp)
    colnames(compvar) <- levels(Y)
    
  } else {
    compvar <- rep(as.numeric(NA), length(comp))
    
    for(k in seq(1, length(comp))) {
      i <- comp[k]
      
      vec <- object$cpca$CPC[, i]
    
      proj <- as.numeric(tcrossprod(vec, X)) # tcrossprod(vec, X) is a fast version of computing a projection `X * vec`
  
      compvar[k] <- var(proj) 
    }    
    
    if(prop) {
      vartot <- sum(diag(var(X)))
      
      compvar <- compvar / vartot
    }
    
    names(compvar) <- paste0("CPC", comp)
  }

  ### return
  if(perc) {
    compvar <- round(100 * compvar, 1)
  }
  
  if(sorted) {
    ord <- order(compvar, decreasing = TRUE)
    compvar <- compvar[ord]
  }
  
  return(compvar)
}

#-------------------------------------------
# Plotting methods for `comprcomp` class
#-------------------------------------------

#' @rdname comprcompClass
#' @export
varplot.comprcomp <- function(object, X, Y, comp = 1:2, facet = TRUE, ...)
{
  S <- scores(object, X, Y, comp = comp, grouping = TRUE)

  vars <- compvar(object, X, Y, comp = comp, grouping = TRUE, perc = TRUE, ...)
  vars.perc <- c(paste0("(", paste(vars[1, ], collapse = ", "), "%)"),
    paste0("(", paste(vars[2, ], collapse = ", "), " %)"))
    
  ### prepare data.frame `df` for plotting
  df <- as.data.frame(S)
  colnames(df) <- paste0("comp", 1:2)

  df$group <- Y
  
  ### plot
  comp1 <- comp2 <- group <- NULL # no visible binding for global variable
  
  p <- ggplot(df, aes(comp1, comp2, color = group)) + geom_point()

  p <- p + labs(x = paste0("CPC", comp[1], " ", vars.perc[1]), 
    y = paste0("CPC", comp[2], " ", vars.perc[2]))
   
  if(facet) {
    p <- p + facet_wrap(~ group)
  }
  
  return(p)
}


#' @rdname comprcompClass
#' @export
scoreplot.comprcomp <- function(object, X, Y, comp = 1:2, ...)
{
  ### args
  stopifnot(length(comp) == 2)
  
  ### get scores
  S <- scores(object, X, Y, comp = comp, ...)
  
  vars <- compvar(object, X, Y, comp = comp, perc = TRUE, ...)
  vars.perc <- paste0("(", vars, "%)")
  
  ### prepare data.frame `df` for plotting
  df <- as.data.frame(S)
  colnames(df) <- paste0("comp", 1:2)

  df$group <- Y
  
  ### plot
  comp1 <- comp2 <- group <- NULL # no visible binding for global variable
  
  p <- ggplot(df, aes(comp1, comp2, color = group)) + geom_point()
  
  # labs
  p <- p + labs(x = paste0("CPC", comp[1], " ", vars.perc[1]), 
    y = paste0("CPC", comp[2], " ", vars.perc[2]))
  
  return(p)
}

#' @rdname comprcompClass
#' @export
scoreplot.prcomp <- function(object, X, Y, comp = 1:2, ...)
{
  ### args
  stopifnot(length(comp) == 2)
  
  ### get scores
  S <- object$x[, comp]
  
  values <- (object$sdev)^2
  vars <- values[comp] / sum(values)

  vars.perc <- round(100 * vars, 1)
  vars.perc <- paste0("(", vars.perc, "%)")
  
  ### prepare data.frame `df` for plotting
  df <- as.data.frame(S)
  colnames(df) <- paste0("comp", 1:2)

  df$group <- Y
  
  ### plot
  comp1 <- comp2 <- group <- NULL # no visible binding for global variable
  
  p <- ggplot(df, aes(comp1, comp2, color = group)) + geom_point()
  
  # labs
  p <- p + labs(x = paste0("PC", comp[1], " ", vars.perc[1]), 
    y = paste0("PC", comp[2], " ", vars.perc[2]))
  
  return(p)
}

#' @rdname comprcompClass
#' @export
scoreplot.lda <- function(object, X, Y, comp = 1:2, ...)
{
  ### args
  stopifnot(length(comp) == 2)
  
  ### get scores
  S <- scores(object, X, Y, comp = comp, ...)
 
  ### prepare data.frame `df` for plotting
  df <- as.data.frame(S)
  colnames(df) <- paste0("comp", 1:2)

  df$group <- Y
  
  ### plot
  comp1 <- comp2 <- group <- NULL # no visible binding for global variable
  
  p <- ggplot(df, aes(comp1, comp2, color = group)) + geom_point()
  
  # labs
  #p <- p + labs(x = paste0("PC", comp[1], " ", vars.perc[1]), 
  #  y = paste0("PC", comp[2], " ", vars.perc[2]))
  p <- p + labs(x = paste0("LD", comp[1]), y = paste0("LD", comp[2]))
    
  return(p)
}

