### inc 
library(plyr)

library(pls)
library(MASS)

library(png)

library(ggplot2)
library(gridExtra)

theme_set(theme_linedraw())

#plot(1:64, 1:64, type = 'n')
#rasterImage(img, 1, 1, 64, 64)

plotFaces <- function(dat, ncol = 10) 
{
  grobs <- llply(1:nrow(dat), function(i) {
    rasterGrob(matrix(dat[i, ], 64, 64))
  })
  p <- marrangeGrob(grobs, nrow = ceiling(nrow(dat) / ncol), ncol = ncol, top = NULL)
  
  return(p)
}

### data
data(faces, package = "RnavGraphImageData")

faces <- faces / 255 # 0-255 values, i.e. 8-bit encoding

faces <- t(faces)

### prepare data
X1 <- faces

persons <- paste0("P", 1:40)
Y1 <- factor(rep(persons, each = 10), levels = persons)

m1 <- prcomp(X1)

p1 <- scoreplot(m1, Y = Y1)

### subset of faces
X2 <- faces[1:50, ]

persons <- paste0("P", 1:5)
Y2 <- factor(rep(persons, each = 10), levels = persons)

m21 <- prcomp(X2)

p21 <- scoreplot(m21, Y = Y2)

m22 <- lda(X2, Y2)
p22 <- scoreplot(m22, Y = Y2)

### show plots
grid.arrange(p21, p22)

