library(plyr)
library(ggplot2)

library(microbenchmark)

hilbert <- function(n) { 
  i <- 1:n
  1 / outer(i - 1, i, "+") 
}


nseq <- seq(50, 500, length = 5)
df <- ldply(nseq, function(n) {
  n <- ceiling(n)
  v <- runif(n)
  M <- hilbert(n)
  
  cat(" * n:", n, "\n")
  
  out <- microbenchmark(
    eigen(M),
    eigen(M, symmetric = TRUE),
    eigenPowerRcppEigen(M),
    eigen_RcppEigen(M),
    eigenSelfAdjoint_RcppEigen(M),
    times = 10)
  
  df <- subset(as.data.frame(summary(out)), select = c("expr", "median"))
  df$n <- n
  
  return(df)
})  

p <- ggplot(df, aes(n, median, color = expr)) + geom_point() + geom_line()
p

