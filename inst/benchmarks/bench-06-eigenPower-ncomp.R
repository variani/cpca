library(plyr)
library(ggplot2)

library(microbenchmark)

hilbert <- function(n) { 
  i <- 1:n
  1 / outer(i - 1, i, "+") 
}


nseq <- seq(100, 1200, length = 5)
df <- ldply(nseq, function(n) {
  n <- ceiling(n)
  v <- runif(n)
  M <- hilbert(n)
  
  cat(" * n:", n, "\n")
  
  out <- microbenchmark(
    eigen(M),
    eigenPower(M), 
    eigenPower(M, ncomp = 5),
    eigenPowerRcppEigen(M),
    eigenPowerRcppEigen(M, ncomp = 5),  
    times = 10)
  
  df <- subset(as.data.frame(summary(out)), select = c("expr", "median"))
  df$n <- n
  
  return(df)
})  

p1 <- ggplot(df, aes(n, median, color = expr)) + geom_point() + geom_line()
p1

p2 <- ggplot(subset(df, n > 500), aes(n, median, color = expr)) + geom_point() + geom_line()
p2

