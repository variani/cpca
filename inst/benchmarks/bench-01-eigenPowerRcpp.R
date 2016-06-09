library(plyr)
library(ggplot2)

library(microbenchmark)

hilbert <- function(n) { 
  i <- 1:n
  1 / outer(i - 1, i, "+") 
}


#nseq <- seq(500, 1500, length = 5)
nseq <- c(500, 1500, 3000)
df <- ldply(nseq, function(n) {
  n <- ceiling(n)
  v <- runif(n)
  M <- hilbert(n)
  
  cat(" * n:", n, "\n")
  
  out <- microbenchmark(
    eigenPower(M, v), 
    #eigenPowerRcpp(M, v),
    eigenPowerRcppArmadillo(M, v),
    eigenPowerRcppEigen(M, v), 
    times = 10)
  
  df <- subset(as.data.frame(summary(out)), select = c("expr", "median"))
  df$n <- n
  
  return(df)
})  

p1 <- ggplot(df, aes(n, median, color = expr)) + geom_point() + geom_line()
p1

p2 <- ggplot(subset(df, expr != "eigenPowerRcpp(M, v)"), aes(n, median, color = expr)) + geom_point() + geom_line()
p2
