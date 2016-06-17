library(plyr)
library(ggplot2)

library(Matrix)

library(microbenchmark)

N <- 1000
K <- 3
density <- 0.5

n <- 1000

ncomp <- 5

pseq <- seq(10, 150, length = 5)
df <- ldply(pseq, function(p) {
  p <- ceiling(p)
  cat(" * p:", p, "\n")

  out1 <- list()
  out2 <- list()
  
  for(i in 1:K) {
    mat <- matrix(runif(n * p), nrow = n, ncol = p)
    cov <- cov(mat)
    
    cov[cov < 0.01] <- 0
    
    out1[[i]] <- cov
    out2[[i]] <- Matrix(cov, sparse = TRUE)
  }

  ng <- rep(100, K)
  
  out <- microbenchmark(
    cpca_stepwise_base(out1, ng, ncomp = ncomp, start = "random"),
    cpca_stepwise_eigen(out2, ng, ncomp = ncomp, start = "random"),
    times = 10)
  
  df <- subset(as.data.frame(summary(out)), select = c("expr", "median"))
  df$p <- p
  
  return(df)
})  

p <- ggplot(df, aes(p, median, color = expr)) + geom_point() + geom_line()
p

