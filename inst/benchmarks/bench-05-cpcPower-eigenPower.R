library(plyr)
library(ggplot2)

library(Matrix)

library(microbenchmark)

N <- 1000
K <- 3
density <- 0.5

ncomp1 <- 1
ncomp2 <- 3

pseq <- seq(10, 500, length = 5)
df <- ldply(pseq, function(p) {
  p <- ceiling(p)
  cat(" * p:", p, "\n")

  stopifnot(ncomp1 <= p)
  stopifnot(ncomp2 <= p)  

  out <- list()
  
  for(i in 1:K) {
    mat <- matrix(runif(N * p), nrow = N, ncol = p)
    cov <- cov(mat)
    
    out[[i]] <- cov
  }

  ng <- rep(100, K)
  
  out <- microbenchmark(
    cpca_stepwise_base(out, ng, ncomp = ncomp1, start = "eigen"),
    cpca_stepwise_base(out, ng, ncomp = ncomp2, start = "eigen"),
    cpca_stepwise_base(out, ng, ncomp = ncomp1, start = "eigenPower"),
    cpca_stepwise_base(out, ng, ncomp = ncomp2, start = "eigenPower"),        
    times = 5)
  
  df <- subset(as.data.frame(summary(out)), select = c("expr", "median"))
  df$p <- p
  
  return(df)
})  

p <- ggplot(df, aes(p, median, color = expr)) + geom_point() + geom_line()
p

