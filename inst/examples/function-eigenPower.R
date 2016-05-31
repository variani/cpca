A <- matrix(c(-261, 209, -49, 
  -530, 422, -98,
  -800, 631, -144),
  ncol = 3, nrow = 3, byrow = TRUE)
  
v0 <- c(1, 2, 3)

out <- eigenPower(A, v0, maxit = 40, verbose = 2)

out[c("lambda", "v")]
