library(testthat)
library(cpca)

pattern <- "^(?!.*?timing).*" # @ http://stackoverflow.com/a/3688210/551589
test_check("cpca", filter = pattern, perl = TRUE)


