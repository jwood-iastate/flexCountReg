# test_that("rpnb Tests", {
# 
#   data("washington_roads")
#   expect_no_warning(nb1.rp <- rpnb(Total_crashes ~ - 1 + lnlength + lnaadt,
#                  rpar_formula = ~ speed50,
#                  data = washington_roads,
#                  ndraws = 100,
#                  correlated = FALSE,
#                  rpardists = c(intercept="u", speed50="t"),
#                  form = 'nb1',
#                  method = "bfgs",
#                  print.level = 0))
# 
# })
# 
# # library(testthat)
# 
# 
# 
