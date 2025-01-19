test_that("rpnb Tests", {

  data("washington_roads")
  expect_no_warning(nbp.rp <- rpnb(Total_crashes ~ - 1 + lnlength + lnaadt,
                 rpar_formula = ~ speed50,
                 data = washington_roads,
                 ndraws = 100,
                 correlated = FALSE,
                 rpardists = c(intercept="n", speed50="n"),
                 form = 'nbp',
                 method = "bfgs"))

})

# # library(testthat)
# 
# 
# 
