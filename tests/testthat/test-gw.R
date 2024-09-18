# PDF ----

test_that("Generalized Waring PDF", {
  
  expect_equal(dgwar(0, mu=1, k=2, rho=3), 0.6)
})

test_that("Generalized Waring CDF", {
  
  pdf_vals <- pgwar(c(0,1,2,3), mu=1, k=2, rho=3)
  expect_true(length(pdf_vals)==4)
})


test_that("Generalized Waring CDF (error)", {
  
  ## Negative P value
  expect_error(pgwar(c(0,1,2,-3), mu=1, k=2, rho=3))
})

# Quantiles ----

test_that("Generalized Waring Quantiles", {
  
  quant <- qgwar(c(0.01,0.1,.5,.8,.99), mu=1, k=2, rho=3)
  
  expect_true(length(quant)==5)
})


test_that("Generalized Waring Quantiles (errors)", {
  
  ## Negative q
  expect_error(qgwar(c(0.01,0.1,.5,.8,-.99), mu=1, k=2, rho=3))
  
  ## Q>1
  expect_error(qgwar(c(0.01,0.1,.5,.8,1.99), mu=1, k=2, rho=3))
})

# Random samples ----

test_that("Generalized Waring Samples", {
  
  set.seed(666)
  gwsamples <- rgwar(100, mu=1, k=2, rho=3)
  expect_true(length(gwsamples)==100)
})




# library(testthat)
# 
# 
# plot(\(x) dtri(x, mode=8, upper=13, lower=1), 0, 14)
# ptri(c(0, 1, 2, 3, 5, 7, 9, 10), mode = 3, upper=9, lower = 1)

# qtri(c(0.1, 0.3, 0.5, 0.9, 0.95), mode = 3, upper = 9, lower = 1)
# rtri(30, mode = 5, sigma = 3)