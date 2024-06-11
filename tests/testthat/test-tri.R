test_that("Triangular PDF", {
  
  ## Riemann sum:
  delta <- 0.1
  x <- seq(0, 14, by = delta)
  y <- dtri(x, mode = 8, upper = 13, lower = 1)
  out <- delta * sum(y) 
  
  expect_equal(out, 1)
})


# library(testthat)
# 
# 
# plot(\(x) dtri(x, mode=8, upper=13, lower=1), 0, 14)
# 
# ptri(c(0, 1, 2, 3, 5, 7, 9, 10), mode = 3, upper=9, lower = 1)
# qtri(c(0.1, 0.3, 0.5, 0.9, 0.95), mode = 3, upper = 9, lower = 1)
# rtri(30, mode = 5, sigma = 3)