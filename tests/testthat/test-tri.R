# PDF ----

test_that("Triangular PDF", {
  
  ## Riemann sum:
  delta <- 0.1
  x <- seq(0, 14, by = delta)
  y <- dtri(x, mode = 8, upper = 13, lower = 1)
  out <- delta * sum(y) 
  
  expect_equal(out, 1)
})

test_that("Triangular PDF (error 1)", {
  
  ## Lower > upper
  expect_error(dtri(8, mode = 8, upper = 13, lower = 14))
})

test_that("Triangular PDF (error 2)", {
  
  ## Mode > upper
  expect_error(dtri(8, mode = 14, upper = 13, lower = 1))
})

test_that("Triangular PDF (error 3)", {
  
  ## Mode < lower
  expect_error(dtri(8, mode = 0, upper = 13, lower = 1))
})

# PDF ----

test_that("Triangular CDF", {
  
  pdf_vals <- ptri(c(0, 1, 3, 9, 10), mode = 3, upper=9, lower = 1)
  pdf_expected <- c(0.00, 0.00, 0.25, 1.00, 1.00)
  out <- all(pdf_expected == pdf_vals)
  
  expect_true(out)
})


test_that("Triangular CDF (errors)", {
  
  ## Lower > upper
  expect_error(ptri(8, mode = 8, upper = 13, lower = 14))
  
  ## Mode > upper
  expect_error(ptri(8, mode = 14, upper = 13, lower = 1))
  
  ## Mode < lower
  expect_error(ptri(8, mode = 0, upper = 13, lower = 1))
})

# Quantiles ----

test_that("Triangular Quantiles", {
  
  quant <- qtri(c(0, 0.5, 1), mode = 5, upper=9, lower = 1)
  quant_expected <- c(1, 5, 9)
  out <- all(quant_expected == quant)
  
  expect_true(out)
})


test_that("Triangular Quantiles (errors)", {
  
  ## Lower > upper
  expect_error(qtri(8, mode = 8, upper = 13, lower = 14))
  
  ## Mode > upper
  expect_error(qtri(8, mode = 14, upper = 13, lower = 1))
  
  ## Mode < lower
  expect_error(qtri(8, mode = 0, upper = 13, lower = 1))
})

# Random samples ----

test_that("Triangular Samples", {
  
  set.seed(666)
  trisamples <- rtri(1e4, mode = 2, upper = 9, lower = 1)
  sample_median <- median(trisamples)
  median <- 9 - sqrt(28)
  diff <- abs(sample_median - median) 
  
  testthat::expect_true(diff < 0.04)
})


test_that("Triangular Samples (errors)", {
  
  ## Lower > upper
  expect_error(rtri(10, mode = 8, upper = 13, lower = 14))
  
  ## Mode > upper
  expect_error(rtri(10, mode = 14, upper = 13, lower = 1))
  
  ## Mode < lower
  expect_error(rtri(10, mode = 0, upper = 13, lower = 1))
})




# library(testthat)
# 
# 
# plot(\(x) dtri(x, mode=8, upper=13, lower=1), 0, 14)
# ptri(c(0, 1, 2, 3, 5, 7, 9, 10), mode = 3, upper=9, lower = 1)

# qtri(c(0.1, 0.3, 0.5, 0.9, 0.95), mode = 3, upper = 9, lower = 1)
# rtri(30, mode = 5, sigma = 3)