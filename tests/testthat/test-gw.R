# PDF ----

test_that("Generalized Waring PDF", {
  
  expect_equal(dgwar(0, mu=1, k=2, rho=3), 0.6)
})

test_that("Generalized Waring CDF", {
  
  pdf_vals <- pgwar(c(0,1,2,3), mu=1, k=2, rho=3)
  expect_true(length(pdf_vals)==4)
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



# Changes to Generalized Waring ----
test_that("Generalized Waring accepts rho = 2", {
  y <- 0:20
  pmf <- dgwar(y, mu = 1, k = 2, rho = 2)
  cdf <- pgwar(y, mu = 1, k = 2, rho = 2)
  
  # rho = 2 has an infinite variance, but the PMF and mean are valid.
  expect_true(all(is.finite(pmf)))
  expect_equal(pmf[1], 0.6857142857142857, tolerance = 1e-12)
  expect_equal(cdf, cumsum(pmf), tolerance = 1e-10)
  expect_true(sum(pmf) < 1) # the omitted tail is positive
})

test_that("Generalized Waring uses one consistent rho domain", {
  # rho > 1 is sufficient for the mean parameterization.
  expect_true(is.finite(dgwar(0, mu = 1, k = 2, rho = 1.01)))
  expect_true(is.finite(pgwar(0, mu = 1, k = 2, rho = 1.01)))
  
  # rho <= 1 does not define the requested mean parameterization.
  expect_error(dgwar(0, mu = 1, k = 2, rho = 1))
  expect_error(pgwar(0, mu = 1, k = 2, rho = 1))
  expect_error(dgwar(0, mu = 1, k = 2, rho = 0.5))
})

test_that("Generalized Waring rejects invalid k before gamma evaluation", {
  for (bad_k in c(0, -1, NA_real_, NaN, Inf, -Inf)) {
    expect_error(dgwar(0, mu = 1, k = bad_k, rho = 3))
    expect_error(pgwar(0, mu = 1, k = bad_k, rho = 3))
  }
})

test_that("Generalized Waring PMF and CDF agree at rho = 2", {
  y <- 0:6
  expect_equal(
    dgwar(y, mu = 1, k = 2, rho = 2),
    diff(c(0, pgwar(y, mu = 1, k = 2, rho = 2))),
    tolerance = 1e-10
  )
})

