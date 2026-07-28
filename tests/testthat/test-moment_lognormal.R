test_that("moment_lognormal calculates lognormal raw moments", {
  expect_equal(
    moment_lognormal(mu = 0, sigma = 1, n = 1),
    exp(1 / 2)
  )
  
  expect_equal(
    moment_lognormal(mu = 0, sigma = 1, n = 2),
    exp(2)
  )
  
  expect_equal(
    moment_lognormal(mu = 2, sigma = 0, n = 1),
    exp(2)
  )
})