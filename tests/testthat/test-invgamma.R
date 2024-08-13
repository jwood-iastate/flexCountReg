test_that("Inverse Gamma Tests", {

  ## Density: Integral should approximate 1
  set.seed(666)
  # plot(dinvgamma, -0.2, 6)
  x <- seq(-1, 6, by = 0.01)
  y <- dinvgamma(x)
  expect_lt(abs(sum(y * 0.01) - 1), 0.005)

  ## PDF: function is increasing 
  p_vec <- pinvgamma(c(-1, -.01, 0, 0.3, 0.5, 10, 30), shape = 3, scale = 2)
  expect_true(all(diff(p_vec) >= 0))
    
  # qinvgamma(c(0.1,0.3,0.5,0.9,0.95), shape=3, scale=2)
  # rinvgamma(30, shape=3, scale=2)
  # # expect_equal(2 * 2, 4)
})

# library(testthat)


