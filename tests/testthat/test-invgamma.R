test_that("Inverse Gamma Tests", {

  set.seed(666)
  alpha <- 4
  beta <- 3

  # plot(\(x) dinvgamma(x, alpha, beta), -0.2, 4)
  delta <- 0.01
  x <- seq(-1, 4, by = delta)
  y <- dinvgamma(x, alpha, beta)
  riemann_sum <- delta * sum(y)
  expect_equal(riemann_sum, 1, tolerance = 0.01)


  ## Density: mode should be close to theoretical
  mode_ig <- beta / (alpha + 1)
  mode_riemann <- x[which.max(y)]
  expect_equal(mode_riemann, mode_ig, tolerance = 0.1)


  ## PDF: function is increasing
  qq <- c(-1, -.01, 0, 0.3, 0.5, 10, 30)
  p_vec <- pinvgamma(qq, shape = 3, scale = 2)
  expect_true(all(diff(p_vec) >= 0))

  ## Quantile: check if the quantiles of the probs. above equal the original
  ## quantile vector (ignoring negative values)
  rr <- qinvgamma(p_vec[3:7], shape = 3, scale = 2)
  expect_equal(rr, qq[3:7])

  ## Random samples:
  ## Check that sample mean and variance approximate the theoretical ones
  mean_ig <- beta / (alpha - 1)
  var_ig <- beta^2 / ((alpha - 1)^2 * (alpha - 2))

  ig_samples <- rinvgamma(1e5, alpha, beta)
  sample_mean <- mean(ig_samples)
  sample_var <- var(ig_samples)

  expect_equal(sample_mean, mean_ig, tolerance = 0.1)
  expect_equal(sample_var, var_ig, tolerance = 0.1)
})

# library(testthat)



