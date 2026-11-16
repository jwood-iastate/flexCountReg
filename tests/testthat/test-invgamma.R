test_that("Inverse Gamma Tests", {

  set.seed(666)
  alpha <- 4
  beta <- 3

  # plot(\(x) dinvgamma(x, alpha, beta), -0.2, 4)
  delta <- 0.01
  x <- seq(-1, 4, by = delta)
  y <- suppressWarnings(dinvgamma(x, alpha, beta))
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




test_that("inverse-gamma functions give NaN for invalid parameters", {
  expect_warning(d <- dinvgamma(1, shape = -1), "'shape' must be positive")
  expect_true(is.nan(d))
  expect_warning(pr <- pinvgamma(c(1, 2), scale = c(-1, 1)),
                 "'scale' must be positive")
  expect_true(is.nan(pr[1]))
  expect_equal(pr[2], pinvgamma(2))
  expect_warning(q <- qinvgamma(c(1.2, 0.5), shape = 2), "must be in \\[0, 1\\]")
  expect_true(is.nan(q[1]))
  expect_equal(q[2], qinvgamma(0.5, shape = 2))
  expect_warning(q <- qinvgamma(0.5, shape = 2, scale = 0), "'scale'")
  expect_true(is.nan(q))
  expect_true(is.na(dinvgamma(NA)))
  expect_true(is.na(pinvgamma(NA)))
  expect_true(is.na(qinvgamma(NA)))
  expect_equal(pinvgamma(c(-1, 0, Inf)), c(0, 0, 1))
  expect_equal(qinvgamma(c(0, 1)), c(0, Inf))
})

test_that("rinvgamma follows the base-R conventions for n and parameters", {
  expect_length(rinvgamma(c(9, 9, 9)), 3)
  expect_length(rinvgamma(0), 0)
  expect_error(rinvgamma(-1), "non-negative")
  expect_warning(r <- rinvgamma(3, shape = c(2, -1, 2)), "NAs produced")
  expect_identical(is.na(r), c(FALSE, TRUE, FALSE))
  set.seed(2)
  r <- rinvgamma(20000, shape = c(3, 5), scale = 2)
  expect_equal(mean(r[c(TRUE, FALSE)]), 2 / (3 - 1), tolerance = 0.03)
  expect_equal(mean(r[c(FALSE, TRUE)]), 2 / (5 - 1), tolerance = 0.03)
})
