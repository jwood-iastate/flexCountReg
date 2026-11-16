# halton_dists(), corr_haltons() and cor2cov(): draws have the requested
# moments, ranges and correlations.

test_that("halton_dists() draws have the stated mean, spread and range", {
  n <- 20000
  x <- halton_dists("n", mean = 2, sdev = 0.5, ndraws = n)
  expect_length(x, n)
  expect_equal(c(mean(x), sd(x)), c(2, 0.5), tolerance = 1e-3)
  u <- halton_dists("u", mean = 8, sdev = 3, ndraws = n)
  expect_equal(range(u), c(5, 11), tolerance = 1e-3)  # half-width 3
  expect_equal(sd(u), 3 / sqrt(3), tolerance = 1e-3)
  tr <- halton_dists("t", mean = 1, sdev = 2, ndraws = n)
  expect_equal(range(tr), c(-1, 3), tolerance = 1e-2)
  expect_equal(mean(tr), 1, tolerance = 1e-3)
  g <- halton_dists("g", mean = 2, sdev = 1, ndraws = n)
  expect_equal(c(mean(g), sd(g)), c(2, 1), tolerance = 5e-3)
  ln <- halton_dists("ln", mean = 0, sdev = 0.5, ndraws = n)
  expect_equal(mean(log(ln)), 0, tolerance = 1e-3)
  h <- randtoolbox::halton(10)
  expect_equal(halton_dists("n", 0, 1, hdraw = h), qnorm(h))
  expect_error(halton_dists("q", 0, 1), "Invalid distribution type")
})

test_that("corr_haltons() reproduces means, SDs and correlations", {
  C <- matrix(c(1, 0.5, -0.3, 0.5, 1, 0.2, -0.3, 0.2, 1), 3)
  s <- c(1, 2, 0.5)
  m <- c(0, 5, -1)
  x <- corr_haltons(m, stdev = s, correlations = C, ndraws = 5000)
  expect_equal(dim(x), c(5000L, 3L))
  expect_equal(colMeans(x), m, tolerance = 0.02)
  expect_equal(apply(x, 2, sd), s, tolerance = 0.02)
  expect_equal(cor(x), C, tolerance = 0.02, ignore_attr = TRUE)
  # The Cholesky route gives the same draws
  y <- corr_haltons(m, cholesky = chol(cor2cov(C, s)), ndraws = 5000)
  expect_equal(x, y)
  # Truncated normal: within the bounds, still positively correlated
  tn <- corr_haltons(c(1, 1), stdev = c(1, 1),
                     correlations = matrix(c(1, 0.6, 0.6, 1), 2),
                     dist = "truncated_normal", lower = 0, upper = 3,
                     ndraws = 5000)
  expect_true(all(tn >= 0 & tn <= 3))
  expect_gt(cor(tn)[1, 2], 0.4)
  # Each column uses its own mean: truncated at 0 from a mean of 1 vs 2
  tn2 <- corr_haltons(c(1, 2), stdev = c(1, 1),
                      correlations = diag(2), dist = "truncated_normal",
                      lower = 0, ndraws = 5000)
  expect_lt(mean(tn2[, 1]), mean(tn2[, 2]))
  expect_equal(mean(tn2[, 1]),
               1 + dnorm(-1) / (1 - pnorm(-1)), tolerance = 0.01)
  hd <- qnorm(randtoolbox::halton(100, 2))
  expect_equal(corr_haltons(c(0, 0), stdev = c(1, 1), correlations = diag(2),
                            hdraws = hd), hd, ignore_attr = TRUE)
  expect_error(corr_haltons(c(0, 0), stdev = c(1, 1), correlations = diag(2),
                            hdraws = hd[, 1, drop = FALSE]), "columns")
  expect_error(corr_haltons(c(0, 0)), "Cholesky")
  expect_error(corr_haltons(c(0, 0), stdev = c(1, 1), correlations = diag(2),
                            dist = "gamma"), "Unsupported")
})
