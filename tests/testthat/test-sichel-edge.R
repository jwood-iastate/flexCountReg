# Edge cases of the Sichel distribution functions.

test_that("Sichel approaches the Poisson distribution as sigma -> 0", {
  # 1/sigma > 700 uses the exponentially scaled Bessel function
  expect_equal(dsichel(0:6, mu = 3, sigma = 1e-4, gamma = 1), dpois(0:6, 3),
               tolerance = 1e-3)
  expect_equal(psichel(0:6, mu = 3, sigma = 1e-4, gamma = -0.5),
               ppois(0:6, 3), tolerance = 1e-3)
})

test_that("Sichel functions handle invalid and missing inputs", {
  expect_equal(dsichel(c(-1, 1.5, Inf), mu = 2), c(0, 0, 0))
  expect_true(is.na(dsichel(NA, mu = 2)))
  expect_warning(d <- dsichel(1, mu = -1), "'mu' must be positive")
  expect_true(is.nan(d))
  expect_warning(d <- dsichel(1, sigma = 0), "'sigma' must be positive")
  expect_true(is.nan(d))
  expect_warning(pr <- psichel(c(-1, 2), mu = c(1, -1)), "'mu' must be")
  expect_identical(pr[1], 0)
  expect_true(is.nan(pr[2]))
  expect_warning(psichel(1, sigma = -1), "'sigma' must be positive")
  expect_true(is.na(psichel(NA, mu = 2)))
  expect_warning(
    d <- dsichel(1, mu = 2, gamma = 1e6),
    "Bessel function computation failed")
  expect_equal(d, 0)
})

test_that("qsichel handles the boundaries, tails and log scale", {
  expect_identical(qsichel(c(0, 1, NA), mu = 2), c(0, Inf, NA))
  expect_warning(q <- qsichel(1.2, mu = 2), "must be in \\[0, 1\\]")
  expect_true(is.nan(q))
  expect_warning(q <- qsichel(0.5, sigma = -1), "'sigma' must be positive")
  expect_true(is.nan(q))
  expect_warning(qsichel(0.5, mu = 0), "'mu' must be positive")
  p <- c(0.1, 0.5, 0.9)
  q <- qsichel(p, mu = 2, sigma = 1, gamma = -0.5)
  expect_identical(qsichel(1 - p, mu = 2, sigma = 1, gamma = -0.5,
                           lower.tail = FALSE), q)
  expect_identical(qsichel(log(p), mu = 2, sigma = 1, gamma = -0.5,
                           log.p = TRUE), q)
  pr <- psichel(0:4, mu = 2, sigma = 1, gamma = -0.5)
  expect_equal(psichel(0:4, mu = 2, sigma = 1, gamma = -0.5,
                       lower.tail = FALSE), 1 - pr)
  expect_equal(psichel(0:4, mu = 2, sigma = 1, gamma = -0.5, log.p = TRUE),
               log(pr))
})

test_that("rsichel validates its inputs", {
  expect_length(rsichel(0), 0)
  expect_length(rsichel(c(5, 5, 5), mu = 2), 3)
  expect_error(rsichel(-1), "non-negative integer")
  expect_warning(rsichel(2, mu = -1), "'mu' must be positive")
  expect_warning(rsichel(2, sigma = -1), "'sigma' must be positive")
})
