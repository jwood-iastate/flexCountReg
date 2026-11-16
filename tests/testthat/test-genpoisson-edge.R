# Boundaries, tails and invalid inputs of the generalized Poisson functions.

test_that("GP1 and GP2 reduce to the Poisson distribution at zero dispersion", {
  expect_equal(dgp1(0:6, mu = 2.5, phi = 0), dpois(0:6, 2.5))
  expect_equal(pgp1(0:6, mu = 2.5, phi = 0), ppois(0:6, 2.5))
  expect_equal(qgp1(c(0.1, 0.5, 0.9), mu = 2.5, phi = 0),
               qpois(c(0.1, 0.5, 0.9), 2.5))
  expect_equal(dgp2(0:6, mu = 2.5, alpha = 0), dpois(0:6, 2.5))
  expect_equal(pgp2(0:6, mu = 2.5, alpha = 0), ppois(0:6, 2.5))
  expect_equal(qgp2(c(0.1, 0.5, 0.9), mu = 2.5, alpha = 0),
               qpois(c(0.1, 0.5, 0.9), 2.5))
  set.seed(1)
  expect_equal(mean(rgp1(5000, mu = 2.5, phi = 0)), 2.5, tolerance = 0.05)
  expect_equal(mean(rgp2(5000, mu = 2.5, alpha = 0)), 2.5, tolerance = 0.05)
})

test_that("GP2 equals GP1 with phi = alpha * mu", {
  for (alpha in c(-0.1, 0.3)) {
    expect_equal(dgp2(0:8, mu = 3, alpha = alpha),
                 dgp1(0:8, mu = 3, phi = alpha * 3), tolerance = 1e-10)
    expect_equal(pgp2(0:8, mu = 3, alpha = alpha),
                 pgp1(0:8, mu = 3, phi = alpha * 3), tolerance = 1e-10)
  }
})

test_that("GP quantiles handle the boundaries of the support", {
  # Underdispersed: finite support, p = 1 gives its upper end
  ymax <- qgp2(1, mu = 3, alpha = -0.1)
  expect_true(is.finite(ymax))
  expect_equal(pgp2(ymax, mu = 3, alpha = -0.1), 1)
  expect_equal(pgp2(ymax + 5, mu = 3, alpha = -0.1), 1)
  expect_identical(qgp2(c(0, 1, NA), mu = 3, alpha = 0.2), c(0, Inf, NA))
  expect_identical(qgp1(c(0, 1, NA), mu = 3, phi = 0.5), c(0, Inf, NA))
  expect_true(is.nan(qgp2(1.5, mu = 3, alpha = 0.1)))
  expect_true(is.nan(qgp1(1.5, mu = 3, phi = 0.1)))
  expect_true(is.na(pgp2(NA, mu = 3, alpha = 0.2)))
  expect_equal(pgp2(-1, mu = 3, alpha = 0.2), 0)
  expect_equal(pgp1(-1, mu = 3, phi = 0.2), 0)
  expect_equal(dgp2(c(-1, 1.5), mu = 3, alpha = -0.1), c(0, 0))
})

test_that("GP tails and log scale are consistent", {
  p <- c(0.2, 0.6, 0.95)
  for (alpha in c(-0.1, 0.2)) {
    q <- qgp2(p, mu = 3, alpha = alpha)
    expect_identical(qgp2(1 - p, mu = 3, alpha = alpha, lower.tail = FALSE), q)
    expect_identical(qgp2(log(p), mu = 3, alpha = alpha, log.p = TRUE), q)
    pr <- pgp2(0:5, mu = 3, alpha = alpha)
    expect_equal(pgp2(0:5, mu = 3, alpha = alpha, lower.tail = FALSE,
                      log.p = TRUE), log(1 - pr))
  }
  for (phi in c(-0.2, 0.5)) {
    q <- qgp1(p, mu = 3, phi = phi)
    expect_identical(qgp1(1 - p, mu = 3, phi = phi, lower.tail = FALSE), q)
    expect_identical(qgp1(log(p), mu = 3, phi = phi, log.p = TRUE), q)
  }
})

test_that("invalid GP parameters warn and give missing values", {
  expect_warning(d <- dgp2(1, mu = 3, alpha = -0.5), "1 \\+ alpha \\* mu")
  expect_true(is.na(d))
  expect_warning(q <- qgp2(0.5, mu = -1, alpha = 0.1), "mu > 0")
  expect_true(is.na(q))
  expect_warning(pgp2(1, mu = 3, alpha = Inf), "Invalid")
  expect_warning(qgp1(0.5, mu = -1, phi = 0.1), "`mu` must be greater")
  expect_warning(qgp1(0.5, mu = 1, phi = -2), "`phi` must be greater")
  expect_warning(r <- rgp2(3, mu = c(1, -1, 1), alpha = 0.1), "invalid")
  expect_identical(is.na(r), c(FALSE, TRUE, FALSE))
  expect_warning(r <- rgp1(2, mu = 1, phi = c(0.1, -3)), "phi")
  expect_identical(is.na(r), c(FALSE, TRUE))
})

test_that("GP1 d and p functions give NaN for invalid parameters", {
  expect_warning(d <- dgp1(1, mu = -1, phi = 0.1), "`mu`")
  expect_true(is.nan(d))
  expect_warning(pr <- pgp1(1, mu = 1, phi = -2), "`phi`")
  expect_true(is.nan(pr))
})
