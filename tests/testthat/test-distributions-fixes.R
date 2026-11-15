# Regression tests for the distribution-function fixes: log-scale PIG and
# inverse-gamma densities, vectorized random generators, and q-functions.

test_that("PIG and Poisson-inverse-gamma densities stay finite at large counts", {
  expect_true(all(is.finite(dpinvgaus(c(100, 200), mu = 80, eta = 0.5))))
  expect_true(all(is.finite(dpinvgaus(0:300, mu = 2, eta = 0.5))))
  expect_true(all(is.finite(dpinvgamma(c(100, 200), mu = 2, eta = 0.5))))
  y <- 0:3000
  p <- dpinvgaus(y, mu = 80, eta = 0.5)
  expect_equal(sum(p), 1, tolerance = 1e-6)
  expect_equal(sum(y * p), 80, tolerance = 1e-4)
  p2 <- dpinvgaus(y, mu = 5, eta = 0.5, form = "Type 2")
  expect_equal(sum(p2), 1, tolerance = 1e-6)
  p3 <- dpinvgamma(0:2000, mu = 2, eta = 0.3)
  expect_equal(sum(p3), 1, tolerance = 1e-4)
  expect_equal(sum((0:2000) * p3), 2, tolerance = 1e-3)
  expect_equal(dpinvgaus(3, 2, 0.5, log = TRUE), log(dpinvgaus(3, 2, 0.5)))
})

test_that("off-support values give 0 with a warning, bad parameters NaN", {
  expect_warning(p <- dpinvgaus(c(-1, 1.5, 2), 1, 1), "whole number")
  expect_equal(p[1:2], c(0, 0))
  expect_warning(expect_true(is.nan(dpinvgamma(1, mu = -1))))
  expect_warning(expect_true(is.nan(dplind(0, mean = -1, theta = 1))))
})

test_that("dplind uses lambda", {
  expect_equal(dplind(0, lambda = 2, theta = 1), dplind(0, mean = 3, theta = 1))
  expect_equal(dplind(0, lambda = 2, theta = 1), 2 / 9)
  y <- 0:400
  expect_equal(sum(y * dplind(y, mean = 2, theta = 1.5)), 2, tolerance = 1e-8)
})

test_that("q-functions handle p = 1 and agree with the CDF", {
  expect_equal(qplind(1, 2, 1.5), Inf)
  expect_equal(qpinvgaus(1, 2, 0.5), Inf)
  expect_equal(qpinvgamma(1, 2, 0.5), Inf)
  expect_equal(qgwar(1, 2, 2, 4), Inf)
  expect_true(is.nan(qplind(1.5, 2, 1.5)))
  for (pr in c(0.1, 0.5, 0.9)) {
    q <- qplind(pr, 2, 1.5)
    expect_gte(pplind(q, 2, 1.5), pr)
    if (q > 0) expect_lt(pplind(q - 1, 2, 1.5), pr)
    q2 <- qpinvgaus(pr, 2, 0.5)
    expect_gte(ppinvgaus(q2, 2, 0.5), pr - 1e-10)
  }
  expect_equal(qplind(c(0.1, 0.5, 0.9), lambda = 4.67, theta = 7),
               qplind(c(0.1, 0.5, 0.9), mean = 4.67 * 9 / 56, theta = 7))
})

test_that("PLL lambda and mean parameterizations agree", {
  lam <- 2 * 7 * 8 / 9
  expect_equal(pplindLnorm(2, lambda = lam, theta = 7, sigma = 0.5),
               pplindLnorm(2, mean = 2, theta = 7, sigma = 0.5))
})

test_that("r-functions accept vector parameters and reproduce the mean", {
  set.seed(10)
  n <- 40000
  mu <- c(1, 4)
  grp <- rep_len(1:2, n)
  check <- function(y, target, tol = 0.05) {
    expect_length(y, n)
    m <- tapply(y, grp, mean)
    expect_equal(as.numeric(m), target, tolerance = tol)
  }
  check(rpLnorm(n, mean = mu, sigma = 0.5), mu * exp(0.125))
  check(rplind(n, mean = mu, theta = 1.5), mu)
  check(rplindLnorm(n, mean = mu, theta = 1.5, sigma = 0.5), mu * exp(0.125))
  check(rpge(n, mean = mu, shape = 2, scale = 1), mu)
  check(rplindGamma(n, mean = mu, theta = 1.5, alpha = 0.5), mu)
  check(rcom(n, mu = mu, nu = 0.7), mu)
  check(rgp1(n, mu = mu, phi = 0.5), mu)
  check(rgp2(n, mu = mu, alpha = 0.3), mu)
  check(rpinvgaus(n, mu = mu, eta = 0.5), mu)
  expect_warning(y <- rpLnorm(3, mean = c(1, -1, 2)), "NAs produced")
  expect_true(is.na(y[2]))
  expect_false(anyNA(y[-2]))
})

test_that("internal objects are not exported", {
  ex <- getNamespaceExports("flexCountReg")
  expect_false("msg1" %in% ex)
  expect_false(".dpge_one" %in% ex)
})
