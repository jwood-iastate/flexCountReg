# Numerical paths of the Poisson-Lindley-Gamma functions: the GSL closed
# form against quadrature, extreme means, and argument checks.

test_that("the GSL closed form and quadrature give the same PMF", {
  x <- 0:8
  fast <- dplindGamma(x, mean = 2, theta = 1.5, alpha = 0.4)
  quad <- dplindGamma(x, mean = 2, theta = 1.5, alpha = 0.4,
                      method = "quadrature")
  expect_equal(fast, quad, tolerance = 1e-8)
})

test_that("a failing GSL evaluation falls back to quadrature", {
  quad <- dplindGamma(0:4, mean = 2, theta = 1.5, alpha = 0.4,
                      method = "quadrature")
  local_mocked_bindings(.plg_gsl_U = function(a, b, z) stop("GSL error"))
  expect_equal(dplindGamma(0:4, mean = 2, theta = 1.5, alpha = 0.4), quad,
               tolerance = 1e-10)
  local_mocked_bindings(.plg_gsl_U = function(a, b, z) {
    warning("GSL warning")
    list(val = c(1, 1), err = c(0, 0), status = c(0, 0))
  })
  expect_equal(dplindGamma(0:4, mean = 2, theta = 1.5, alpha = 0.4), quad,
               tolerance = 1e-10)
  local_mocked_bindings(.plg_gsl_U = function(a, b, z) {
    list(val = c(1, 1), err = c(0, 0), status = c(1, 0))
  })
  expect_equal(dplindGamma(0:4, mean = 2, theta = 1.5, alpha = 0.4), quad,
               tolerance = 1e-10)
  local_mocked_bindings(.plg_gsl_U = function(a, b, z) {
    list(val = c(1, 1), err = c(0.5, 0.5), status = c(0, 0))
  })
  expect_equal(dplindGamma(0:4, mean = 2, theta = 1.5, alpha = 0.4), quad,
               tolerance = 1e-10)
})

test_that("tiny means: P(Y = 0) -> 1 and P(Y >= 1) is about the mean", {
  for (alpha in c(0, 0.4)) {
    lp <- dplindGamma(0:1, mean = 1e-310, theta = 1.5, alpha = alpha,
                      log = TRUE)
    expect_equal(lp[1], 0, tolerance = 1e-12)
    expect_equal(lp[2], log(1e-310), tolerance = 1e-6)
    lt <- pplindGamma(0, mean = 1e-310, theta = 1.5, alpha = alpha,
                      lower.tail = FALSE, log.p = TRUE)
    expect_equal(lt, log(1e-310), tolerance = 1e-6)
  }
})

test_that("huge means give a flat PMF over small counts", {
  for (alpha in c(0, 0.4)) {
    lp <- dplindGamma(c(0, 5, 50), mean = 1e305, theta = 1.5, alpha = alpha,
                      log = TRUE)
    expect_true(all(is.finite(lp)))
    expect_lt(diff(range(lp)), 1e-6)
  }
  # Beyond what can be resolved, the CDF warns and returns NaN
  expect_warning(pr <- pplindGamma(5, mean = 1e305, theta = 1.5, alpha = 0.4),
                 "numerical evaluation failed")
  expect_true(is.nan(pr))
})

test_that("the CDF tails are complementary and handle the support ends", {
  q <- c(0, 2, 10)
  lo <- pplindGamma(q, 2, 1.5, 0.4)
  hi <- pplindGamma(q, 2, 1.5, 0.4, lower.tail = FALSE)
  expect_equal(lo + hi, rep(1, 3), tolerance = 1e-8)
  expect_equal(pplindGamma(q, 2, 1.5, 0.4, log.p = TRUE), log(lo),
               tolerance = 1e-8)
  expect_equal(pplindGamma(c(-1, Inf), 2, 1.5, 0.4), c(0, 1))
  expect_equal(pplindGamma(c(-1, Inf), 2, 1.5, 0.4, lower.tail = FALSE),
               c(1, 0))
  expect_equal(pplindGamma(3, mean = 0, theta = 1.5, alpha = 0.4), 1)
})

test_that("PLG functions check their arguments", {
  expect_identical(dplindGamma(numeric(0)), numeric(0))
  expect_error(dplindGamma(1, rel.tol = 0.1), "'rel.tol'")
  expect_error(dplindGamma(1, subdivisions = 10), "'subdivisions'")
  expect_error(dplindGamma(1, log = NA), "'log' must be TRUE or FALSE")
  expect_error(pplindGamma(1, lower.tail = "yes"), "'lower.tail'")
  expect_error(dplindGamma("1"), "must be numeric")
  expect_warning(d <- dplindGamma(1, mean = -1), "Invalid distribution")
  expect_true(is.nan(d))
  expect_warning(d <- dplindGamma(1.5, mean = 2), "Noninteger")
  expect_equal(d, 0)
  expect_true(is.nan(dplindGamma(NaN)))
  expect_true(is.na(dplindGamma(NA_real_)))
  expect_equal(dplindGamma(c(-1, Inf), 2), c(0, 0))
  expect_named(dplindGamma(c(a = 0, b = 1), 2), c("a", "b"))
  expect_warning(r <- rplindGamma(2, mean = c(1, -1)), "NAs produced")
  expect_identical(is.na(r), c(FALSE, TRUE))
})
