# Regression tests for Error 06: Poisson-Generalized-Exponential integration.
# The fixed Halton vector makes these tests deterministic and inexpensive.

test_that("dpge returns one probability per observation", {
  h <- c(0.25, 0.50, 0.75)

  p_scalar <- dpge(
    0,
    mean = 2,
    shape = 2,
    scale = 1,
    haltons = h
  )
  expect_length(p_scalar, 1)

  p_vector <- dpge(
    0:3,
    mean = 2,
    shape = 2,
    scale = 1,
    haltons = h
  )
  expect_length(p_vector, 4)
  expect_true(all(is.finite(p_vector)))
})

test_that("dpge averages all integration draws for each observation", {
  h <- c(0.25, 0.50, 0.75)
  q_unit <- -log(-expm1(log(h) / 2))
  mu_draw <- 2 * q_unit / (digamma(3) - digamma(1))

  expected <- vapply(0:3, function(y) {
    mean(stats::dpois(y, mu_draw))
  }, numeric(1))

  actual <- dpge(
    0:3,
    mean = 2,
    shape = 2,
    scale = 1,
    haltons = h
  )

  expect_equal(actual, expected, tolerance = 1e-12)
  expect_equal(actual[1], 0.21997184, tolerance = 1e-8)
})

test_that("dpge vectorizes observation parameters but not halton draws", {
  h <- c(0.25, 0.50, 0.75)
  means <- c(1, 2, 3, 4)
  shapes <- c(1, 2, 3, 4)
  scales <- c(0.5, 1, 2, 3)

  actual <- dpge(
    0:3,
    mean = means,
    shape = shapes,
    scale = scales,
    haltons = h
  )

  expected <- vapply(seq_along(means), function(i) {
    q_unit <- -log(-expm1(log(h) / shapes[i]))
    mu_draw <- means[i] * q_unit /
      (digamma(shapes[i] + 1) - digamma(1))
    mean(stats::dpois(i - 1, mu_draw))
  }, numeric(1))

  expect_length(actual, length(means))
  expect_equal(actual, expected, tolerance = 1e-12)
})

test_that("dpge log probabilities correspond to ordinary probabilities", {
  h <- c(0.25, 0.50, 0.75)
  p <- dpge(0:5, mean = 2, shape = 2, scale = 1, haltons = h)
  lp <- dpge(0:5, mean = 2, shape = 2, scale = 1,
             haltons = h, log = TRUE)

  expect_equal(exp(lp), p, tolerance = 1e-12)
  expect_true(all(is.finite(lp)))
})

test_that("dpge rejects invalid integration draws", {
  expect_error(dpge(0, mean = 2, shape = 2, scale = 1,
                    haltons = c(0, 0.5, 0.75)))
  expect_error(dpge(0, mean = 2, shape = 2, scale = 1,
                    haltons = c(0.25, 1, 0.75)))
  expect_error(dpge(0, mean = 2, shape = 2, scale = 1,
                    haltons = numeric()))
})

test_that("dpge rejects nonpositive parameters", {
  expect_error(dpge(0, mean = 0, shape = 2, scale = 1,
                    haltons = c(0.25, 0.5, 0.75)))
  expect_error(dpge(0, mean = 2, shape = 0, scale = 1,
                    haltons = c(0.25, 0.5, 0.75)))
  expect_error(dpge(0, mean = 2, shape = 2, scale = 0,
                    haltons = c(0.25, 0.5, 0.75)))
})
