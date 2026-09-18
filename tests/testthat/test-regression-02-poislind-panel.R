# Error 02: Poisson-Lindley panel likelihood, not the RENB likelihood.
# Uses maxLik's public objectiveFn() accessor, with no mocking or copied body.

reg02_data <- function() {
  data.frame(
    id = rep(seq_len(16), each = 3),
    x = rep(c(-1, 0, 1), 16),
    y = c(0,0,1, 1,2,4, 0,1,0, 3,4,8, 1,0,2, 2,4,5,
          0,0,0, 2,1,4, 5,7,9, 0,1,2, 1,3,2, 0,0,1,
          3,2,6, 1,1,3, 0,2,1, 4,5,8)
  )
}

# One real fit, cached only within this test file. No random seed is needed.
reg02_fit <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      cached <<- poisLind.re(y ~ x, group_var = "id", data = reg02_data(),
                            method = "NM", max.iters = 3000, print.level = 0)
    }
    cached
  }
})

reg02_objective <- function() maxLik::objectiveFn(reg02_fit()$model)

# Independent oracle: Lindley = a mixture of Gamma(1,theta) and Gamma(2,theta).
# Panel total has a mixture of negative-binomial distributions; conditional
# allocation across observations is multinomial. This does not repeat the
# package's gamma/log-sum-exp likelihood expression or use its private helpers.
reg02_reference <- function(y, mean, theta, group) {
  lambda <- mean * theta * (theta + 1) / (theta + 2)
  panels <- split(seq_along(y), group, drop = TRUE)
  vapply(panels, function(i) {
    total <- sum(y[i])
    exposure <- sum(lambda[i])
    components <- c(
      log(theta / (theta + 1)) +
        dnbinom(total, size = 1, mu = exposure / theta, log = TRUE),
      -log1p(theta) +
        dnbinom(total, size = 2, mu = 2 * exposure / theta, log = TRUE)
    )
    largest <- max(components)
    mixture <- largest + log(sum(exp(components - largest)))
    mixture + dmultinom(y[i], prob = lambda[i] / exposure, log = TRUE)
  }, numeric(1))
}

test_that("error 02: singleton panels agree with the univariate distribution", {
  objective <- reg02_objective()
  y <- 0:12
  X <- matrix(1, nrow = length(y), ncol = 1)
  for (mu in c(0.5, 2, 8)) {
    for (theta in c(0.5, 3, 10)) {
      actual <- objective(c(log(mu), log(theta)), y, X,
                          group = seq_along(y), offset_vec = rep(0, length(y)))
      expected <- dplind(y, mean = mu, theta = theta, log = TRUE)
      expect_equal(sort(as.numeric(actual)), sort(as.numeric(expected)),
                   tolerance = 1e-10)
    }
  }
  # The original bug's concrete example: y=0, marginal mean=2, theta=3.
  actual <- objective(c(log(2), log(3)), 0, matrix(1), 1, 0)
  expect_equal(exp(as.numeric(actual)), 0.3254437869822485, tolerance = 1e-12)
})

test_that("error 02: multi-observation panels match an independent mixture calculation", {
  objective <- reg02_objective()
  y <- c(0, 2, 4, 1, 3, 0)
  X <- cbind(1, c(-1, 0, 1, -0.5, 0.5, 1.5))
  group <- c(1, 1, 1, 2, 2, 2)
  offset <- c(0, 0.2, -0.1, 0.3, 0, -0.2)
  beta <- c(log(2), 0.25)
  theta <- 3
  actual <- objective(c(beta, log(theta)), y, X, group, offset)
  expected <- reg02_reference(y, exp(drop(X %*% beta) + offset), theta, group)
  expect_length(actual, 2)
  expect_equal(sort(as.numeric(actual)), sort(as.numeric(expected)),
               tolerance = 1e-10)
})

test_that("error 02: zero-count panels retain positive finite probabilities", {
  objective <- reg02_objective()
  y <- rep(0, 8)
  group <- rep(1:2, each = 4)
  actual <- objective(c(log(2), log(3)), y, matrix(1, 8, 1), group, rep(0, 8))
  expected <- reg02_reference(y, rep(2, 8), 3, group)
  expect_true(all(is.finite(actual) & actual <= 0))
  expect_equal(as.numeric(actual), as.numeric(expected), tolerance = 1e-10)
})

test_that("error 02: duplicating independent panels doubles total log likelihood", {
  objective <- reg02_objective()
  y <- c(0, 2, 1, 4)
  X <- cbind(1, c(-1, 0, 0.5, 1))
  group <- c(1, 1, 2, 2)
  offset <- c(0, 0.1, -0.2, 0.3)
  par <- c(log(2), 0.2, log(3))
  original <- objective(par, y, X, group, offset)
  doubled <- objective(par, rep(y, 2), rbind(X, X), c(group, group + 2),
                       rep(offset, 2))
  expect_length(doubled, 4)
  expect_equal(sum(doubled), 2 * sum(original), tolerance = 1e-10)
})

test_that("error 02: likelihood is invariant to row order and group labels", {
  objective <- reg02_objective()
  y <- c(0, 2, 1, 4, 3, 1)
  X <- cbind(1, c(-1, 0, 0.5, 1, -0.5, 0.2))
  group <- c(1, 1, 2, 2, 3, 3)
  offset <- c(0, 0.1, -0.2, 0.3, 0.2, -0.1)
  par <- c(log(2), 0.2, log(3))
  original <- objective(par, y, X, group, offset)
  ix <- c(6, 3, 1, 5, 2, 4)
  shuffled <- objective(par, y[ix], X[ix, , drop = FALSE],
                        c("z", "a", "m")[group[ix]], offset[ix])
  expect_equal(sort(as.numeric(shuffled)), sort(as.numeric(original)),
               tolerance = 1e-10)
})

test_that("error 02: explicit offsets change the marginal mean exactly once", {
  objective <- reg02_objective()
  y <- c(0, 2, 4, 1)
  X <- matrix(1, 4, 1)
  group <- c(1, 1, 2, 2)
  par <- c(log(2), log(3))
  off <- c(log(0.5), log(2), log(1.5), log(3))
  actual <- objective(par, y, X, group, off)
  expected <- reg02_reference(y, 2 * exp(off), 3, group)
  expect_equal(as.numeric(actual), as.numeric(expected), tolerance = 1e-10)
  expect_gt(abs(sum(actual) - sum(objective(par, y, X, group, rep(0, 4)))), 1e-4)
})

test_that("error 02: large counts remain finite on the log scale", {
  objective <- reg02_objective()
  y <- c(200, 300, 450, 600)
  group <- c(1, 1, 2, 2)
  actual <- objective(c(log(250), log(3)), y, matrix(1, 4, 1), group, rep(0, 4))
  expected <- reg02_reference(y, rep(250, 4), 3, group)
  expect_true(all(is.finite(actual)))
  expect_equal(as.numeric(actual), as.numeric(expected), tolerance = 1e-9)
})

test_that("error 02: invalid optimizer proposals produce negative infinity", {
  objective <- reg02_objective()
  for (par in list(c(log(2), 1000), c(log(2), -1000), c(Inf, log(3)))) {
    actual <- objective(par, c(0, 1, 2, 3), matrix(1, 4, 1),
                        c(1, 1, 2, 2), rep(0, 4))
    expect_equal(as.numeric(actual), rep(-Inf, 2))
  }
})

test_that("error 02: the public fit reports the corrected log likelihood", {
  model <- reg02_fit()
  d <- reg02_data()
  expect_s3_class(model, "flexCountReg")
  expect_identical(model$model$modelType, "poisLindRE")
  expect_true(all(is.finite(model$model$estimate)))
  expect_true(is.finite(model$model$theta) && model$model$theta > 0)
  X <- model.matrix(y ~ x, d)
  mean <- exp(drop(X %*% model$model$beta_pred))
  expected <- sum(reg02_reference(d$y, mean, model$model$theta, d$id))
  expect_equal(as.numeric(model$model$maximum), expected, tolerance = 1e-8)
  expect_equal(as.numeric(model$model$LL), expected, tolerance = 1e-8)
})

test_that("error 02: the public fitter passes offsets to the optimizer", {
  d <- reg02_data()
  d$off <- rep(c(-0.2, 0.1, 0.3), 16)
  fit <- poisLind.re(y ~ x, group_var = "id", data = d, offset = "off",
                     method = "NM", max.iters = 3000, print.level = 0)
  X <- model.matrix(y ~ x, d)
  mean <- exp(drop(X %*% fit$model$beta_pred) + d$off)
  expected <- sum(reg02_reference(d$y, mean, fit$model$theta, d$id))
  expect_equal(as.numeric(fit$model$maximum), expected, tolerance = 1e-8)
})
