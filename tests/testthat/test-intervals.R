# Confidence, prediction and tolerance intervals for countreg() models.

pois_data <- function(n = 40, seed = 1) {
  set.seed(seed)
  data.frame(y = stats::rpois(n, 2), e = stats::runif(n, 0.5, 2))
}

test_that("exact Poisson tolerance limits match poisson.test and qpois", {
  d <- pois_data()
  fit <- countreg(y ~ 1 + offset(log(e)), data = d, family = "POISSON",
                  method = "BFGS")
  new <- data.frame(e = c(0.5, 1, 4))
  S <- sum(d$y)
  E <- sum(d$e)
  for (side in c("two-sided", "upper", "lower")) {
    ti <- tolerance_interval(fit, new, content = 0.9, confidence = 0.95,
                             side = side)
    expect_identical(attr(ti, "tolerance")$method, "exact")
    alt <- switch(side, "two-sided" = "two.sided", upper = "less",
                  lower = "greater")
    ci <- stats::poisson.test(S, E, conf.level = 0.95,
                              alternative = alt)$conf.int
    lo <- switch(side, "two-sided" = stats::qpois(0.05, new$e * ci[1]),
                 upper = rep(0, 3), lower = stats::qpois(0.1, new$e * ci[1]))
    hi <- switch(side, "two-sided" = stats::qpois(0.95, new$e * ci[2]),
                 upper = stats::qpois(0.9, new$e * ci[2]), lower = rep(Inf, 3))
    expect_equal(ti$lower, lo, label = side)
    expect_equal(ti$upper, hi, label = side)
  }
})

test_that("the exact Poisson construction has at least its confidence", {
  # Enumerate the sufficient statistic S ~ Poisson(n lambda)
  conf_of <- function(n, lambda, p, g) {
    s <- 0:stats::qpois(1 - 1e-14, n * lambda)
    ok <- vapply(s, function(si) {
      ci <- stats::poisson.test(si, n, conf.level = g)$conf.int
      L <- stats::qpois((1 - p) / 2, ci[1])
      U <- stats::qpois((1 + p) / 2, ci[2])
      stats::ppois(U, lambda) - stats::ppois(L - 1, lambda) >= p
    }, logical(1))
    sum(stats::dpois(s, n * lambda)[ok])
  }
  for (lam in c(0.1, 1, 10)) for (n in c(5, 40)) {
    expect_gte(conf_of(n, lam, 0.9, 0.95), 0.95)
  }
})

test_that("confidence intervals for the mean agree with the delta method", {
  data("washington_roads", package = "flexCountReg", envir = environment())
  fit <- countreg(Total_crashes ~ lnaadt + lnlength, data = washington_roads,
                  family = "NB2", method = "BFGS")
  scen <- data.frame(lnaadt = log(c(2000, 20000)), lnlength = 0)
  ci <- predict(fit, scen, interval = "confidence")
  expect_named(ci, c("fit", "lwr", "upr"))
  expect_equal(ci$fit, predict(fit, scen))
  X <- cbind(1, scen$lnaadt, scen$lnlength)
  V <- stats::vcov(fit)[1:3, 1:3]
  eta <- as.vector(X %*% stats::coef(fit)[1:3])
  se <- sqrt(diag(X %*% V %*% t(X)))
  expect_equal(ci$lwr, exp(eta - stats::qnorm(0.975) * se), tolerance = 0.01)
  expect_equal(ci$upr, exp(eta + stats::qnorm(0.975) * se), tolerance = 0.01)
  # Halton draws: identical on repeat
  expect_identical(predict(fit, scen, interval = "confidence"), ci)
  expect_identical(attr(ci, "interval")$uncertainty, "normal")
})

test_that("prediction intervals average the PMF over Halton draws", {
  d <- pois_data()
  fit <- countreg(y ~ 1, data = d, family = "POISSON", method = "BFGS")
  pi <- predict(fit, data.frame(z = 1), interval = "prediction",
                level = 0.9, ndraws = 300)
  # Independent calculation: theta* = b + se * qnorm(Halton)
  b <- stats::coef(fit)[[1]]
  se <- sqrt(stats::vcov(fit)[1, 1])
  th <- b + se * stats::qnorm(randtoolbox::halton(300))
  pm <- rowMeans(vapply(th, function(t) stats::dpois(0:100, exp(t)),
                        numeric(101)))
  cum <- cumsum(pm)
  expect_equal(pi$lwr, which(cum >= 0.05)[1] - 1)
  expect_equal(pi$upr, which(cum >= 0.95)[1] - 1)
  expect_true(pi$lwr <= pi$fit && pi$fit <= pi$upr)
})

test_that("tolerance intervals have the requested content on the set", {
  data("washington_roads", package = "flexCountReg", envir = environment())
  fit <- countreg(Total_crashes ~ lnaadt + lnlength, data = washington_roads,
                  family = "NB2", method = "BFGS")
  scen <- data.frame(lnaadt = log(c(2000, 10000, 30000)), lnlength = 0)
  joint <- tolerance_interval(fit, scen, content = 0.9)
  env <- tolerance_interval(fit, scen, content = 0.9, method = "envelope")
  simu <- tolerance_interval(fit, scen, content = 0.9,
                             scope = "simultaneous")
  up <- tolerance_interval(fit, scen, content = 0.9, side = "upper")
  lo <- tolerance_interval(fit, scen, content = 0.9, side = "lower")
  expect_s3_class(joint, "flexCountReg_tolerance")
  expect_true(all(joint$min_content >= 0.9))
  expect_true(all(env$min_content >= 0.9))
  expect_true(all(joint$upper - joint$lower <= env$upper - env$lower))
  expect_true(all(simu$upper - simu$lower >= joint$upper - joint$lower))
  expect_true(all(up$lower == 0) && all(lo$upper == Inf))
  expect_true(all(up$min_content >= 0.9) && all(lo$min_content >= 0.9))
  # Wider than the plug-in quantile interval of the fitted distribution
  mu <- predict(fit, scen)
  size <- 1 / exp(stats::coef(fit)[["ln(alpha)"]])
  plug_w <- stats::qnbinom(0.95, size, mu = mu) -
    stats::qnbinom(0.05, size, mu = mu)
  expect_true(all(joint$upper - joint$lower >= plug_w))
  expect_output(print(joint), "content 0.9, confidence 0.95")
})

test_that("intervals handle zero models, missing rows and errors", {
  data("washington_roads", package = "flexCountReg", envir = environment())
  zi <- countreg(Animal ~ lnaadt + lnlength, data = washington_roads,
                 family = "NB2", zero_formula = ~ speed50, method = "BFGS")
  s2 <- data.frame(lnaadt = c(log(5000), NA), lnlength = 0, speed50 = 0:1)
  ti <- tolerance_interval(zi, s2, content = 0.9)
  expect_true(ti$min_content[1] >= 0.9)
  expect_true(is.na(ti$lower[2]) && ti$status[2] == "missing values")
  pr <- predict(zi, s2, interval = "prediction")
  expect_true(is.na(pr$lwr[2]) && !is.na(pr$lwr[1]))
  expect_error(tolerance_interval(zi, s2, content = 1.2), "between 0 and 1")
  expect_error(tolerance_interval(zi, s2, method = "exact"), "intercept-only")
  expect_error(predict(zi, s2, interval = "prediction", type = "true"),
               "observed count")
})

test_that("the parametric bootstrap is reproducible with a seed", {
  skip_on_cran()
  d <- pois_data(60, seed = 4)
  d$x <- stats::rnorm(60)
  fit <- countreg(y ~ x, data = d, family = "NB2", method = "BFGS")
  scen <- data.frame(x = c(0, 1))
  a <- predict(fit, scen, interval = "prediction", uncertainty = "bootstrap",
               bootstraps = 30, seed = 9)
  b <- predict(fit, scen, interval = "prediction", uncertainty = "bootstrap",
               bootstraps = 30, seed = 9)
  expect_identical(a, b)
  expect_identical(attr(a, "interval")$uncertainty, "bootstrap")
  ti <- tolerance_interval(fit, scen, content = 0.8,
                           uncertainty = "bootstrap", bootstraps = 40,
                           seed = 9)
  expect_true(all(ti$min_content >= 0.8))
  expect_identical(attr(ti, "tolerance")$uncertainty, "bootstrap")
  expect_lte(attr(ti, "tolerance")$accepted, 40)
})

test_that("a non-positive-definite covariance is detected", {
  expect_true(flexCountReg:::.fcr_is_pd(diag(2)))
  expect_false(flexCountReg:::.fcr_is_pd(matrix(c(1, 1, 1, 1), 2)))
  expect_false(flexCountReg:::.fcr_is_pd(matrix(c(1, NA, NA, 1), 2)))
})
