# Regression tests for corrections from the external (ChatGPT) development
# review, section 1, not covered by test-review-p0-fixes.R.

test_that("quantile residuals are computed on the log scale (1.8)", {
  dd <- data.frame(y = c(rep(0, 24), 25))
  mp <- countreg(y ~ 1, data = dd, family = "Poisson", method = "BFGS")
  r <- .fcr_quantile_residuals(mp, u = rep(0.5, 25))
  expect_equal(r[25], 10.61415007, tolerance = 1e-8)
  expect_equal(r[1], stats::qnorm(0.5 * stats::dpois(0, 1)), tolerance = 1e-8)
  # Agreement with the direct formula where it is accurate (NB2)
  data("washington_roads", package = "flexCountReg")
  m <- countreg(Total_crashes ~ lnaadt, data = washington_roads,
                family = "NB2", method = "BFGS")
  set.seed(2)
  u <- stats::runif(nrow(washington_roads))
  y <- m$model$y
  mu <- predict(m)
  a <- exp(coef(m)[["ln(alpha)"]])
  F1 <- stats::pnbinom(y - 1, size = 1 / a, mu = mu)
  F2 <- stats::pnbinom(y, size = 1 / a, mu = mu)
  expect_equal(.fcr_quantile_residuals(m, u), stats::qnorm(F1 + u * (F2 - F1)),
               tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("Empirical Bayes uses the evaluated response (1.6)", {
  skip_on_cran()
  data("washington_roads", package = "flexCountReg")
  w <- washington_roads
  w$c2 <- 2 * w$Total_crashes
  fit <- function(f) suppressWarnings(countreg_rp(
    f, rpar_formula = ~ -1 + speed50, data = w, group_var = "ID",
    ndraws = 50, method = "BFGS"))
  a <- fit(I(c2 / 2) ~ lnaadt)
  b <- fit(Total_crashes ~ lnaadt)
  expect_equal(as.vector(predict(a, newdata = w, method = "Individual")),
               as.vector(predict(b, newdata = w, method = "Individual")),
               tolerance = 1e-8)
  w2 <- w
  w2$c2 <- NULL
  expect_error(predict(a, newdata = w2, method = "Individual"), "c2")
})

test_that("comparison baselines keep the formula's environment (1.5)", {
  data("washington_roads", package = "flexCountReg")
  d <- washington_roads
  f <- function() {
    mylog <- function(v) log(v)
    fit <- countreg(Total_crashes ~ 0 + mylog(AADT), data = d,
                    family = "NBP", method = "BFGS")
    r <- suppressMessages(suppressWarnings(
      regCompTest(fit, basemodel = "NB2", variables = TRUE)))
    ref <- countreg(Total_crashes ~ 0 + mylog(AADT), data = d,
                    family = "NB2", method = "BFGS")
    c(r$LLbase, as.numeric(logLik(ref)))
  }
  ll <- f()
  expect_equal(ll[1], ll[2], tolerance = 1e-4)
})

test_that("the inverse-Gaussian root is accurate and overflow-safe (1.11)", {
  m <- 1
  v <- 1
  lambda <- 1e-12
  t <- m * v / (2 * lambda)
  expect_equal(m / (1 + t + sqrt(t) * sqrt(t + 2)), 1e-12, tolerance = 1e-6)
  set.seed(1)
  x <- rinvgaus(1000, mean = 1e200, shape = 1)
  expect_true(all(is.finite(x) & x > 0))
})

test_that("the Weibull shape is solved for any coefficient of variation (1.11)", {
  p <- calculate_params(mean_value = 1, sd_value = 0.01)
  expect_equal(p[1], 127.53015, tolerance = 1e-6)
  for (s in c(0.01, 0.5, 3, 20)) {
    p <- calculate_params(mean_value = 2, sd_value = s)
    expect_equal(p[2] * gamma(1 + 1 / p[1]), 2, tolerance = 1e-8)
    expect_equal(p[2] * sqrt(gamma(1 + 2 / p[1]) - gamma(1 + 1 / p[1])^2), s,
                 tolerance = 1e-6)
  }
  expect_error(calculate_params(1, -1), "positive")
})

test_that("regCompTable names unnamed models (1.11)", {
  data("washington_roads", package = "flexCountReg")
  d <- washington_roads
  a <- countreg(Total_crashes ~ lnaadt, data = d, family = "NB2",
                method = "BFGS")
  b <- countreg(Total_crashes ~ lnaadt, data = d, family = "POISSON",
                method = "BFGS")
  expect_identical(names(regCompTable(list(a, b), RSquare = FALSE)),
                   c("Parameter", "Model 1", "Model 2"))
  expect_identical(names(regCompTable(list(x = a, b), RSquare = FALSE)),
                   c("Parameter", "x", "Model 2"))
})

test_that("the family validity hook catches GP2's coupled domain", {
  expect_identical(.fcr_family_valid("GP2", c(5, 5, 1), c(0.5, 1.5, 0.5)),
                   c(FALSE, TRUE, TRUE))
  expect_identical(.fcr_family_valid("NB2", c(1, NA, Inf), 1),
                   c(TRUE, FALSE, FALSE))
  # A valid underdispersed GP2 is a proper distribution with finite support
  expect_equal(sum(dgp2(0:100, 1, -0.5)), 1, tolerance = 1e-12)
})

test_that("intervals reject invalid draws and explain unresolved tails (1.4, 1.11)", {
  data("washington_roads", package = "flexCountReg")
  d <- washington_roads
  set.seed(2)
  d$z <- rnorm(nrow(d))
  m <- countreg(Total_crashes ~ lnaadt, data = d, family = "NB2",
                dis_param_formula_1 = ~ z, method = "BFGS")
  nd <- d[1:3, ]
  nd$z[2] <- NA
  pi <- expect_no_warning(predict(m, newdata = nd, interval = "prediction",
                                  ndraws = 50))
  expect_true(is.na(pi$lwr[2]) && is.na(pi$upr[2]))
  expect_false(anyNA(pi[c(1, 3), c("lwr", "upr")]))
  expect_equal(attr(pi, "interval")$accepted_by_row, c(50, 0, 50))
  # A max_count below the upper limit gives an explained NA
  expect_warning(r <- predict(m, newdata = d[1:2, ], interval = "prediction",
                              ndraws = 20, max_count = 2), "max_count")
  expect_true(all(is.na(r$upr)))

  # GP2 draws outside the coupled domain are rejected and counted
  set.seed(3)
  n <- 300
  x <- runif(n)
  dg <- data.frame(y = rbinom(n, 6, plogis(0.5 + x)), x = x)
  g <- countreg(y ~ x, data = dg, family = "GP2", method = "BFGS")
  r <- suppressWarnings(predict(g, newdata = data.frame(x = c(0, 1, 30)),
                                interval = "prediction", ndraws = 100))
  acc <- attr(r, "interval")$accepted_by_row
  expect_equal(acc[1:2], c(100, 100))
  expect_equal(acc[3], 0)
  expect_true(is.na(r$upr[3]))
})

test_that("tolerance intervals handle a degenerate predictor direction (1.11)", {
  set.seed(1)
  dd <- data.frame(x = runif(200, 0, 2))
  dd$y <- rpois(200, exp(0.8 * dd$x))
  m0 <- countreg(y ~ 0 + x, data = dd, family = "POISSON", method = "BFGS")
  nd <- data.frame(x = c(0, 1))
  ti <- tolerance_interval(m0, nd, content = 0.9)
  expect_identical(ti$status, c("ok", "ok"))
  expect_equal(ti$fit[1], 1)
  expect_gte(ti$min_content[1], 0.9)
  tb <- tolerance_interval(m0, nd, content = 0.9, uncertainty = "bootstrap",
                           bootstraps = 30, seed = 1)
  expect_identical(tb$status, c("ok", "ok"))
  expect_equal(.fcr_mahalanobis_psd(matrix(c(1, 2, 0, 0), 2), c(0, 0),
                                    diag(c(1, 0))), c(1, 4))
  expect_equal(.fcr_mahalanobis_psd(matrix(1:4, 2), c(0, 0),
                                    matrix(0, 2, 2)), c(0, 0))
})
