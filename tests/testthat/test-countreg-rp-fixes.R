# Regression tests for countreg.rp(): draw conventions, rp_mgf(), Exact vs
# Simulated predictions, heterogeneity in means, input validation and
# compatibility with fits from earlier versions.

rp_data <- function() {
  d <- washington_roads
  d$x <- d$lnaadt - mean(d$lnaadt)
  d
}

test_that("rp_mgf matches simulation from generate_random_draws", {
  set.seed(5)
  h <- matrix(runif(2e5), ncol = 1)
  for (dist in c("n", "u", "t", "g", "ln")) {
    for (x in c(-2, 0.7, 2.5)) {
      if (dist == "ln" && (x > 0 || !requireNamespace("lamW", quietly = TRUE)))
        next
      # gamma at x = 2.5 has infinite Monte Carlo variance (2 x s^2 / mu > 1)
      if (dist == "g" && x > 2) next
      mu <- if (dist == "g") 0.6 else 0.3
      sim <- mean(exp(generate_random_draws(h, mu, 0.4, dist, "x",
                                            matrix(x))$xb_rand_mat))
      expect_equal(rp_mgf(dist, x, mu, 0.4), sim, tolerance = 0.01)
    }
  }
  expect_equal(rp_mgf("ln", 1, 0, 0.5), Inf)
  expect_equal(rp_mgf("g", 10, 0.5, 1), Inf)
  expect_equal(rp_mgf("u", 0, 0.3, 0.4), 1)
})

test_that("uniform draws span mu +/- sigma", {
  h <- matrix(c(0, 0.5, 1), ncol = 1)
  b <- generate_random_draws(h, 1, 0.5, "u", "x", matrix(1))$xb_rand_mat
  expect_equal(as.vector(b), c(0.5, 1, 1.5))
})

test_that("Exact and Simulated agree with heterogeneity, for each distribution", {
  d <- rp_data()
  for (dist in c("n", "u", "t", "g")) {
    m <- suppressWarnings(countreg.rp(
      Total_crashes ~ lnlength, rpar_formula = ~ -1 + x, data = d,
      rpardists = c(x = dist), ndraws = 100, method = "BFGS",
      het_mean_formula = ~ speed50, het_var_formula = ~ ShouldWidth04))
    m$model$ndraws <- 5000
    ps <- predict(m, newdata = d, method = "Simulated")
    pe <- predict(m, newdata = d, method = "Exact")
    expect_lt(max(abs(ps / pe - 1)), 0.01)
  }
})

test_that("heterogeneity in means is additive with one effect per parameter", {
  d <- rp_data()
  m <- countreg.rp(Total_crashes ~ lnlength, rpar_formula = ~ -1 + x + speed50,
                   data = d, rpardists = c(x = "n", speed50 = "n"),
                   het_mean_formula = ~ ShouldWidth04, ndraws = 50,
                   method = "BFGS", max.iters = 200)
  nm <- names(m$model$estimate)
  expect_true(all(c("HetMean:x:ShouldWidth04", "HetMean:speed50:ShouldWidth04")
                  %in% nm))
  # Additive form, checked through the exact prediction for the normal case
  cf <- m$model$estimate
  nd <- d[1:5, ]
  shift <- function(k) cf[paste0("HetMean:", k, ":ShouldWidth04")] *
    nd$ShouldWidth04
  manual <- exp(cf["(Intercept)"] + cf["lnlength"] * nd$lnlength +
                  nd$x * (cf["x:Mean"] + shift("x")) +
                  0.5 * nd$x^2 * cf["x:St.Dev"]^2 +
                  nd$speed50 * (cf["speed50:Mean"] + shift("speed50")) +
                  0.5 * nd$speed50^2 * cf["speed50:St.Dev"]^2)
  expect_equal(predict(m, newdata = nd, method = "Exact"), unname(manual))
})

test_that("rpardists are matched by name and validated", {
  d <- rp_data()
  m <- countreg.rp(Total_crashes ~ lnlength, rpar_formula = ~ -1 + x + speed50,
                   data = d, rpardists = c(speed50 = "u", x = "n"),
                   ndraws = 30, method = "BFGS", max.iters = 50)
  expect_equal(m$model$rpardists, c(x = "n", speed50 = "u"))
  expect_error(countreg.rp(Total_crashes ~ lnlength,
                           rpar_formula = ~ -1 + x + speed50, data = d,
                           rpardists = c(speed50 = "n")), "Missing: x")
  expect_error(countreg.rp(Total_crashes ~ lnlength, rpar_formula = ~ -1 + x,
                           data = d, rpardists = c(x = "q")),
               "Unknown random-parameter distribution")
})

test_that("invalid specifications fail early", {
  d <- rp_data()
  expect_error(countreg.rp(Total_crashes ~ lnlength, rpar_formula = ~ x,
                           data = d, rpardists = c(Intercept = "n", x = "n")),
               "intercept is in both")
  expect_error(countreg.rp(Total_crashes ~ lnlength, rpar_formula = ~ -1 + x,
                           data = d, group_var = "nope"), "not a column")
  expect_error(countreg.rp(Total_crashes ~ lnlength, rpar_formula = ~ -1 + x,
                           data = d, family = "Poisson",
                           dis_param_formula_1 = ~ speed50),
               "no distribution parameter")
  expect_error(countreg.rp(Total_crashes ~ lnlength, rpar_formula = ~ -1 + x,
                           data = d, method = "BFSG"), "`method` must be")
})

test_that("an implicit intercept in rpar_formula is not a random parameter", {
  d <- rp_data()
  m <- countreg.rp(Total_crashes ~ lnlength, rpar_formula = ~ x, data = d,
                   ndraws = 30, method = "BFGS", max.iters = 50)
  expect_equal(m$model$rpar, "x")
  expect_false("(Intercept):Mean" %in% names(m$model$estimate))
})

test_that("gamma random parameters start from a positive mean", {
  d <- rp_data()
  d$negx <- -d$speed50   # Poisson start for this coefficient is negative
  m <- suppressWarnings(countreg.rp(
    Total_crashes ~ lnaadt + lnlength, rpar_formula = ~ -1 + negx, data = d,
    rpardists = c(negx = "g"), ndraws = 50, method = "BFGS", max.iters = 50))
  expect_true(is.finite(m$model$maximum))
})

test_that("formula offsets, missing values and weights are handled", {
  d <- rp_data()
  a <- countreg.rp(Total_crashes ~ lnaadt + offset(lnlength),
                   rpar_formula = ~ -1 + speed50, data = d, ndraws = 50,
                   method = "BFGS")
  b <- countreg.rp(Total_crashes ~ lnaadt, rpar_formula = ~ -1 + speed50,
                   data = d, offset = "lnlength", ndraws = 50, method = "BFGS")
  expect_equal(a$model$estimate, b$model$estimate, tolerance = 1e-6)

  d2 <- d
  d2$speed50[1:4] <- NA
  m <- suppressMessages(countreg.rp(Total_crashes ~ lnaadt,
                                    rpar_formula = ~ -1 + speed50, data = d2,
                                    ndraws = 30, method = "BFGS",
                                    max.iters = 50))
  expect_equal(m$model$nobs, nrow(d) - 4)

  # Doubling every weight doubles the log-likelihood at the same estimates
  d$w2 <- 2
  w <- countreg.rp(Total_crashes ~ lnaadt, rpar_formula = ~ -1 + speed50,
                   data = d, weights = "w2", group_var = "ID", ndraws = 30,
                   method = "BFGS",
                   start.vals = b$model$estimate[names(b$model$estimate)],
                   max.iters = 0)
  u <- countreg.rp(Total_crashes ~ lnaadt, rpar_formula = ~ -1 + speed50,
                   data = d, group_var = "ID", ndraws = 30, method = "BFGS",
                   start.vals = b$model$estimate, max.iters = 0)
  expect_equal(w$model$maximum, 2 * u$model$maximum)
})

test_that("fits from earlier versions predict with their own conventions", {
  d <- rp_data()
  m <- countreg.rp(Total_crashes ~ lnlength, rpar_formula = ~ -1 + x,
                   data = d, rpardists = c(x = "u"), ndraws = 50,
                   method = "BFGS")
  old <- m
  old$model$draw_version <- NULL
  old$model$rpar <- NULL
  old$model$estimate["x:St.Dev"] <- 2 * m$model$estimate["x:St.Dev"]
  expect_equal(suppressMessages(predict(old, newdata = d, method = "Exact")),
               predict(m, newdata = d, method = "Exact"))
})

test_that("lognormal Exact predictions with positive x are Inf with a warning", {
  d <- rp_data()
  m <- suppressWarnings(countreg.rp(
    Total_crashes ~ lnlength, rpar_formula = ~ -1 + speed50, data = d,
    rpardists = c(speed50 = "ln"), ndraws = 30, method = "BFGS",
    max.iters = 30))
  expect_warning(p <- predict(m, newdata = d, method = "Exact"), "infinite")
  expect_true(all(is.infinite(p[d$speed50 > 0])))
  expect_true(all(is.finite(p[d$speed50 == 0])))
})
