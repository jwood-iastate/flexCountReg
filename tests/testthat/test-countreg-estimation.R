# Regression tests for countreg() estimation fixes (data handling, offsets,
# standard errors, underreporting sign and input validation).

test_that("formula offset equals offset argument and matches MASS::glm.nb", {
  d <- washington_roads
  a <- countreg(Total_crashes ~ lnaadt + offset(lnlength), d, method = "BFGS")
  b <- countreg(Total_crashes ~ lnaadt, d, offset = "lnlength",
                method = "BFGS")
  m <- MASS::glm.nb(Total_crashes ~ lnaadt + offset(lnlength), d)
  expect_equal(a$model$estimate[1:2], b$model$estimate[1:2],
               tolerance = 1e-6)
  expect_equal(unname(a$model$estimate[1:2]), unname(coef(m)),
               tolerance = 1e-3)
  expect_equal(a$model$maximum, as.numeric(logLik(m)), tolerance = 1e-5)
  expect_equal(predict(a, newdata = d), as.numeric(fitted(m)), tolerance = 1e-3)
})

test_that("giving both a formula offset and the offset argument is an error", {
  expect_error(
    countreg(Total_crashes ~ lnaadt + offset(lnlength), washington_roads,
             offset = "lnlength"),
    "not both")
})

test_that("several offset columns are summed", {
  d <- washington_roads
  d$half <- d$lnlength / 2
  a <- countreg(Total_crashes ~ lnaadt, d, offset = c("half", "half"),
                method = "BFGS")
  b <- countreg(Total_crashes ~ lnaadt, d, offset = "lnlength",
                method = "BFGS")
  expect_equal(a$model$estimate, b$model$estimate, tolerance = 1e-6)
})

test_that("weights as a column or a vector match MASS::glm.nb", {
  d <- washington_roads
  d$w <- rep(1:3, length.out = nrow(d))
  a <- countreg(Total_crashes ~ lnaadt, d, weights = "w", method = "BFGS")
  b <- countreg(Total_crashes ~ lnaadt, d, weights = d$w, method = "BFGS")
  m <- MASS::glm.nb(Total_crashes ~ lnaadt, d, weights = w)
  expect_equal(a$model$estimate, b$model$estimate, tolerance = 1e-6)
  expect_equal(unname(a$model$estimate[1:2]), unname(coef(m)),
               tolerance = 1e-3)
})

test_that("missing values in an auxiliary formula drop the same rows", {
  d <- washington_roads
  d$ShouldWidth04[1:5] <- NA
  m <- suppressMessages(
    countreg(Total_crashes ~ lnaadt, d, dis_param_formula_1 = ~ShouldWidth04,
             method = "BFGS"))
  expect_equal(m$model$nobs, nrow(d) - 5)
  expect_equal(m$model$n_dropped, 5)
  ref <- countreg(Total_crashes ~ lnaadt, d[-(1:5), ],
                  dis_param_formula_1 = ~ShouldWidth04, method = "BFGS")
  expect_equal(m$model$estimate, ref$model$estimate, tolerance = 1e-6)
})

test_that("bootstrapped SEs are stored by parameter name", {
  set.seed(42)
  m <- countreg(Total_crashes ~ lnaadt + lnlength, washington_roads,
                method = "BFGS", stderr = "boot", bootstraps = 20)
  expect_identical(names(m$model$se), names(m$model$estimate))
  expect_identical(rownames(m$model$vcov), names(m$model$estimate))
  expect_equal(m$model$vcov_type, "bootstrap")
  expect_equal(m$model$n_boot_ok, 20)
  h <- sqrt(diag(solve(-m$model$hessian)))
  # Each bootstrap SE is close to its own Hessian SE (not another term's)
  expect_true(all(abs(m$model$se / h - 1) < 0.5))
  s <- utils::capture.output(tab <- summary(m))
  expect_equal(tab$`Std. Err.`, unname(round(m$model$se, 3)))
})

test_that("stderr options are validated and stored", {
  d <- washington_roads
  r <- countreg(Total_crashes ~ lnaadt, d, method = "BFGS", stderr = "Robust")
  expect_equal(r$model$stderr, "robust")
  expect_equal(r$model$vcov_type, "robust")
  expect_equal(unname(r$model$se),
               unname(sqrt(diag(sandwich::sandwich(r$model)))),
               tolerance = 1e-8)
  expect_error(countreg(Total_crashes ~ lnaadt, d, stderr = "bootstrapped"))
  expect_message(
    countreg(Total_crashes ~ lnaadt, d[1:200, ], method = "BFGS",
             stderr = "boot", max.iters = 50),
    "200 bootstrap")
})

test_that("probit and logit underreporting share the same sign convention", {
  set.seed(3)
  n <- 3000
  x <- rnorm(n)
  z <- rbinom(n, 1, 0.5)
  mu <- exp(0.5 + 0.3 * x)
  p_rep <- plogis(-0.5 + 2 * z)
  y <- rbinom(n, rnbinom(n, size = 2, mu = mu), p_rep)
  sim <- data.frame(y, x, z)
  lg <- countreg(y ~ x, sim, underreport_formula = ~z, method = "BFGS")
  pr <- countreg(y ~ x, sim, underreport_formula = ~z,
                 underreport_family = "probit", method = "BFGS")
  b_lg <- lg$model$estimate["Underreporting:z"]
  b_pr <- pr$model$estimate["Underreporting:z"]
  expect_gt(b_lg, 0)
  expect_gt(b_pr, 0)
  expect_error(countreg(y ~ x, sim, underreport_formula = ~z,
                        underreport_family = "probti"))
})

test_that("invalid inputs fail early with clear errors", {
  d <- washington_roads
  expect_error(countreg(Total_crashes ~ lnaadt, d, family = "NB3"),
               "Unknown `family`")
  expect_error(countreg(Total_crashes ~ lnaadt, d, family = "Poisson",
                        dis_param_formula_1 = ~speed50),
               "no distribution parameter")
  expect_error(countreg(Total_crashes ~ lnaadt, d, family = "NB2",
                        dis_param_formula_2 = ~speed50),
               "no second distribution parameter")
  expect_error(countreg(Total_crashes ~ lnaadt, d, start.vals = c(0, 1)),
               "has 3 parameter")
  expect_error(countreg(Total_crashes ~ lnaadt, d,
                        start.vals = c(a = 0, b = 1, c = 0)),
               "do not match")
  expect_error(countreg(Total_crashes ~ lnaadt, d, method = "BFSG"),
               "`method` must be one of")
})

test_that("named start values are matched by name", {
  d <- washington_roads
  sv <- c(`ln(alpha)` = 0, lnaadt = 0.5, `(Intercept)` = -4)
  m <- countreg(Total_crashes ~ lnaadt, d, start.vals = sv, method = "BFGS")
  ref <- countreg(Total_crashes ~ lnaadt, d, method = "BFGS")
  expect_equal(m$model$estimate, ref$model$estimate, tolerance = 1e-3)
})

test_that("convergence is reported", {
  m <- countreg(Total_crashes ~ lnaadt, washington_roads, method = "BFGS")
  expect_true(m$model$converged)
  expect_warning(
    countreg(Total_crashes ~ lnaadt, washington_roads, method = "NM",
             max.iters = 3),
    "did not report convergence")
})
