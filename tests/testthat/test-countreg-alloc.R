# Total-plus-allocation models: countreg_alloc()

sim_alloc <- function(n = 1500, seed = 42) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n), v = rbinom(n, 1, 0.4),
                  len = runif(n, 0.5, 2), id = rep(seq_len(n / 3), each = 3),
                  yr = rep(2016:2018, n / 3))
  mu <- exp(0.6 + 0.5 * d$x + log(d$len))
  tot <- stats::rnbinom(n, size = 2, mu = mu)
  psi1 <- -0.3 + 0.4 * d$v
  psi2 <- psi1 + exp(0.2 - 0.3 * d$v)
  xb <- 0.7 * d$z
  pr <- cbind(stats::pnorm(psi1 - xb),
              stats::pnorm(psi2 - xb) - stats::pnorm(psi1 - xb),
              1 - stats::pnorm(psi2 - xb))
  Y <- t(vapply(seq_len(n), function(i)
    as.numeric(stats::rmultinom(1, tot[i], pr[i, ])), numeric(3)))
  d$a <- Y[, 1]
  d$b <- Y[, 2]
  d$c <- Y[, 3]
  d
}

test_that("the ordered model recovers the truth and matches polr and countreg", {
  d <- sim_alloc()
  m <- countreg_alloc(cbind(a, b, c) ~ x + offset(log(len)), data = d,
                      family = "NB2", alloc_formula = ~ z,
                      threshold_formula = ~ v)
  expect_s3_class(m, "countreg_alloc")
  expect_true(m$converged)
  cf <- coef(m)
  expect_equal(unname(cf[c("x", "Allocation:z", "Threshold1:v")]),
               c(0.5, 0.7, 0.4), tolerance = 0.15)
  expect_equal(m$categories, c("a", "b", "c"))
  # The parts are separable without random parameters: the total is the
  # countreg() fit, the split the ordered probit of polr()
  cr <- countreg(I(a + b + c) ~ x + offset(log(len)), data = d,
                 family = "NB2", method = "BFGS")
  expect_equal(unname(cf[1:3]), unname(coef(cr)), tolerance = 1e-4)
  skip_if_not_installed("MASS")
  m0 <- countreg_alloc(cbind(a, b, c) ~ x + offset(log(len)), data = d,
                       family = "NB2", alloc_formula = ~ z)
  n <- nrow(d)
  long <- data.frame(
    cat = factor(rep(c("a", "b", "c"), each = n), ordered = TRUE),
    z = rep(d$z, 3), wt = c(d$a, d$b, d$c))
  long <- long[long$wt > 0, ]
  pm <- suppressMessages(MASS::polr(cat ~ z, data = long, weights = wt,
                                    method = "probit", Hess = TRUE))
  c0 <- coef(m0)
  expect_equal(unname(c0["Allocation:z"]), unname(coef(pm)),
               tolerance = 1e-4)
  expect_equal(unname(c(c0["Threshold1:(Intercept)"],
                        c0["Threshold1:(Intercept)"] +
                          exp(c0["Threshold2:(Intercept)"]))),
               unname(pm$zeta), tolerance = 1e-4)
  expect_equal(sqrt(vcov(m0)["Allocation:z", "Allocation:z"]),
               unname(sqrt(diag(vcov(pm)))[1]), tolerance = 1e-3)
  Y <- as.matrix(d[, c("a", "b", "c")])
  const <- sum(lgamma(rowSums(Y) + 1) - rowSums(lgamma(Y + 1)))
  expect_equal(as.numeric(logLik(m0)),
               as.numeric(logLik(pm) + logLik(cr)) + const,
               tolerance = 1e-6)
})

test_that("the analytic gradient matches numeric differences", {
  d <- sim_alloc(600)
  for (link in c("ordered_probit", "ordered_logit", "mnl")) {
    thr <- if (link == "mnl") ~ 1 else ~ v
    fit_at <- function(p) {
      suppressWarnings(suppressMessages(countreg_alloc(
        cbind(a, b, c) ~ x, data = d, family = "NB2", alloc_formula = ~ z,
        threshold_formula = thr, link = link, start.vals = p,
        max.iters = 0)))
    }
    m <- countreg_alloc(cbind(a, b, c) ~ x, data = d, family = "NB2",
                        alloc_formula = ~ z, threshold_formula = thr,
                        link = link)
    p <- coef(m) + 0.05
    num <- numDeriv::grad(function(q) fit_at(q)$loglik, p)
    ana <- colSums(fit_at(p)$gradientObs)
    expect_equal(ana, num, tolerance = 1e-5, ignore_attr = TRUE)
  }
  # Truncated total
  dp <- d[d$a + d$b + d$c > 0, ]
  fit_t <- function(p) {
    suppressWarnings(countreg_alloc(
      cbind(a, b, c) ~ x, data = dp, family = "NB2", alloc_formula = ~ z,
      truncated = TRUE, start.vals = p, max.iters = 0))
  }
  p <- coef(countreg_alloc(cbind(a, b, c) ~ x, data = dp, family = "NB2",
                           alloc_formula = ~ z, truncated = TRUE)) + 0.05
  expect_equal(colSums(fit_t(p)$gradientObs),
               numDeriv::grad(function(q) fit_t(q)$loglik, p),
               tolerance = 1e-5, ignore_attr = TRUE)
})

test_that("a Poisson total with mnl is independent Poissons; NB2 is negmultinom", {
  set.seed(7)
  n <- 1200
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  tot <- stats::rpois(n, exp(0.3 + 0.5 * d$x))
  eta <- cbind(0, 0.5 + 0.8 * d$z, -0.5 - 0.4 * d$z)
  pr <- exp(eta) / rowSums(exp(eta))
  Y <- t(vapply(seq_len(n), function(i)
    as.numeric(stats::rmultinom(1, tot[i], pr[i, ])), numeric(3)))
  d$a <- Y[, 1]
  d$b <- Y[, 2]
  d$c <- Y[, 3]
  m <- countreg_alloc(cbind(a, b, c) ~ x, data = d, family = "POISSON",
                      alloc_formula = ~ z, link = "mnl")
  expect_equal(as.numeric(logLik(m)),
               sum(stats::dpois(Y, predict(m), log = TRUE)),
               tolerance = 1e-10)
  expect_equal(names(coef(m))[3:4],
               c("Allocation[b]:(Intercept)", "Allocation[b]:z"))
  skip_if_not_installed("nnet")
  keep <- rowSums(Y) > 0
  mn <- nnet::multinom(Y[keep, ] ~ z, data = d[keep, ], trace = FALSE)
  expect_equal(unname(coef(m)[c(3, 5, 4, 6)]), c(coef(mn)),
               tolerance = 1e-4)
  # NB2 total with no covariates: the negative multinomial
  tot2 <- stats::rnbinom(n, size = 2, mu = 3)
  Y2 <- t(vapply(seq_len(n), function(i)
    as.numeric(stats::rmultinom(1, tot2[i], c(0.5, 0.3, 0.2))), numeric(3)))
  d2 <- data.frame(a = Y2[, 1], b = Y2[, 2], c = Y2[, 3])
  m2 <- countreg_alloc(cbind(a, b, c) ~ 1, data = d2, family = "NB2",
                       link = "mnl")
  long <- data.frame(y = c(Y2), cat = factor(rep(c("a", "b", "c"), each = n)),
                     unit = rep(seq_len(n), 3))
  nm <- negmultinom(y ~ cat, group_var = "unit", data = long)
  expect_equal(as.numeric(logLik(m2)), as.numeric(logLik(nm)),
               tolerance = 1e-6)
})

test_that("zero totals, truncation and the fractional likelihood", {
  d <- sim_alloc()
  # Zero totals are kept (not truncated)
  m <- countreg_alloc(cbind(a, b, c) ~ x, data = d, family = "NB2",
                      alloc_formula = ~ z)
  expect_equal(m$nobs, nrow(d))
  expect_true(any(m$total == 0))
  dp <- d[d$a + d$b + d$c > 0, ]
  expect_message(countreg_alloc(cbind(a, b, c) ~ x, data = dp,
                                family = "NB2", alloc_formula = ~ z),
                 "truncated = TRUE")
  expect_error(countreg_alloc(cbind(a, b, c) ~ x, data = d, family = "NB2",
                              truncated = TRUE), "zero total")
  mt <- countreg_alloc(cbind(a, b, c) ~ x, data = dp, family = "NB2",
                       alloc_formula = ~ z, truncated = TRUE)
  # The truncated fit estimates the untruncated total's parameters
  expect_equal(unname(coef(mt)["x"]), unname(coef(m)["x"]), tolerance = 0.1)
  nd <- dp[1:4, ]
  tp <- predict(mt, newdata = nd, type = "total_positive")
  expect_true(all(tp > predict(mt, newdata = nd, type = "total")))
  expect_equal(rowSums(predict(mt, newdata = nd, conditional = TRUE)), tp,
               tolerance = 1e-10, ignore_attr = TRUE)
  # Fractional: sandwich standard errors, no likelihood
  mf <- countreg_alloc(cbind(a, b, c) ~ x, data = d, family = "NB2",
                       alloc_formula = ~ z, likelihood = "fractional",
                       group_var = "id")
  expect_equal(mf$stderr, "cluster")
  expect_error(AIC(mf), "quasi-likelihood")
  expect_true(is.na(glance(mf)$AIC))
  expect_warning(countreg_alloc(cbind(a, b, c) ~ x, data = d,
                                family = "NB2", likelihood = "fractional",
                                stderr = "normal"), "not valid")
})

test_that("predictions, simulation, methods and tv()", {
  d <- sim_alloc()
  m <- countreg_alloc(cbind(a, b, c) ~ x + offset(log(len)), data = d,
                      family = "NB2", alloc_formula = ~ tv(z, time = yr),
                      threshold_formula = ~ v)
  pr <- predict(m, type = "prob")
  expect_equal(dim(pr), c(nrow(d), 3))
  expect_equal(unname(rowSums(pr)), rep(1, nrow(d)), tolerance = 1e-12)
  expect_equal(predict(m), fitted(m))
  expect_equal(unname(rowSums(predict(m))), predict(m, type = "total"),
               tolerance = 1e-10)
  expect_equal(predict(m, newdata = d[5:7, ]), fitted(m)[5:7, ],
               ignore_attr = TRUE)
  s <- simulate(m, nsim = 2, seed = 1)
  expect_length(s, 2)
  expect_equal(dim(s[[1]]), c(nrow(d), 3))
  expect_equal(colnames(s[[1]]), c("a", "b", "c"))
  te <- temporal_effects(m)
  expect_equal(unique(te$part), "Allocation")
  expect_equal(nrow(te), 3)
  td <- tidy(m, conf.int = TRUE)
  expect_setequal(unique(td$component),
                  c("total", "distribution", "allocation", "thresholds"))
  expect_equal(glance(m)$npar, length(coef(m)))
  expect_output(print(summary(m)), "zero total")
  expect_equal(confint(m)[, 1] < coef(m), rep(TRUE, length(coef(m))),
               ignore_attr = TRUE)
})

test_that("input errors", {
  d <- sim_alloc(300)
  expect_error(countreg_alloc(a ~ x, data = d), "two or more")
  expect_error(countreg_alloc(cbind(a, b, c) ~ x, data = d, link = "mnl",
                              threshold_formula = ~ v), "ordered links")
  expect_error(countreg_alloc(cbind(a, b, c) ~ x, data = d,
                              threshold_formula = ~ 0 + v), "intercept")
  expect_error(countreg_alloc(cbind(a, b, c) ~ x, data = d,
                              stderr = "cluster"), "group_var")
  d$a[1] <- 0.5
  expect_error(countreg_alloc(cbind(a, b, c) ~ x, data = d), "whole numbers")
})
