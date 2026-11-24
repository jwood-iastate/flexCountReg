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
  skip_if_not_installed("numDeriv")
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
  expect_output(print(summary(m)), "unordered; base: a")
  # Another base category: the same fit, coefficients relative to it
  mb <- countreg_alloc(cbind(b, a, c) ~ x, data = d, family = "POISSON",
                       alloc_formula = ~ z, link = "mnl")
  expect_equal(mb$loglik, m$loglik, tolerance = 1e-8)
  expect_equal(unname(coef(mb)[c("Allocation[a]:z", "Allocation[c]:z")]),
               unname(c(-coef(m)["Allocation[b]:z"],
                        coef(m)["Allocation[c]:z"] -
                          coef(m)["Allocation[b]:z"])), tolerance = 1e-4)
  expect_equal(predict(mb)[, c("a", "b", "c")], predict(m), tolerance = 1e-4)
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

# ---- Random parameters ------------------------------------------------------

sim_alloc_rp <- function(G = 200, Tn = 4, rho = 0.5, seed = 1) {
  set.seed(seed)
  n <- G * Tn
  d <- data.frame(site = rep(seq_len(G), each = Tn), x = rnorm(n),
                  z = rnorm(n))
  S <- matrix(c(0.25, rho * 0.3, rho * 0.3, 0.36), 2)
  b <- matrix(rnorm(2 * G), G, 2) %*% chol(S)
  tot <- stats::rnbinom(n, size = 3, mu = exp(0.8 + 0.5 * d$x + b[d$site, 1]))
  xb <- 0.6 * d$z + b[d$site, 2]
  pr <- cbind(stats::pnorm(-0.2 - xb),
              stats::pnorm(0.8 - xb) - stats::pnorm(-0.2 - xb),
              1 - stats::pnorm(0.8 - xb))
  Y <- t(vapply(seq_len(n), function(i)
    as.numeric(stats::rmultinom(1, tot[i], pr[i, ])), numeric(3)))
  d$a <- Y[, 1]
  d$b <- Y[, 2]
  d$c <- Y[, 3]
  d
}

test_that("the simulated likelihood has the right analytic gradient", {
  skip_if_not_installed("numDeriv")
  d <- sim_alloc_rp(G = 40, Tn = 3)
  cases <- list(
    list(rpar_total = ~ 1, rpar_alloc = ~ 1, correlation = "across"),
    list(rpar_total = ~ 1 + x, rpar_alloc = ~ 1, correlation = "within"),
    list(rpar_alloc = ~ 1 + z, rpar_threshold = ~ 1,
         correlation = "across"),
    list(rpar_total = ~ 1, rpar_alloc = ~ 1, link = "mnl",
         correlation = "across"),
    list(rpar_total = ~ 1, rpar_alloc = ~ 1, likelihood = "fractional",
         correlation = "across", stderr = "robust"),
    list(rpar_total = ~ 1, truncated = FALSE, family = "POISSON"))
  for (cs in cases) {
    args <- utils::modifyList(list(
      formula = cbind(a, b, c) ~ x, data = d, family = "NB2",
      alloc_formula = ~ z, group_var = "site", rpar_draws = 15, seed = 3,
      .setup_only = TRUE), cs)
    E <- do.call(countreg_alloc, args)
    set.seed(2)
    p <- flexCountReg:::.fcr_alloc_start(E) +
      stats::runif(length(E$S$names), -0.1, 0.1)
    f <- function(q) sum(flexCountReg:::.fcr_alloc_sim_ll(q, E)$ll)
    ana <- colSums(flexCountReg:::.fcr_alloc_sim_ll(p, E, grad = TRUE)$G)
    expect_equal(ana, numDeriv::grad(f, p), tolerance = 1e-6,
                 info = paste(names(cs), collapse = ","))
  }
})

test_that("the simulated likelihood approaches the exact integral", {
  d <- sim_alloc_rp(G = 150)
  m <- countreg_alloc(cbind(a, b, c) ~ x, data = d, family = "NB2",
                      alloc_formula = ~ z, group_var = "site",
                      rpar_alloc = ~ 1, rpar_draws = 300, seed = 3)
  p <- coef(m)
  E <- countreg_alloc(cbind(a, b, c) ~ x, data = d, family = "NB2",
                      alloc_formula = ~ z, .setup_only = TRUE)
  k <- length(p)
  sdv <- abs(p[k])
  gh <- flexCountReg:::.fcr_gh_rule(60)
  LLq <- vapply(seq_along(gh$z), function(q) {
    ll <- flexCountReg:::.fcr_alloc_rowparts(
      p[-k], E, list(xb = rep(sdv * gh$z[q], nrow(d))), grad = FALSE)$ll
    rowsum(ll, d$site)[, 1]
  }, numeric(150))
  mx <- apply(LLq, 1, max)
  exact <- sum(mx + log(colSums(t(exp(LLq - mx)) * gh$w)))
  expect_equal(m$loglik, exact, tolerance = 0.2 / abs(exact))
})

test_that("correlated random parameters are recovered", {
  skip_on_cran()
  d <- sim_alloc_rp(G = 400, Tn = 5, seed = 11)
  m <- countreg_alloc(cbind(a, b, c) ~ x, data = d, family = "NB2",
                      alloc_formula = ~ z, group_var = "site",
                      rpar_total = ~ 1, rpar_alloc = ~ 1,
                      correlation = "across", rpar_draws = 200, seed = 5)
  expect_true(m$converged)
  vc <- flexCountReg:::.fcr_alloc_varcorr(m)
  expect_equal(sqrt(diag(vc)), c(0.5, 0.6), tolerance = 0.2,
               ignore_attr = TRUE)
  expect_equal(stats::cov2cor(vc)[2, 1], 0.5, tolerance = 0.4)
  expect_equal(unname(coef(m)[c("x", "Allocation:z")]), c(0.5, 0.6),
               tolerance = 0.1)
  expect_equal(m$n_units, 400)
  expect_output(print(summary(m)), "standard deviations and correlations")
  expect_true("random" %in% tidy(m)$component)
  # Predictions: marginal for new units, posterior for estimated ones
  pm <- predict(m, newdata = d[1:6, ])
  pp <- predict(m, newdata = d[1:6, ], random = "posterior")
  expect_equal(dim(pp), c(6, 3))
  expect_false(isTRUE(all.equal(pm, pp)))
  nd <- d[1:2, ]
  nd$site <- -1
  expect_equal(predict(m, newdata = nd, random = "posterior"),
               predict(m, newdata = nd))
  expect_error(predict(m, newdata = d[1:2, c("x", "z")],
                       random = "posterior"), "grouping column")
  s <- simulate(m, seed = 1)[[1]]
  expect_equal(colMeans(s), colMeans(d[, c("a", "b", "c")]),
               tolerance = 0.15, ignore_attr = TRUE)
})

test_that("random-parameter input errors", {
  d <- sim_alloc_rp(G = 20)
  expect_error(countreg_alloc(cbind(a, b, c) ~ x, data = d, link = "mnl",
                              rpar_threshold = ~ 1), "ordered links")
  expect_error(countreg_alloc(cbind(a, b) ~ x, data = d,
                              rpar_threshold = ~ 1), "three or more")
  expect_error(countreg_alloc(cbind(a, b, c) ~ x, data = d,
                              rpar_total = ~ offset(x)), "offset")
})

# ---- Shared site-year term ---------------------------------------------------

sim_alloc_shared <- function(n = 3000, lc = 0.6, la = -0.5, seed = 21) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n), s = rbinom(n, 1, 0.5),
                  site = rep(seq_len(n / 3), each = 3))
  u <- rnorm(n)
  tot <- stats::rpois(n, exp(0.6 + 0.5 * d$x + lc * u))
  xb <- 0.5 * d$z + la * u
  pr <- cbind(stats::pnorm(-xb), stats::pnorm(0.8 - xb) - stats::pnorm(-xb),
              1 - stats::pnorm(0.8 - xb))
  Y <- t(vapply(seq_len(n), function(i)
    as.numeric(stats::rmultinom(1, tot[i], pr[i, ])), numeric(3)))
  d$a <- Y[, 1]
  d$b <- Y[, 2]
  d$c <- Y[, 3]
  d
}

test_that("the shared term has the right analytic gradient", {
  skip_if_not_installed("numDeriv")
  d <- sim_alloc_shared(n = 300)
  d$wt <- rep(1:3, length.out = nrow(d))
  pos <- d[d$a + d$b + d$c > 0, ]
  cases <- list(
    list(shared_formula = ~ 1),
    list(shared_formula = ~ s, shared_in = c("total", "alloc", "threshold"),
         family = "NB2"),
    list(shared_formula = ~ 1, group_var = "site", rpar_total = ~ 1,
         rpar_alloc = ~ 1, correlation = "across"),
    list(shared_formula = ~ s, data = pos, truncated = TRUE),
    list(shared_formula = ~ 1, data = pos, truncated = TRUE,
         group_var = "site", rpar_alloc = ~ 1),
    list(shared_formula = ~ 1, link = "mnl"),
    list(shared_formula = ~ 1, likelihood = "fractional", weights = "wt",
         shared_in = c("alloc", "threshold"), stderr = "robust"))
  for (cs in cases) {
    # (data frames are lists: modifyList() would merge them by column)
    dat <- if (is.null(cs$data)) d else cs$data
    cs$data <- NULL
    args <- utils::modifyList(list(
      formula = cbind(a, b, c) ~ x, family = "POISSON",
      alloc_formula = ~ z, rpar_draws = 10, seed = 3, shared_nodes = 10,
      .setup_only = TRUE), cs)
    args$data <- dat
    E <- suppressMessages(do.call(countreg_alloc, args))
    set.seed(2)
    p <- flexCountReg:::.fcr_alloc_start(E) +
      stats::runif(length(E$S$names), -0.1, 0.1)
    llf <- if (is.null(E$RS)) flexCountReg:::.fcr_alloc_ll else
      flexCountReg:::.fcr_alloc_sim_ll
    f <- function(q) sum(llf(q, E)$ll)
    ana <- colSums(llf(p, E, grad = TRUE)$G)
    expect_equal(ana, numDeriv::grad(f, p), tolerance = 1e-6,
                 info = paste(names(cs), collapse = ","))
  }
})

test_that("a shared term in a Poisson total is the Poisson-lognormal total", {
  d <- sim_alloc_shared(n = 900, la = 0)
  m1 <- countreg_alloc(cbind(a, b, c) ~ x, data = d, family = "POISSON",
                       alloc_formula = ~ z, shared_formula = ~ 1,
                       shared_in = "total")
  m2 <- countreg_alloc(cbind(a, b, c) ~ x, data = d, family = "PLN",
                       alloc_formula = ~ z)
  expect_equal(m1$loglik, m2$loglik, tolerance = 1e-6)
  expect_equal(unname(coef(m1)[1:5]), unname(coef(m2)[c(1:2, 4:6)]),
               tolerance = 1e-4)
  expect_equal(abs(unname(coef(m1)["shared:Total:(Intercept)"])),
               exp(unname(coef(m2)["ln(sigma)"])), tolerance = 1e-4)
  expect_lt(abs(m1$shared$check), 1e-3)
})

test_that("the shared term is recovered, also from truncated data", {
  skip_on_cran()
  d <- sim_alloc_shared(n = 3000, lc = 0.6, la = -0.5)
  m <- countreg_alloc(cbind(a, b, c) ~ x, data = d, family = "POISSON",
                      alloc_formula = ~ z, shared_formula = ~ 1)
  expect_true(m$converged)
  sh <- coef(m)[c("shared:Total:(Intercept)",
                  "shared:Allocation:(Intercept)")]
  # Only the relative sign of the loadings is identified
  expect_equal(unname(sh * sign(sh[1])), c(0.6, -0.5), tolerance = 0.1)
  expect_equal(unname(coef(m)[c("x", "Allocation:z")]), c(0.5, 0.5),
               tolerance = 0.1)
  m0 <- countreg_alloc(cbind(a, b, c) ~ x, data = d, family = "POISSON",
                       alloc_formula = ~ z)
  expect_gt(m$loglik - m0$loglik, 50)
  expect_output(print(summary(m)), "Shared site-year term")
  expect_true("shared" %in% tidy(m)$component)
  # Predictions average over the term: with a positive loading in the
  # total, the expected total is exp(x'b + lambda^2 / 2)
  pt <- predict(m, newdata = d[1:5, ], type = "total")
  lam <- unname(sh[1])
  expect_equal(pt, exp(coef(m)[1] + coef(m)[2] * d$x[1:5] + lam^2 / 2),
               tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(rowSums(predict(m, newdata = d[1:5, ], type = "prob")),
               rep(1, 5))
  expect_error(predict(m, newdata = d[1:5, ], random = "posterior"),
               "unit-level random parameters")
  s <- simulate(m, seed = 1)[[1]]
  expect_equal(colMeans(s), colMeans(d[, c("a", "b", "c")]),
               tolerance = 0.1, ignore_attr = TRUE)
  # Truncated data: the rows with a positive total
  pos <- sim_alloc_shared(n = 6000, lc = 0.6, la = 0.5, seed = 12)
  pos <- pos[pos$a + pos$b + pos$c > 0, ]
  mt <- countreg_alloc(cbind(a, b, c) ~ x, data = pos, family = "POISSON",
                       alloc_formula = ~ z, truncated = TRUE,
                       shared_formula = ~ 1)
  sh <- coef(mt)[c("shared:Total:(Intercept)",
                   "shared:Allocation:(Intercept)")]
  expect_equal(unname(sh * sign(sh[1])), c(0.6, 0.5), tolerance = 0.1)
  expect_equal(unname(coef(mt)[1:2]), c(0.6, 0.5), tolerance = 0.1)
})

test_that("shared-term input errors", {
  d <- sim_alloc_shared(n = 90)
  expect_error(countreg_alloc(cbind(a, b, c) ~ x, data = d, link = "mnl",
                              shared_formula = ~ 1,
                              shared_in = "threshold"), "ordered links")
  expect_error(countreg_alloc(cbind(a, b) ~ x, data = d,
                              shared_formula = ~ 1,
                              shared_in = "threshold"), "three or more")
  expect_error(countreg_alloc(cbind(a, b, c) ~ x, data = d,
                              shared_formula = ~ 1, shared_nodes = 1),
               "shared_nodes")
  expect_error(countreg_alloc(cbind(a, b, c) ~ x, data = d,
                              shared_formula = ~ offset(x)), "offset")
  expect_error(countreg_alloc(cbind(a, b, c) ~ x, data = d,
                              shared_formula = ~ 1, shared_in = "count"))
})
