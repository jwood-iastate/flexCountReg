# Regression tests for the follow-up review (tickets R01-R12) and panel
# frequency weights.

# ---- R10a: random-parameter distribution codes --------------------------------

test_that("R10a: invalid scalar random distributions are rejected", {
  f <- getFromNamespace(".fcr_match_rpardists", "flexCountReg")
  expect_error(f("invalid", "x", FALSE), "Unknown random-parameter")
  expect_error(f(c(x = "invalid"), "x", FALSE), "Unknown random-parameter")
  expect_error(f(NA_character_, "x", FALSE), "without missing")
  expect_error(f(character(0), "x", FALSE), "without missing")
  expect_error(f(1, "x", FALSE), "character vector")
  expect_equal(f("n", c("x", "z"), FALSE), c(x = "n", z = "n"))
  # Named values are reordered to the random parameters
  expect_equal(f(c(z = "ln", x = "g"), c("x", "z"), FALSE),
               c(x = "g", z = "ln"))
  expect_error(f(c(x = "n"), c("x", "z"), FALSE), "Missing: z")
  expect_error(f(c(x = "n", z = "n", w = "n"), c("x", "z"), FALSE),
               "Unknown: w")
  expect_error(f(c(x = "n", x = "u"), "x", FALSE), "unique")
  # Correlated models use normal random parameters
  expect_equal(f(NULL, c("x", "z"), TRUE), c(x = "n", z = "n"))
})

# ---- R10b: bootstrap acceptance ----------------------------------------------------

test_that("R10b: unknown convergence does not yield accepted bootstrap draws", {
  boot <- getFromNamespace(".fcr_bootstrap", "flexCountReg")
  draws <- NULL
  expect_warning({
    draws <- boot(n = 3, B = 3, coef_names = "b", refit = function(idx) {
      structure(c(b = 1), converged = NA)
    })
  }, "Fewer than 2")
  expect_true(all(is.na(draws)))
  expect_equal(as.vector(attr(draws, "failures")["not converged"]), 3)
})

test_that("R10b: each kind of failed replicate is rejected with its reason", {
  chk <- getFromNamespace(".fcr_boot_check_replicate", "flexCountReg")
  ok <- structure(c(a = 1, b = 2), converged = TRUE)
  expect_identical(chk(ok, c("a", "b")), "")
  expect_identical(chk(ok, c("b", "a")), "")            # order is free
  expect_identical(chk(NULL, "a"), "error")
  expect_identical(chk(c(a = 1, b = 2), c("a", "b")), "not converged")
  expect_identical(chk(structure(c(a = 1, b = 2), converged = FALSE),
                       c("a", "b")), "not converged")
  expect_identical(chk(structure(c(a = 1, b = 2), converged = NA),
                       c("a", "b")), "not converged")
  expect_identical(chk(structure(c(a = 1, a = 2), converged = TRUE),
                       c("a", "b")), "wrong coefficients")
  expect_identical(chk(structure(c(a = 1), converged = TRUE), c("a", "b")),
                   "wrong coefficients")
  expect_identical(chk(structure(c(a = 1, c = 2), converged = TRUE),
                       c("a", "b")), "wrong coefficients")
  expect_identical(chk(structure(c(1, 2), converged = TRUE), c("a", "b")),
                   "wrong coefficients")
  expect_identical(chk(structure(c(a = 1, b = 2, c = 3), converged = TRUE),
                       c("a", "b")), "wrong coefficients")
  expect_identical(chk(structure(c(a = 1, b = NaN), converged = TRUE),
                       c("a", "b")), "non-finite estimates")
  # Refit results carry an explicit logical status
  rr <- getFromNamespace(".fcr_refit_result", "flexCountReg")
  expect_false(attr(rr(c(a = 1), NA), "converged"))
  expect_false(attr(rr(c(a = 1), NULL), "converged"))
  expect_true(attr(rr(c(a = 1), TRUE), "converged"))
})

test_that("R10b: covariance uses the same accepted replicates for every entry", {
  boot <- getFromNamespace(".fcr_bootstrap", "flexCountReg")
  bcov <- getFromNamespace(".fcr_boot_cov", "flexCountReg")
  set.seed(3)
  k <- 0
  draws <- suppressWarnings(boot(n = 5, B = 20, coef_names = c("a", "b"),
    refit = function(idx) {
      k <<- k + 1
      if (k %% 4 == 0) stop("refit failed")
      structure(c(a = mean(idx), b = k %% 3), converged = k %% 5 != 0)
    }))
  ok <- stats::complete.cases(draws)
  expect_equal(sum(ok), 20 - 5 - 4 + 1)   # failures at k %% 4 and k %% 5
  expect_equal(bcov(draws), stats::cov(draws[ok, ]), ignore_attr = TRUE)
  f <- attr(draws, "failures")
  expect_equal(as.vector(f[c("error", "not converged")]), c(5, 3))
})

# ---- R10c: COM numerical controls -------------------------------------------------

test_that("R10c: the COM rate solver honours its tolerance and iteration limit", {
  diag_fun <- getFromNamespace("find_lambda_diag_cpp", "flexCountReg")
  find_fun <- getFromNamespace("find_lambda_cpp", "flexCountReg")
  expect_fun <- getFromNamespace("com_expect_cpp", "flexCountReg")
  # A case away from the Poisson shortcut (nu != 1)
  mu <- 7.3; nu <- 0.45
  r <- diag_fun(mu, nu, tol = 1e-12, max_iter = 500)
  expect_identical(r$status, "converged")
  expect_lte(abs(r$residual), 1e-12 * mu)
  expect_equal(expect_fun(r$lambda, nu), mu, tolerance = 1e-12)
  # Too few iterations: reported as not converged, NaN with a warning
  r1 <- diag_fun(mu, nu, tol = 1e-12, max_iter = 1)
  expect_identical(r1$status, "iteration limit")
  expect_true(is.nan(r1$lambda))
  expect_gt(abs(r1$residual), 1e-12 * mu)
  expect_warning(v <- find_fun(mu, nu, tol = 1e-12, max_iter = 1),
                 "did not reach the mean tolerance")
  expect_true(is.nan(v))
  # A looser tolerance needs fewer iterations and still meets its target
  r2 <- diag_fun(mu, nu, tol = 1e-4, max_iter = 500)
  expect_identical(r2$status, "converged")
  expect_lte(r2$iterations, r$iterations)
  expect_lte(abs(r2$residual), 1e-4 * mu)
  expect_error(find_fun(mu, nu, tol = 0), "positive")
  expect_identical(diag_fun(-1, nu)$status, "invalid input")
})

test_that("R10c: COM pmf, cdf, quantiles and draws agree at the solved rate", {
  mu <- 4.2; nu <- 1.8
  x <- 0:60
  p <- dcom(x, mu = mu, nu = nu)
  expect_equal(sum(p), 1, tolerance = 1e-12)
  expect_equal(sum(x * p), mu, tolerance = 1e-9)
  expect_equal(pcom(x, mu = mu, nu = nu), cumsum(p), tolerance = 1e-12)
  u <- c(0.01, 0.3, 0.5, 0.77, 0.999)
  q <- qcom(u, mu = mu, nu = nu)
  cdf <- cumsum(p)
  expect_equal(q, vapply(u, function(v) min(x[cdf >= v - 1e-14]),
                         numeric(1)))
  set.seed(11)
  expect_equal(mean(rcom(20000, mu = mu, nu = nu)), mu, tolerance = 0.02)
})

# ---- R01: rootogram ---------------------------------------------------------------------

test_that("R01: rootogram expected frequencies do not depend on the cutoff", {
  d <- data.frame(y = c(0L, 10L))
  fit <- suppressWarnings(countreg(y ~ 1, data = d, family = "POISSON",
                                   method = "BFGS"))
  expect_true(isTRUE(fit$model$converged))
  full <- rootogram(fit, max_count = 10, plot = FALSE)
  cut <- rootogram(fit, max_count = 5, plot = FALSE)
  m <- match(cut$count, full$count)
  expect_equal(cut$expected, full$expected[m], tolerance = 1e-10)
  mu <- predict(fit, type = "response")
  expect_equal(cut$expected[cut$count == 0], sum(dpois(0, mu)),
               tolerance = 1e-10)
  # The observed outcome 10 is in the tail, not in a shown bin
  expect_equal(sum(cut$observed), 1)
  expect_equal(attr(cut, "tail")[["observed"]], 1)
  expect_equal(attr(cut, "tail")[["expected"]],
               sum(ppois(5, mu, lower.tail = FALSE)), tolerance = 1e-10)
  # With an overflow bin both frequencies add up to the total
  ov <- rootogram(fit, max_count = 5, overflow = TRUE, plot = FALSE)
  expect_equal(nrow(ov), 7)
  expect_equal(sum(ov$observed), 2)
  expect_equal(sum(ov$expected), 2, tolerance = 1e-10)
  expect_equal(ov$expected[1:6], cut$expected)
})

test_that("R01: rootogram with weights, new data and missing outcomes", {
  set.seed(2)
  d <- data.frame(x = rnorm(60))
  d$y <- rpois(60, exp(0.5 + 0.4 * d$x))
  d$w <- sample(1:3, 60, replace = TRUE)
  fit <- countreg(y ~ x, data = d, family = "POISSON", weights = "w",
                  method = "BFGS")
  mu <- predict(fit, type = "response")
  r <- rootogram(fit, max_count = 2, overflow = TRUE, plot = FALSE)
  expect_equal(r$expected[1:3],
               vapply(0:2, function(k) sum(d$w * dpois(k, mu)), numeric(1)),
               tolerance = 1e-8)
  expect_equal(sum(r$observed), sum(d$w))
  expect_equal(sum(r$expected), sum(d$w), tolerance = 1e-8)
  # New data: unweighted rows; a missing outcome is left out of both
  nd <- d[1:10, ]
  nd$y[3] <- NA
  rn <- rootogram(fit, newdata = nd, max_count = 1, overflow = TRUE,
                  plot = FALSE)
  mu_n <- predict(fit, newdata = nd[-3, ])
  expect_equal(sum(rn$observed), 9)
  expect_equal(rn$expected[1], sum(dpois(0, mu_n)), tolerance = 1e-8)
  expect_equal(sum(rn$expected), 9, tolerance = 1e-8)
})

test_that("R01: rows with non-finite fitted probabilities are reported", {
  set.seed(5)
  d <- data.frame(x = rnorm(30))
  d$y <- rpois(30, exp(0.2 + 0.3 * d$x))
  fit <- countreg(y ~ x, data = d, family = "POISSON", method = "BFGS")
  nd <- d
  nd$x[1] <- 1e6   # exp overflows: probabilities are not finite
  expect_warning(r <- rootogram(fit, newdata = nd, max_count = 3,
                                overflow = TRUE, plot = FALSE),
                 "non-finite fitted means")
  expect_equal(sum(r$observed), 29)
  expect_error(rootogram(fit, max_count = -1, plot = FALSE), "max_count")
})

# ---- R02: frequency-weight bootstrap ----------------------------------------------

test_that("R02: bootstrap frequencies resample the represented observations", {
  freq_fun <- getFromNamespace(".fcr_boot_frequencies", "flexCountReg")
  boot <- getFromNamespace(".fcr_bootstrap", "flexCountReg")
  y <- c(1, 9)
  w <- c(100L, 1L)
  # Coupled draws: a resample of the 101 expanded observations gives the
  # compressed frequencies by tabulation, with the same weighted mean
  set.seed(20261116)
  expanded_row <- rep(seq_along(y), times = w)
  sampled_row <- sample(expanded_row, length(expanded_row), replace = TRUE)
  new_w <- tabulate(sampled_row, nbins = length(y))
  expect_equal(weighted.mean(y, new_w), mean(y[sampled_row]))
  # The package's draws have the multinomial distribution: the variance of
  # the weighted mean is that of the expanded bootstrap (about 0.0062), not
  # the two-row bootstrap (about 11.8)
  set.seed(1)
  seen <- list()
  draws <- boot(n = 2, B = 4000, coef_names = "m", freq = w,
                refit = function(idx, wt) {
                  seen[[length(seen) + 1]] <<- list(idx = idx, wt = wt)
                  structure(c(m = sum(y[idx] * wt) / sum(wt)),
                            converged = TRUE)
                })
  v <- var(draws[, "m"])
  expect_equal(v, 0.00621178, tolerance = 0.1)
  # Each refit gets positive new frequencies that replace the weights and
  # sum to the represented sample size
  expect_true(all(vapply(seen, function(s) sum(s$wt) == 101 &&
                           all(s$wt > 0), logical(1))))
  expect_equal(sum(freq_fun(c(3, 0, 2))), 5)
  set.seed(9)
  expect_true(all(replicate(50, freq_fun(c(3, 0, 2))[2]) == 0))
})

test_that("R02: non-frequency weights are rejected for bootstrap SEs", {
  chk <- getFromNamespace(".fcr_check_boot_freq", "flexCountReg")
  expect_error(chk(c(1.5, 2), 2), "frequency weights")
  expect_error(chk(c(0, 0), 2), "frequency weights")
  expect_error(chk(c(-1, 2), 2), "frequency weights")
  expect_true(chk(c(0, 2), 2))
  d <- data.frame(x = rnorm(20), w = runif(20))
  d$y <- rpois(20, 2)
  expect_error(countreg(y ~ x, data = d, family = "POISSON", weights = "w",
                        stderr = "boot", bootstraps = 5),
               "frequency weights")
})

test_that("R02: a compressed replicate refit equals the expanded refit", {
  set.seed(7)
  d <- data.frame(f = factor(sample(c("a", "b", "c"), 40, TRUE)),
                  x = rnorm(40))
  d$y <- rpois(40, exp(0.3 + 0.5 * d$x + (d$f == "b")))
  d$w <- sample(0:4, 40, replace = TRUE)
  d$w[d$f == "a"][1] <- 2   # every level keeps positive weight
  freq_fun <- getFromNamespace(".fcr_boot_frequencies", "flexCountReg")
  set.seed(8)
  repeat {
    nw <- freq_fun(d$w)
    if (all(tapply(nw, d$f, sum) > 0)) break
  }
  pos <- which(nw > 0)
  comp <- countreg(y ~ f + x, data = d[pos, ], family = "POISSON",
                   weights = nw[pos], method = "BFGS")
  ex <- d[rep(seq_len(nrow(d)), nw), ]
  expd <- countreg(y ~ f + x, data = ex, family = "POISSON",
                   method = "BFGS")
  expect_equal(coef(comp), coef(expd), tolerance = 1e-5)
  expect_equal(as.numeric(logLik(comp)), as.numeric(logLik(expd)),
               tolerance = 1e-8)
  expect_equal(attr(logLik(comp), "nobs"), nrow(ex))
})

test_that("R02: countreg bootstrap with weights replaces the weights", {
  set.seed(12)
  d <- data.frame(x = rnorm(30))
  d$y <- rpois(30, exp(0.4 + 0.3 * d$x))
  d$w <- sample(1:5, 30, replace = TRUE)
  set.seed(4)
  m <- countreg(y ~ x, data = d, family = "POISSON", weights = "w",
                method = "BFGS", stderr = "boot", bootstraps = 30)
  expect_equal(m$model$n_boot_ok, 30)
  # The bootstrap SEs are close to those of the expanded data, not to the
  # row bootstrap of the compressed data
  ex <- d[rep(seq_len(30), d$w), ]
  h <- countreg(y ~ x, data = ex, family = "POISSON", method = "BFGS")
  ratio <- sqrt(diag(vcov(m))) / sqrt(diag(vcov(h)))
  expect_true(all(ratio > 0.6 & ratio < 1.6))
})

# ---- Panel frequency weights for every estimator --------------------------------

.panel_weight_data <- function(G = 60, Tn = 4, seed = 1) {
  set.seed(seed)
  d <- data.frame(id = rep(seq_len(G), each = Tn), x = rnorm(G * Tn))
  u <- rep(stats::rgamma(G, 2, 2), each = Tn)
  d$y <- stats::rnbinom(G * Tn, mu = exp(0.8 + 0.5 * d$x) * u, size = 1.5)
  wp <- sample(1:3, G, replace = TRUE)
  d$w <- rep(wp, each = Tn)
  ex <- do.call(rbind, lapply(seq_len(G), function(g) {
    b <- d[d$id == g, ]
    do.call(rbind, lapply(seq_len(wp[g]), function(k) {
      b$id <- g * 100 + k
      b
    }))
  }))
  list(d = d, ex = ex, n_ex = sum(wp))
}

test_that("panel models: weights equal expanded panels", {
  pw <- .panel_weight_data()
  for (fn in c("repois", "renb", "poisLind.re", "negmultinom")) {
    f <- get(fn)
    a <- suppressWarnings(f(y ~ x, group_var = "id", data = pw$d,
                            weights = "w", method = "BFGS"))
    b <- suppressWarnings(f(y ~ x, group_var = "id", data = pw$ex,
                            method = "BFGS"))
    expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)),
                 tolerance = 1e-7, label = paste(fn, "logLik"))
    expect_equal(coef(a), coef(b), tolerance = 1e-3,
                 label = paste(fn, "coefficients"))
    expect_equal(attr(logLik(a), "nobs"), pw$n_ex)
    expect_equal(BIC(a), BIC(b), tolerance = 1e-6)
  }
  # Robust covariance (analytic or well-conditioned numeric Hessians)
  for (fn in c("repois", "negmultinom")) {
    f <- get(fn)
    a <- f(y ~ x, group_var = "id", data = pw$d, weights = "w",
           method = "BFGS", stderr = "robust")
    b <- f(y ~ x, group_var = "id", data = pw$ex, method = "BFGS",
           stderr = "robust")
    expect_equal(sqrt(diag(vcov(a))), sqrt(diag(vcov(b))), tolerance = 0.02,
                 label = paste(fn, "robust SE"))
  }
})

test_that("panel models: weights must be constant within panels", {
  pw <- .panel_weight_data(G = 10)
  d <- pw$d
  d$w[1] <- 7
  for (fn in c("repois", "renb", "poisLind.re", "negmultinom")) {
    expect_error(get(fn)(y ~ x, group_var = "id", data = d, weights = "w",
                         method = "BFGS"),
                 "constant within each panel")
  }
  expect_error(countreg.rp(y ~ 1, rpar_formula = ~ x, data = d,
                           group_var = "id", weights = "w", ndraws = 20),
               "constant within each panel")
})

test_that("panel models: zero-weight panels drop out of the likelihood", {
  pw <- .panel_weight_data(G = 30, seed = 4)
  d <- pw$d
  d$w <- 1
  d$w[d$id %in% 1:5] <- 0
  a <- negmultinom(y ~ x, group_var = "id", data = d, weights = "w")
  b <- negmultinom(y ~ x, group_var = "id", data = d[!d$id %in% 1:5, ])
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)), tolerance = 1e-8)
  expect_equal(coef(a), coef(b), tolerance = 1e-5)
})

test_that("panel bootstrap draws panel frequencies and strict acceptance", {
  pw <- .panel_weight_data(G = 30, seed = 6)
  set.seed(2)
  m <- repois(y ~ x, group_var = "id", data = pw$d, weights = "w",
              method = "BFGS", stderr = "boot", bootstraps = 10)
  expect_true(m$model$n_boot_ok >= 2)
  expect_true(all(is.finite(sqrt(diag(vcov(m))))))
  set.seed(2)
  nm <- negmultinom(y ~ x, group_var = "id", data = pw$d, weights = "w",
                    stderr = "boot", bootstraps = 10)
  expect_true(nm$model$n_boot_ok >= 2)
  expect_s3_class(nm$model$boot_failures, "table")
})

test_that("renb log-likelihood stays finite and below zero for large a, b", {
  # Before the regrouping, lgamma terms of size ~1e203 cancelled and the
  # optimizer found a spurious log-likelihood of about +1e203
  ll <- getFromNamespace("renb_ll", "flexCountReg")
  y <- c(3, 2, 1, 2)
  a <- 4.4e118; b <- 2.3e200
  mu <- rep(5, 4) * (a - 1) / b
  v <- ll(y, mu, a, b, rep(1, 4))
  expect_true(is.finite(v))
  expect_lt(v, 0)
  # The value changes smoothly with a (no rounding noise of size ~1e187)
  v2 <- ll(y, mu * (1 + 1e-6), a * (1 + 1e-6), b, rep(1, 4))
  expect_lt(abs(v2 - v), 1e-3)
  # Huge panel means (M ~ 1e300, as tried by optimizers): finite and
  # negative, with no cancellation of ~1e300 terms
  v3 <- ll(y, c(1e299, 2e300, 5e298, 1e299), 2.7e5, 1.2e-3, rep(1, 4))
  expect_true(is.finite(v3) && v3 < 0)
  v4 <- ll(y, c(1e299, 2e300, 5e298, 1e299) * (1 + 1e-6), 2.7e5, 1.2e-3,
           rep(1, 4))
  expect_lt(abs(v4 - v3), 1e-3 * abs(v3))
  # Moderate values agree with the direct formula
  a <- 3; b <- 2; mu <- c(1.5, 2, 0.7, 1.1)
  direct <- lgamma(a + b) + lgamma(a + sum(mu)) + lgamma(b + sum(y)) -
    lgamma(a) - lgamma(b) - lgamma(a + b + sum(mu) + sum(y)) +
    sum(lgamma(mu + y) - lgamma(mu) - lgamma(y + 1))
  expect_equal(ll(y, mu, a, b, rep(1, 4)), direct, tolerance = 1e-12)
})

# ---- R05: GP domains ----------------------------------------------------------------

test_that("R05: GP moments are NaN outside the parameter space", {
  mom <- getFromNamespace(".gp1_moments", "flexCountReg")
  # GP2 with a = -0.5 and mu = 3: 1 + a mu = -0.5, no distribution
  m <- mom(3, -0.5 * 3)
  expect_true(is.nan(m$mean) && is.nan(m$var))
  # At and just above the bound (1 + a mu = 0, slightly positive)
  expect_true(is.nan(mom(2, -1)$mean))
  ok <- mom(2, -0.999)
  expect_true(is.finite(ok$mean) && is.finite(ok$var))
  expect_true(is.nan(mom(-1, 0.2)$mean))
  expect_true(is.na(mom(NA, 0.2)$mean))
  # Valid negative dispersion: moments of the normalized PMF by direct sums
  for (pr in list(c(3, -0.2), c(0.7, -0.45), c(12, -0.05))) {
    y <- 0:400
    p <- dgp1(y, mu = pr[1], phi = pr[2])
    mm <- mom(pr[1], pr[2])
    expect_equal(sum(p), 1, tolerance = 1e-10)
    expect_equal(mm$mean, sum(y * p), tolerance = 1e-8)
    expect_equal(mm$var, sum(y^2 * p) - sum(y * p)^2, tolerance = 1e-8)
  }
  # GP2 in the registry: rows past the bound are NaN, with a warning, and
  # the other rows keep their (aligned) values
  fm <- getFromNamespace(".fcr_family_mean", "flexCountReg")
  expect_warning(v <- fm("GP2", c(1, 3, 1.5), alpha = 0.5),   # a = -0.5
                 "outside the parameter space")
  expect_true(is.finite(v[1]) && is.nan(v[2]) && is.finite(v[3]))
})

test_that("R05: GP2 random-parameter models require a nonnegative dispersion", {
  set.seed(4)
  d <- data.frame(x = rnorm(150))
  d$y <- rpois(150, exp(0.5 + 0.3 * d$x))
  m <- suppressWarnings(countreg.rp(y ~ 1, rpar_formula = ~ x, data = d,
                                    family = "GP2", ndraws = 30,
                                    method = "NM", max.iters = 400))
  # alpha = a + 1 >= 1: the estimate stays in the allowed region
  expect_gte(exp(coef(m)[["ln(alpha+1)"]]), 1 - 1e-8)
})

# ---- R09: invalid probabilities ---------------------------------------------------

test_that("R09a: zero-truncated draws fail on invalid probabilities", {
  rt <- getFromNamespace(".fcr_rtrunc_positive", "flexCountReg")
  obj <- list(model = list(family = "TEST"))
  parts <- list(mu_link = c(1, 2, 3, 4))
  # Rows 1-4: valid Poisson(2); a NaN pmf at k = 2; P(Y = 0) = 1 (no mass
  # above zero); an Inf pmf at k = 1
  dens <- function(object, parts, k, log = FALSE) {
    id <- parts$mu_link
    v <- dpois(k, 2, log = TRUE) + 0 * id
    if (k == 0) v[id == 3] <- 0
    if (k == 2) v[id == 2] <- NaN
    if (k == 1) v[id == 4] <- Inf
    v
  }
  set.seed(1)
  expect_warning(out <- rt(obj, parts, 1:4, dens = dens),
                 "TEST draws failed for 3 row")
  expect_true(out[1] >= 1 && out[1] == round(out[1]))
  expect_true(all(is.na(out[2:4])) || (out[2] == 1 && all(is.na(out[3:4]))))
  # The count limit is a failure, not an invented count
  dens0 <- function(object, parts, k, log = FALSE)
    rep(if (k == 0) log(0.5) else -Inf, length(parts$mu_link))
  expect_warning(out <- rt(obj, list(mu_link = 1), 1, max_count = 5,
                           dens = dens0), "no draw within 5")
  expect_true(is.na(out))
})

test_that("R09a: positive Poisson draws follow the truncated distribution", {
  rt <- getFromNamespace(".fcr_rtrunc_positive", "flexCountReg")
  obj <- list(model = list(family = "POISSON"))
  lam <- 0.3
  dens <- function(object, parts, k, log = FALSE)
    dpois(k, lam, log = TRUE) + 0 * parts$mu_link
  set.seed(5)
  n <- 20000
  out <- rt(obj, list(mu_link = rep(1, n)), seq_len(n), dens = dens)
  expect_true(all(out >= 1 & out == round(out)))
  ref <- dpois(1:4, lam) / (1 - dpois(0, lam))
  expect_equal(as.vector(table(factor(out, 1:4))) / n, ref, tolerance = 0.03)
  # Tiny positive-count mass: still positive draws
  lam <- 1e-10
  out <- rt(obj, list(mu_link = rep(1, 50)), 1:50, dens = dens)
  expect_true(all(out == 1))
})

test_that("R09b: CDF endpoints are not valid outputs for invalid parameters", {
  q <- c(-1, Inf, 2)
  expect_equal(ppois(q, 2), c(0, 1, ppois(2, 2)))
  for (f in list(
    function(q, v) pgp1(q, mu = v, phi = 0.2),
    function(q, v) pgp2(q, mu = v, alpha = 0.1),
    function(q, v) pcom(q, mu = v, nu = 1.5),
    function(q, v) pgwar(q, mu = v, k = 2, rho = 3),
    function(q, v) ppLnorm(q, mean = v, sigma = 0.5),
    function(q, v) ppge(q, mean = v, shape = 1.5),
    function(q, v) ppinvgamma(q, mu = v, eta = 0.5),
    function(q, v) ppinvgaus(q, mu = v, eta = 0.5),
    function(q, v) pplind(q, mean = v, theta = 1.5),
    function(q, v) pplindLnorm(q, mean = v, theta = 1.5, sigma = 0.5),
    function(q, v) psichel(q, mu = v, sigma = 1, gamma = -0.5),
    function(q, v) plindley(q, theta = v))) {
    valid <- f(q, 3)
    expect_equal(valid[1:2], c(0, 1))
    expect_true(valid[3] > 0 && valid[3] < 1)
    bad <- suppressWarnings(f(q, -1))
    expect_true(all(is.nan(bad) | is.na(bad)))
    expect_true(all(is.na(f(q, NA))))
    # Aligned: one invalid element does not affect the others
    mixed <- suppressWarnings(f(c(2, 2), c(3, -1)))
    expect_equal(mixed[1], valid[3])
    expect_true(is.na(mixed[2]))
  }
})

test_that("R09b: Lindley tails keep their precision", {
  th <- 2; q <- 30
  # P(X > q) = (1 + theta q / (1 + theta)) exp(-theta q)
  ref <- log1p(th * q / (1 + th)) - th * q
  expect_equal(plindley(q, th, lower.tail = FALSE, log.p = TRUE), ref)
  expect_equal(plindley(1e-12, th, log.p = TRUE),
               log(th^2 / (1 + th) * 1e-12), tolerance = 1e-6)
  expect_equal(plindley(c(-1, 0, Inf), th), c(0, 0, 1))
  expect_warning(expect_true(is.nan(plindley(1, -2))), "NaN")
})

# ---- R04: mixture log probabilities -----------------------------------------------

test_that("R04: PGE with shape 1 is geometric for every count and scale", {
  y <- c(0, 1, 7, 100, 1000, 20000)
  for (mu in c(0.05, 1, 30, 2500)) {
    for (sc in c(0.5, 1, 2)) {
      expect_equal(dpge(y, mean = mu, shape = 1, scale = sc, log = TRUE),
                   dgeom(y, prob = 1 / (1 + mu), log = TRUE),
                   tolerance = 1e-9)
    }
  }
  # The review's example: y = 1000 with mean 1. Halton draws give about
  # -3847; the exact value is about -693.84
  expect_equal(dpge(1000, mean = 1, shape = 1, log = TRUE),
               dgeom(1000, 0.5, log = TRUE), tolerance = 1e-10)
  old <- dpge(1000, mean = 1, shape = 1, log = TRUE, engine = "halton")
  expect_lt(old, -3000)
  # Tails and quantiles from the same PMF
  expect_equal(ppge(c(5, 2000), mean = 4, shape = 1, scale = 2,
                    lower.tail = FALSE, log.p = TRUE),
               pgeom(c(5, 2000), 1 / 5, lower.tail = FALSE, log.p = TRUE),
               tolerance = 1e-9)
  expect_equal(qpge(c(0.1, 0.5, 0.99), mean = 4, shape = 1, scale = 2),
               qgeom(c(0.1, 0.5, 0.99), 1 / 5))
})

test_that("R04: PW with alpha 1 is geometric; PGE P(0) is a beta function", {
  y <- c(0, 2, 40, 3000)
  # Weibull(1, scale) mixing is exponential: geometric with the same mean
  expect_equal(dpoisweibull(y, lambda = 2, alpha = 1, sigma = 1.5,
                            log = TRUE),
               dgeom(y, 1 / (1 + 3), log = TRUE), tolerance = 1e-9)
  # P(Y = 0) = E[exp(-c X)] = a B(a, c + 1) for X ~ GE(a, 1)
  for (a in c(0.3, 2.5)) {
    for (mu in c(0.1, 40)) {
      cc <- mu / (digamma(a + 1) - digamma(1))
      expect_equal(dpge(0, mean = mu, shape = a, log = TRUE),
                   log(a) + lbeta(a, cc + 1), tolerance = 1e-10)
    }
  }
})

test_that("R04: PLN and PLL match independent references in the tails", {
  # Reference: the integral over the normal mixing variable, on a fine grid
  ref <- function(y, lam, sig, kern) {
    z <- seq(-45, 45, length.out = 300001)
    lp <- kern(y, lam * exp(sig * z)) + dnorm(z, log = TRUE)
    m <- max(lp)
    m + log(sum(exp(lp - m)) * (z[2] - z[1]))
  }
  pk <- function(y, l) dpois(y, l, log = TRUE)
  for (y in c(0, 3, 60, 800)) {
    for (pr in list(c(2, 0.4), c(0.05, 1.5), c(300, 0.8))) {
      expect_equal(dpLnorm(y, pr[1], pr[2], log = TRUE),
                   ref(y, pr[1], pr[2], pk), tolerance = 1e-8)
    }
  }
  th <- 1.7
  plk <- function(y, l) {
    lam <- l * th * (th + 1) / (th + 2)
    2 * log(th) - log1p(th) + y * log(lam) + log(th + lam + y + 1) -
      (y + 2) * log(th + lam)
  }
  for (y in c(0, 5, 400)) {
    expect_equal(dplindLnorm(y, 2.5, th, 0.6, log = TRUE),
                 ref(y, 2.5, 0.6, plk), tolerance = 1e-8)
  }
  # sigma = 0: Poisson and Poisson-Lindley
  expect_equal(dpLnorm(0:5, 2, 0), dpois(0:5, 2), tolerance = 1e-12)
  expect_equal(dplindLnorm(0:5, 2, th, 0), dplind(0:5, mean = 2, theta = th),
               tolerance = 1e-10)
})

test_that("R04: mixtures have unit mass, the right mean and a consistent CDF", {
  y <- 0:3000
  cases <- list(
    PLN = function(...) dpLnorm(y, 3, 0.7, ...),
    PGE = function(...) dpge(y, 3, 0.6, ...),
    PW = function(...) dpoisweibull(y, lambda = 3 / gamma(1 + 1 / 0.8),
                                    alpha = 0.8, sigma = 1, ...),
    PLL = function(...) dplindLnorm(y, 3, 1.2, 0.5, ...))
  means <- c(PLN = 3 * exp(0.7^2 / 2), PGE = 3, PW = 3,
             PLL = 3 * exp(0.5^2 / 2))
  for (nm in names(cases)) {
    p <- cases[[nm]]()
    expect_equal(sum(p), 1, tolerance = 1e-9, label = nm)
    expect_equal(sum(y * p), means[[nm]], tolerance = 1e-7, label = nm)
  }
  expect_equal(ppLnorm(0:20, 3, 0.7), cumsum(dpLnorm(0:20, 3, 0.7)),
               tolerance = 1e-12)
  expect_equal(pplindLnorm(0:20, 3, 1.2, sigma = 0.5),
               cumsum(dplindLnorm(0:20, 3, 1.2, 0.5)), tolerance = 1e-12)
  expect_equal(ppoisweibull(0:20, lambda = 2, alpha = 0.8, sigma = 1),
               cumsum(dpoisweibull(0:20, lambda = 2, alpha = 0.8, sigma = 1)),
               tolerance = 1e-12)
  # Simulated moments agree
  set.seed(3)
  expect_equal(mean(rpge(40000, 3, 0.6)), 3, tolerance = 0.03)
  expect_equal(mean(rplindLnorm(40000, 3, 1.2, 0.5)), means[["PLL"]],
               tolerance = 0.03)
})

test_that("R04: no underflow to -Inf or probability floor for extreme counts", {
  lp <- dpLnorm(c(1e4, 1e5), mean = 1, sigma = 0.3, log = TRUE)
  expect_true(all(is.finite(lp)))
  # Reference: the integral on a grid around the mode
  for (i in 1:2) {
    y <- c(1e4, 1e5)[i]
    z <- seq(20, 50, length.out = 300001)
    h <- dpois(y, exp(0.3 * z), log = TRUE) + dnorm(z, log = TRUE)
    m <- max(h)
    expect_equal(lp[i], m + log(sum(exp(h - m)) * (z[2] - z[1])),
                 tolerance = 1e-8)
  }
  expect_true(is.finite(dpge(5e4, mean = 0.5, shape = 0.4, log = TRUE)))
  expect_true(is.finite(dplindLnorm(1e4, 1, 2, 0.3, log = TRUE)))
  # Missing and invalid values stay aligned
  expect_warning(v <- dpLnorm(c(1, 1, 1), mean = c(2, -1, NA), sigma = 0.5),
                 "NaN")
  expect_true(is.finite(v[1]) && is.nan(v[2]) && is.na(v[3]))
  expect_equal(dpLnorm(c(-1, 1.5), 2, 0.5), c(0, 0))
})

test_that("R04: fixed and random-parameter likelihoods use the quadrature", {
  set.seed(8)
  d <- data.frame(x = rnorm(120))
  d$y <- rpge(120, mean = exp(0.4 + 0.3 * d$x), shape = 2)
  m <- countreg(y ~ x, data = d, family = "PGE", method = "BFGS")
  mu <- exp(drop(cbind(1, d$x) %*% coef(m)[1:2]))
  ll <- sum(dpge(d$y, mean = mu, shape = exp(coef(m)[["ln(shape)"]]),
                 log = TRUE))
  expect_equal(as.numeric(logLik(m)), ll, tolerance = 1e-8)
  expect_identical(m$model$engine, "quadrature")
  # With a zero random-parameter SD, the random-parameter likelihood is the
  # fixed one
  dens <- getFromNamespace(".fcr_family_registry", "flexCountReg")$PGE$dens
  expect_equal(dens(d$y, mu, exp(coef(m)[["ln(shape)"]]), NULL, NULL, NULL,
                    "quadrature", log = TRUE),
               dpge(d$y, mean = mu, shape = exp(coef(m)[["ln(shape)"]]),
                    log = TRUE))
})

# ---- R07: moment existence -------------------------------------------------------

test_that("R07: existence of E[exp(x beta)] is decided without the MGF", {
  ex <- getFromNamespace(".fcr_rp_mean_exists", "flexCountReg")
  expect_equal(ex("ln", c(-2, 0, 1), 0.3, 0.5), c(TRUE, TRUE, FALSE))
  expect_equal(ex("ln", 1, 0.3, 0), TRUE)           # fixed coefficient
  expect_equal(ex("g", c(1, 5), 1, 0.5), c(TRUE, FALSE))  # t s^2 / mu < 1
  expect_equal(ex("g", 5, 1, 0), TRUE)
  expect_equal(ex("g", 1, -1, 0.5), FALSE)
  expect_equal(ex("n", 1e6, 3, 2), TRUE)            # finite, may overflow
  expect_true(is.na(ex("n", NA, 3, 2)))
  lm <- getFromNamespace("rp_log_mgf", "flexCountReg")
  # A normal coefficient with a huge finite MGF: finite on the log scale
  expect_equal(lm("n", 400, 2, 1), 400 * 2 + 0.5 * 400^2)
  expect_equal(lm("g", 0.5, 1, 0.5), -4 * log(1 - 0.5 * 0.25 / 1))
})

test_that("R07: Simulated lognormal predictions do not need lamW", {
  set.seed(9)
  d <- data.frame(x = -abs(rnorm(150)), id = rep(1:50, each = 3))
  d$y <- rpois(150, exp(0.5 - 0.4 * d$x))
  m <- suppressWarnings(countreg.rp(y ~ 1, rpar_formula = ~ -1 + x,
                                    data = d, family = "POISSON",
                                    rpardists = c(x = "ln"), ndraws = 50,
                                    method = "BFGS"))
  # Make the MGF unavailable: Simulated must not call it
  testthat::local_mocked_bindings(
    rp_log_mgf = function(...) stop("MGF evaluated"),
    .package = "flexCountReg")
  p <- predict(m, newdata = d[1:5, ], method = "Simulated")
  expect_true(all(is.finite(p)))
  expect_error(predict(m, newdata = d[1:5, ], method = "Exact"),
               "MGF evaluated")
  # A positive x gives an infinite expectation, with a warning
  nd <- d[1:2, ]
  nd$x <- c(-0.5, 0.5)
  expect_warning(p <- predict(m, newdata = nd, method = "Simulated"),
                 "expected value does not exist")
  expect_true(is.finite(p[1]) && is.infinite(p[2]))
})

test_that("R07: overflow of a finite expectation is not called nonexistence", {
  set.seed(10)
  d <- data.frame(x = rnorm(120))
  d$y <- rpois(120, exp(0.3 + 0.2 * d$x))
  m <- suppressWarnings(countreg.rp(y ~ 1, rpar_formula = ~ -1 + x,
                                    data = d, family = "POISSON",
                                    ndraws = 30, method = "BFGS"))
  nd <- data.frame(x = c(1, 1e4))
  expect_warning(p <- predict(m, newdata = nd, method = "Exact"), "overflow")
  expect_true(is.finite(p[1]) && is.infinite(p[2]))
})

# ---- R06: marginal-effect bootstrap refits ------------------------------------

test_that("R06: a dot-formula refit keeps the original design", {
  set.seed(21)
  d <- data.frame(y = rpois(80, 2), x1 = rnorm(80), x2 = rnorm(80))
  m <- countreg(y ~ ., data = d, family = "POISSON", method = "BFGS")
  info <- .margEff_boot_units(m, .fcr_training_rows(m), NULL, TRUE)
  rs <- .margEff_resample(info)
  rs$data$.fcr_boot_id <- 1   # a helper column must not become a predictor
  r <- .margEff_refit(m, rs$data)
  expect_identical(names(coef(r)), names(coef(m)))
})

test_that("R06: superclusters drawn twice keep their panels distinct", {
  # Supercluster A holds panels 1 and 2 (3 rows each); B holds panel 3
  d <- data.frame(sc = rep(c("A", "A", "B"), each = 3),
                  pid = rep(1:3, each = 3), t = rep(1:3, 3))
  info <- list(rows = d, panel_index = d$pid,
               unit_rows = split(seq_len(9), factor(d$sc)),
               unit_w = NULL, panel_col = ".fcr_boot_id")
  # Force A to be drawn twice
  testthat::local_mocked_bindings(
    sample.int = function(n, size, replace) c(1L, 1L), .package = "base")
  rs <- .margEff_resample(info)
  p <- rs$data$.fcr_boot_id
  expect_equal(length(unique(p)), 4)
  # Each new panel holds the rows of one original panel, in order
  by_p <- split(rs$data, p)
  expect_true(all(vapply(by_p, function(b) length(unique(b$pid)) == 1 &&
                           identical(b$t, 1:3), NA)))
})

test_that("R06: units, weights and panels follow the analysis rows", {
  set.seed(22)
  d <- data.frame(id = rep(1:30, each = 3), x = rnorm(90))
  d$y <- rpois(90, exp(0.2 + 0.4 * d$x))
  d$x[5] <- NA   # dropped row
  w <- rep(sample(1:3, 30, TRUE), each = 3)
  m <- suppressMessages(repois(y ~ x, group_var = "id", data = d,
                               weights = w, method = "BFGS"))
  info <- .margEff_boot_units(m, .fcr_training_rows(m), NULL, TRUE)
  expect_equal(nrow(info$rows), 89)
  expect_length(info$unit_rows, 30)
  expect_equal(info$unit_w, w[!duplicated(d$id)])
  set.seed(3)
  rs <- .margEff_resample(info)
  # The represented panels are preserved: sum of panel frequencies
  pw <- tapply(rs$weights, rs$data[[info$panel_col]], `[`, 1)
  expect_equal(sum(pw), sum(w[!duplicated(d$id)]))
  # cluster_var with weights is not supported
  expect_error(.margEff_boot_units(m, info$rows, "id", TRUE),
               "cannot be combined with frequency weights")
  # A panel split across clusters is rejected
  m2 <- suppressMessages(repois(y ~ x, group_var = "id", data = d,
                                method = "BFGS"))
  rows2 <- .fcr_training_rows(m2)
  rows2$cl <- seq_len(nrow(rows2)) %% 2
  expect_error(.margEff_boot_units(m2, rows2, "cl", TRUE),
               "within one `cluster_var`")
})

test_that("R06: refits use stored values, not changed external objects", {
  set.seed(23)
  d <- data.frame(x = rnorm(100))
  d$y <- rnbinom(100, mu = exp(0.5 + 0.3 * d$x), size = 2)
  fam <- "NB2"
  m <- countreg(y ~ x, data = d, family = fam, method = "BFGS")
  fam <- "POISSON"
  r <- .margEff_refit(m, d)
  expect_identical(r$model$family, "NB2")
  expect_equal(coef(r), coef(m), tolerance = 1e-6)
})

test_that("R06: the panel column name does not overwrite user columns", {
  set.seed(24)
  d <- data.frame(id = rep(1:20, each = 3), x = rnorm(60),
                  .fcr_boot_id = 99)
  d$y <- rpois(60, exp(0.3 + 0.3 * d$x))
  m <- repois(y ~ x, group_var = "id", data = d, method = "BFGS")
  info <- .margEff_boot_units(m, .fcr_training_rows(m), NULL, TRUE)
  expect_identical(info$panel_col, ".fcr_boot_id_1")
  rs <- .margEff_resample(info)
  expect_true(all(rs$data$.fcr_boot_id == 99))
})

test_that("R06: a converged panel bootstrap and weighted effects", {
  set.seed(25)
  d <- data.frame(id = rep(1:40, each = 3), x = rnorm(120))
  d$y <- rpois(120, exp(0.3 + 0.5 * d$x) * rep(rgamma(40, 3, 3), each = 3))
  w <- rep(sample(1:3, 40, TRUE), each = 3)
  m <- repois(y ~ x, group_var = "id", data = d, weights = w,
              method = "BFGS")
  set.seed(1)
  me <- margEffTable(m, vars = "x", se = "bootstrap", bootstraps = 8,
                     digits = 10)
  expect_equal(me$successful_bootstraps, 8)
  expect_true(is.finite(me$std_error) && me$std_error > 0)
  # The effect is the frequency-weighted average over the analysis rows,
  # i.e. the average over the expanded data
  ex <- d[rep(seq_len(120), w), ]
  mu <- function(dd, h) exp(drop(cbind(1, dd$x + h) %*% coef(m)[1:2]))
  hh <- 1e-5
  ame_ex <- mean((mu(ex, hh) - mu(ex, -hh)) / (2 * hh))
  expect_equal(me$estimate, ame_ex, tolerance = 1e-5)
})

# ---- R03: LR eligibility ------------------------------------------------------------

.lr_data <- function(n = 200, seed = 31) {
  set.seed(seed)
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n),
                  id = rep(seq_len(n / 4), each = 4))
  d$y <- rpois(n, exp(0.4 + 0.3 * d$x1 + 0.2 * d$x2))
  d
}

test_that("R03: a nested pair on the same rows gives the textbook LR", {
  d <- .lr_data()
  r <- countreg(y ~ x1, d, family = "POISSON", method = "BFGS")
  f <- countreg(y ~ x1 + x2, d, family = "POISSON", method = "BFGS")
  t <- regCompTest(r, model2 = f)
  lr <- 2 * (as.numeric(logLik(f)) - as.numeric(logLik(r)))
  expect_equal(t$LR, lr, tolerance = 1e-6)
  expect_equal(t$LRdof, 1)
  expect_equal(t$LR_pvalue, pchisq(lr, 1, lower.tail = FALSE),
               tolerance = 1e-6)
  # An algebraically equivalent basis of the restricted design is nested
  r2 <- countreg(y ~ I(2 * x1 + 1), d, family = "POISSON", method = "BFGS")
  t2 <- regCompTest(r2, model2 = f)
  expect_equal(t2$LR, lr, tolerance = 1e-5)
  # A design outside the column space is not
  r3 <- countreg(y ~ x3, d, family = "POISSON", method = "BFGS")
  expect_warning(t3 <- regCompTest(r3, model2 = f), "column space")
  expect_true(is.na(t3$LR_pvalue))
  expect_match(t3$LR_reason, "column space")
})

test_that("R03: the same outcome with changed predictor data is rejected", {
  d <- .lr_data()
  d2 <- d
  d2$x1 <- d2$x1 + 0.5
  r <- countreg(y ~ x1, d, family = "POISSON", method = "BFGS")
  f <- countreg(y ~ x1, d2, family = "NB2", method = "BFGS")
  expect_error(regCompTest(r, model2 = f), "different values of `x1`")
  # Variables only one model uses are not compared
  d3 <- d
  d3$x3 <- 0
  f3 <- countreg(y ~ x1 + x2, d3, family = "POISSON", method = "BFGS")
  expect_silent(regCompTest(r, model2 = f3))
})

test_that("R03: changed panels and random-parameter distributions", {
  d <- .lr_data()
  d$id2 <- rev(d$id)
  d$id3 <- rep(seq_len(25), each = 8)
  a <- repois(y ~ x1, group_var = "id", data = d, method = "BFGS")
  b <- repois(y ~ x1 + x2, group_var = "id3", data = d, method = "BFGS")
  expect_error(regCompTest(a, model2 = b), "different panels")
  # The same partition with other labels is the same set of panels
  b2 <- repois(y ~ x1 + x2, group_var = "id2", data = d, method = "BFGS")
  expect_equal(regCompTest(a, model2 = b2)$LRdof, 1)
  pn <- suppressWarnings(countreg.rp(y ~ x2, rpar_formula = ~ -1 + x1,
                                     data = d, family = "POISSON",
                                     rpardists = c(x1 = "u"), ndraws = 30,
                                     method = "BFGS"))
  nb <- suppressWarnings(countreg.rp(y ~ x2, rpar_formula = ~ -1 + x1,
                                     data = d, family = "NB2",
                                     rpardists = c(x1 = "n"), ndraws = 30,
                                     method = "BFGS"))
  expect_warning(t <- regCompTest(pn, model2 = nb), "distributions")
  expect_true(is.na(t$LR))
})

test_that("R03: unknown or failed convergence and negative LR are refused", {
  d <- .lr_data()
  r <- countreg(y ~ x1, d, family = "POISSON", method = "BFGS")
  f <- countreg(y ~ x1 + x2, d, family = "POISSON", method = "BFGS")
  for (cv in list(NA, FALSE, NULL)) {
    f2 <- f
    f2$model["converged"] <- list(cv)
    expect_warning(t <- regCompTest(r, model2 = f2), "convergence")
    expect_true(is.na(t$LR_pvalue))
  }
  # A materially worse full model: unavailable, not a p-value of 1
  f3 <- f
  f3$model$maximum <- as.numeric(logLik(r)) - 2
  expect_warning(t <- regCompTest(r, model2 = f3), "lower log-likelihood")
  expect_true(is.na(t$LR) && is.na(t$LR_pvalue))
  lp <- getFromNamespace(".fcr_lr_pvalue", "flexCountReg")
  nest <- list(nested = TRUE, type = "interior", df = 1)
  expect_true(is.na(lp(-100, -99, nest)$LR))
  expect_equal(lp(-100, -100 + 1e-9, nest)$LR, 0)   # rounding only
  # Summary: no test against Poisson when convergence is not established
  nbm <- countreg(y ~ x1, d, family = "NB2", method = "BFGS")
  expect_false(is.null(summary(nbm)$lr_test))
  nbm$model$converged <- NA
  expect_null(summary(nbm)$lr_test)
})

test_that("R03: the generated base model must be nested", {
  d <- .lr_data()
  m <- countreg(y ~ -1 + x1, d, family = "NB2", method = "BFGS")
  t <- suppressMessages(regCompTest(m, basemodel = "Poisson"))
  expect_true(is.na(t$LR))
  expect_match(t$LR_reason, "intercept")
  t2 <- regCompTest(m, basemodel = "Poisson", variables = TRUE)
  expect_true(is.finite(t2$LR))
})

test_that("R03: older objects without sample evidence fail conservatively", {
  d <- .lr_data()
  r <- countreg(y ~ x1, d, family = "POISSON", method = "BFGS")
  f <- countreg(y ~ x1 + x2, d, family = "POISSON", method = "BFGS")
  f$model$lr_sample <- NULL
  expect_warning(t <- regCompTest(r, model2 = f), "older version")
  expect_true(is.na(t$LR))
})

# ---- R08: quantile residuals without a probability grid ------------------------

# The earlier algorithm: a dense probability grid up to each row's count
.old_qres <- function(object, u) {
  pm <- getFromNamespace(".fcr_prob_matrix", "flexCountReg")
  d <- .fcr_training_rows(object)
  y <- as.numeric(object$model$y)
  P <- pm(object, d, max(y))
  cdf <- t(apply(P, 1, cumsum))
  r <- seq_along(y)
  up <- cdf[cbind(r, y + 1)]
  lo <- ifelse(y > 0, cdf[cbind(r, pmax(y, 1))], 0)
  qnorm(lo + u * (pmin(up, 1) - lo))
}

test_that("R08: direct CDF residuals equal the grid residuals", {
  set.seed(41)
  d <- data.frame(x = rnorm(150), z = rnorm(150))
  d$y <- rnbinom(150, mu = exp(0.6 + 0.4 * d$x), size = 2)
  d$y[sample(150, 30)] <- 0
  qr <- getFromNamespace(".fcr_quantile_residuals", "flexCountReg")
  fits <- list(
    pois = countreg(y ~ x, d, family = "POISSON", method = "BFGS"),
    nb2 = countreg(y ~ x, d, family = "NB2", method = "BFGS"),
    nb1_und = suppressWarnings(countreg(y ~ x, d, family = "NB1",
                                        underreport_formula = ~ z,
                                        method = "BFGS")),
    pl = countreg(y ~ x, d, family = "PL", method = "BFGS"),   # streaming
    zinb = countreg(y ~ x, d, family = "NB2", zero_formula = ~ z,
                    zero_type = "inflated", method = "BFGS"),
    hnb = countreg(y ~ x, d, family = "NB2", zero_formula = ~ z,
                   zero_type = "hurdle", method = "BFGS"))
  u <- runif(150)
  for (nm in names(fits)) {
    expect_equal(qr(fits[[nm]], u = u), .old_qres(fits[[nm]], u),
                 tolerance = 1e-7, label = nm)
  }
  # The public method returns standard normal residuals for a good model
  set.seed(1)
  r <- residuals(fits$nb2, type = "quantile")
  expect_length(r, 150)
  expect_true(all(is.finite(r)))
})

test_that("R08: hurdle CDF is accurate when P(Y = 0) is close to 1", {
  fc <- getFromNamespace(".fcr_fitted_log_cdf", "flexCountReg")
  obj <- list(model = list(family = "POISSON"))
  lam <- 1e-9
  parts <- list(family = "POISSON", mu_link = lam, p_report = 1,
                zero_prob = 0.4, zero_type = "hurdle",
                lf0 = dpois(0, lam, log = TRUE))
  # F(0) = 1 - z, and F(1) = 1 - z P(Y > 1) / P(Y > 0)
  expect_equal(exp(fc(obj, parts, 0)), 0.6, tolerance = 1e-12)
  s1 <- ppois(1, lam, lower.tail = FALSE) / ppois(0, lam, lower.tail = FALSE)
  expect_equal(1 - exp(fc(obj, parts, 1)), 0.4 * s1, tolerance = 1e-6)
  expect_equal(fc(obj, parts, -1), -Inf)
})

test_that("R08: million-scale counts need no probability grid", {
  set.seed(42)
  n <- 300
  d <- data.frame(x = rnorm(n))
  d$y <- rpois(n, exp(13.8 + 0.01 * d$x))   # about 1e6
  m <- countreg(y ~ x, d, family = "POISSON", method = "NR")
  t <- system.time(r <- residuals(m, type = "quantile"))[["elapsed"]]
  expect_true(all(is.finite(r)))
  expect_lt(t, 10)
  expect_lt(abs(mean(r)), 0.3)
  # Overdispersed counts for NB2 (on Poisson data NB2 sits at its alpha = 0
  # boundary, where the Hessian is singular)
  d$y2 <- rnbinom(n, mu = exp(13.8 + 0.01 * d$x), size = 50)
  mn <- countreg(y2 ~ x, d, family = "NB2", method = "BFGS")
  expect_true(all(is.finite(residuals(mn, type = "quantile"))))
})

test_that("R08: families without a closed-form CDF stop before the work", {
  pc <- getFromNamespace(".fcr_parts_cdf", "flexCountReg")
  obj <- list(model = list(family = "PL", engine = "quadrature"))
  parts <- list(family = "PL", mu_link = rep(1e6, 1000), p_report = 1,
                alpha = rep(1, 1000))
  t <- system.time(
    expect_error(pc(obj, parts, rep(1e6, 1000)), "budget"))[["elapsed"]]
  expect_lt(t, 1)
})

# ---- R12: PGE parameterization and object schema --------------------------------

test_that("R12: the PGE rate cancels and the moments match the formulas", {
  a <- 1.7; mu <- 2.5
  y <- 0:2000
  p1 <- dpge(y, mean = mu, shape = a, scale = 1)
  p2 <- dpge(y, mean = mu, shape = a, scale = 2)
  expect_equal(p1, p2, tolerance = 1e-12)
  A <- digamma(a + 1) - digamma(1)
  B <- trigamma(1) - trigamma(a + 1)
  expect_equal(sum(y * p2), mu, tolerance = 1e-8)
  expect_equal(sum(y^2 * p2) - mu^2, mu + B / A^2 * mu^2, tolerance = 1e-7)
  # Simulation with rate 2 has the same mean
  set.seed(6)
  expect_equal(mean(rpge(40000, mean = mu, shape = a, scale = 2)), mu,
               tolerance = 0.03)
  # The registry variance is the formula
  v <- getFromNamespace(".fcr_family_registry", "flexCountReg")$PGE$variance
  expect_equal(v(mu, a, NULL), mu + B / A^2 * mu^2)
})

test_that("R12: objects carry a schema; older objects degrade safely", {
  expect_equal(getFromNamespace(".fcr_schema_version", "flexCountReg"), 2L)
  set.seed(51)
  d <- data.frame(x = rnorm(80))
  d$y <- rpge(80, mean = exp(0.3 + 0.3 * d$x), shape = 2)
  m <- countreg(y ~ x, d, family = "PGE", method = "BFGS")
  expect_equal(m$model$fcr_schema, 2L)
  # An older (schema 1) object: no schema, sample evidence, panels or engine
  old <- m
  old$model[c("fcr_schema", "lr_sample", "panel_index", "engine")] <- NULL
  expect_equal(predict(old, newdata = d[1:5, ]), predict(m, newdata = d[1:5, ]))
  # Probabilities use the engine it was fit with (Halton draws)
  pr <- predict(old, newdata = d[1:3, ], type = "prob", max_count = 2)
  mu <- predict(m, newdata = d[1:3, ])
  sh <- exp(coef(m)[["ln(shape)"]])
  h <- randtoolbox::halton(m$model$ndraws %||% 1500)
  expect_equal(unname(pr[, 1]),
               dpge(0, mean = mu, shape = sh, haltons = h,
                    engine = "halton"), tolerance = 1e-10)
  # LR tests need the sample evidence
  p <- countreg(y ~ x, d, family = "POISSON", method = "BFGS")
  expect_warning(t <- regCompTest(p, model2 = old), "older version")
  expect_true(is.na(t$LR))
  # Panel bootstraps need the panels
  dp <- data.frame(id = rep(1:20, each = 4), x = rnorm(80))
  dp$y <- rpois(80, exp(0.2 + 0.3 * dp$x))
  rp <- repois(y ~ x, group_var = "id", data = dp, method = "BFGS")
  rp$model$panel_index <- NULL
  expect_error(margEffTable(rp, vars = "x", se = "bootstrap", bootstraps = 2),
               "refit")
})

test_that("renb: analytic scores give an accurate, negative definite Hessian", {
  d <- washington_roads
  f <- renb(Total_crashes ~ lnaadt + speed50, d, group_var = "ID",
            method = "BFGS")
  p0 <- coef(f)
  X <- cbind(1, d$lnaadt, d$speed50)
  g <- factor(d$ID)
  ll_fun <- getFromNamespace("renb_ll", "flexCountReg")
  ll <- function(p) {
    a <- 1 + exp(p[4]); b <- exp(p[5])
    mu <- exp(drop(X %*% p[1:3])) * (a - 1) / b
    sum(ll_fun(d$Total_crashes, mu, a, b, g))
  }
  # Independent central-difference Hessian (the log-likelihood is smooth at
  # this step)
  k <- 5; h <- 1e-3; H <- matrix(0, k, k)
  for (i in 1:k) for (j in 1:k) {
    ei <- diag(k)[i, ] * h; ej <- diag(k)[j, ] * h
    H[i, j] <- (ll(p0 + ei + ej) - ll(p0 + ei - ej) - ll(p0 - ei + ej) +
                  ll(p0 - ei - ej)) / (4 * h^2)
  }
  expect_equal(unname(f$model$hessian), H, tolerance = 1e-3)
  expect_true(all(eigen(f$model$hessian)$values < 0))
  expect_true(all(diag(vcov(f)) > 0))
  # The scores vanish at the optimum
  expect_lt(max(abs(colSums(f$model$gradientObs))), 0.05)
})
