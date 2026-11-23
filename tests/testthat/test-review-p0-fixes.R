# Regression tests for the P0 defects of the November 2026 review: each can
# change results without an error.

test_that("PIG PMFs stay normalized near the Poisson limit", {
  for (eta in c(1e-8, 1e-12, 1e-16, 1e-300)) {
    for (form in c("Type 1", "Type 2")) {
      p <- dpinvgaus(0:60, mu = 3, eta = eta, form = form)
      expect_equal(sum(p), 1, tolerance = 1e-10)
      expect_equal(p[1], dpois(0, 3), tolerance = 1e-6)
    }
  }
  expect_equal(ppinvgaus(5, mu = 3, eta = 1e-16), ppois(5, 3),
               tolerance = 1e-12)
})

test_that("PIG1 fits to Poisson data never report a positive log-likelihood", {
  skip_on_cran()
  set.seed(3)
  n <- 500
  x <- rnorm(n)
  d <- data.frame(y = rpois(n, exp(-0.3 + 0.4 * x)), x = x)
  ll_pois <- as.numeric(logLik(glm(y ~ x, poisson, data = d)))
  for (meth in c("NM", "BFGS", "NR", "BHHH")) {
    f <- suppressWarnings(countreg(y ~ x, data = d, family = "PIG1",
                                   method = meth, max.iters = 2000))
    ll <- as.numeric(logLik(f))
    expect_lte(ll, 0)
    # The PIG1 likelihood is at least the Poisson one (its limit) and not
    # far above it for Poisson data
    expect_lt(ll - ll_pois, 5)
  }
})

test_that("count quantiles are correct for tiny probabilities", {
  expect_identical(qpinvgaus(1e-15, mu = 100, eta = 1e-6),
                   qpois(1e-15, 100))
  expect_identical(qpinvgaus(1e-10, mu = 100, eta = 1e-6),
                   qpois(1e-10, 100))
})

test_that("Sichel variance, density and quantiles are finite at extremes", {
  v <- .fcr_family_registry$SI$variance
  vv <- v(mu = 2, alpha = exp(1), sigma = c(1e-8, 0.001, 0.5))
  expect_true(all(is.finite(vv)))
  expect_equal(vv[1], 2, tolerance = 1e-6)        # Poisson limit
  # Variance against the moments of the PMF
  y <- 0:3000
  p <- dsichel(y, mu = 3, sigma = 0.8, gamma = 1.3)
  m1 <- sum(p * y)
  expect_equal(v(3, exp(1.3), 0.8), sum(p * (y - m1)^2), tolerance = 1e-8)

  expect_true(is.finite(dsichel(500, mu = 500, sigma = 1, gamma = 1,
                                log = TRUE)))
  expect_no_warning(q <- qsichel(0.999999, mu = 5000, sigma = 5,
                                 gamma = -0.5))
  expect_gte(psichel(q, 5000, 5, -0.5), 0.999999 - 1e-12)
  expect_lt(psichel(q - 1, 5000, 5, -0.5), 0.999999)
  # The order recurrence agrees with direct evaluation
  for (g in c(-3.7, -0.5, 0, 2.3)) {
    expect_equal(.fcr_sichel_logpmf_seq(40, 4, 0.7, g),
                 dsichel(0:40, 4, 0.7, g, log = TRUE), tolerance = 1e-12)
  }
})

test_that("the NB1/NB2 analytic scores are finite and accurate near alpha = 0", {
  data("washington_roads", package = "flexCountReg")
  d <- washington_roads
  X <- model.matrix(~ lnaadt, d)
  y <- d$Total_crashes
  b <- c(-5, 0.6)
  mu <- as.vector(exp(X %*% b))
  n <- nrow(d)
  ref <- function(la) {
    r <- 1 / exp(la)
    sum(vapply(seq_along(y), function(i) {
      s <- if (y[i] > 0) sum(1 / (r + 0:(y[i] - 1))) else 0
      -r * (s - log1p(mu[i] / r) + (mu[i] - y[i]) / (r + mu[i]))
    }, numeric(1)))
  }
  for (la in c(-2, -10, -18)) {
    G <- .fcr_analytic_grad("NB2", c(b, la), y, X, rep(0, n),
                            matrix(1, n, 1), rep(1, n))
    expect_equal(sum(G[, 3]), ref(la), tolerance = 1e-6)
  }
  for (fam in c("NB1", "NB2")) {
    G <- .fcr_analytic_grad(fam, c(b, -800), y, X, rep(0, n),
                            matrix(1, n, 1), rep(1, n))
    expect_false(anyNA(G))
    # Poisson limit of the mean score
    expect_equal(G[, 1], y - mu, tolerance = 1e-10, ignore_attr = TRUE)
  }
  # NB1 ln(alpha) score against its leading-order series
  a <- exp(-12)
  G <- .fcr_analytic_grad("NB1", c(b, -12), y, X, rep(0, n),
                          matrix(1, n, 1), rep(1, n))
  expect_equal(sum(G[, 3]), sum(a * ((y - mu)^2 - y) / (2 * mu)),
               tolerance = 1e-4)
})

test_that("vectorized digamma and lgamma differences match the finite sums", {
  set.seed(1)
  x <- 10^runif(400, -3, 14)
  d <- rpois(400, 10^runif(400, -1, 2.5))
  dg <- mapply(function(a, k) if (k == 0) 0 else sum(1 / (a + 0:(k - 1))),
               x, d)
  lg <- mapply(function(a, k) if (k == 0) 0 else sum(log1p((0:(k - 1)) / a)),
               x, d)
  expect_equal(.fcr_digamma_diff(x, d), dg, tolerance = 1e-11)
  expect_lt(max(abs(.fcr_lgamma_diff_rel(x, d) - lg)), 1e-8)
})

test_that("a missing ancillary covariate gives NA probabilities", {
  data("washington_roads", package = "flexCountReg")
  d <- washington_roads
  set.seed(2)
  d$z <- rnorm(nrow(d))
  m <- countreg(Total_crashes ~ lnaadt + lnlength, data = d, family = "NB2",
                dis_param_formula_1 = ~ z, method = "BFGS")
  nd <- d[1:3, ]
  nd$z[2] <- NA
  P <- predict(m, newdata = nd, type = "prob", max_count = 4,
               overflow = TRUE)
  expect_true(all(is.na(P[2, ])))
  expect_equal(unname(rowSums(P[-2, ])), c(1, 1), tolerance = 1e-12)
})

test_that("the comparison baseline keeps a model's missing intercept", {
  data("washington_roads", package = "flexCountReg")
  d <- washington_roads
  m0 <- countreg(Total_crashes ~ 0 + lnaadt + lnlength, data = d,
                 family = "NBP", method = "BFGS")
  nb2 <- countreg(Total_crashes ~ 0 + lnaadt + lnlength, data = d,
                  family = "NB2", method = "BFGS")
  r <- suppressWarnings(suppressMessages(
    regCompTest(m0, basemodel = "NB2", variables = TRUE)))
  expect_equal(r$LLbase, as.numeric(logLik(nb2)), tolerance = 1e-4)
})

test_that("Empirical Bayes predictions match panels by id, not row order", {
  skip_on_cran()
  data("washington_roads", package = "flexCountReg")
  d <- washington_roads
  gvec <- rep(seq_len(nrow(d) %/% 3 + 1), each = 3)[seq_len(nrow(d))]
  rp <- suppressWarnings(countreg_rp(
    Total_crashes ~ lnaadt + lnlength, rpar_formula = ~ -1 + speed50,
    data = d, family = "NB2", ndraws = 50, method = "BFGS",
    group_var = gvec))
  expect_identical(rp$model$group_var, ".fcr_panel_id")
  p1 <- predict(rp, method = "Individual")
  set.seed(4)
  o <- sample(nrow(d))
  nd <- rp$data[o, ]
  p2 <- predict(rp, newdata = nd, method = "Individual")
  expect_equal(as.vector(p2), as.vector(p1)[o], tolerance = 1e-12)
  expect_equal(attr(p2, "eb_coefficients"),
               attr(p1, "eb_coefficients")[o, , drop = FALSE],
               tolerance = 1e-12)
  # Data without the stored ids cannot be matched: an error, not a guess
  expect_error(predict(rp, newdata = d[o, ], method = "Individual"),
               "group_var")
})

test_that("quantile residuals stay finite far in the upper tail", {
  dd <- data.frame(y = c(rep(0, 24), 25))
  mp <- countreg(y ~ 1, data = dd, family = "Poisson", method = "BFGS")
  set.seed(1)
  r <- residuals(mp, type = "quantile")
  expect_true(all(is.finite(r)))
  expect_gt(r[25], 8)
})

test_that("default marginal effects include dispersion variables of hurdles", {
  data("washington_roads", package = "flexCountReg")
  mh <- countreg(Total_crashes ~ lnaadt + lnlength,
                 data = washington_roads, family = "NB2",
                 dis_param_formula_1 = ~ speed50, zero_formula = ~ 1,
                 zero_type = "hurdle", method = "BFGS")
  me <- margEffTable(mh, digits = 6)
  expect_true("speed50" %in% me$term)
  expect_true(.fcr_mean_uses_dispersion("PLN"))
  expect_false(.fcr_mean_uses_dispersion("NB2"))
})

test_that("delta-method marginal-effect SEs match closed forms (NB2)", {
  data("washington_roads", package = "flexCountReg")
  d <- washington_roads
  m <- countreg(Total_crashes ~ lnaadt + lnlength + speed50, data = d,
                family = "NB2", method = "BFGS")
  b <- coef(m)
  V <- vcov(m)
  X <- model.matrix(~ lnaadt + lnlength + speed50, d)
  mu <- as.vector(exp(X %*% b[1:4]))
  g <- numeric(length(b))
  g[1:4] <- b[2] * colMeans(mu * X)
  g[2] <- g[2] + mean(mu)
  me <- margEffTable(m, vars = "lnaadt", measure = "ame", digits = 10)
  expect_equal(me$estimate, unname(b[2] * mean(mu)), tolerance = 1e-6)
  expect_equal(me$std_error, sqrt(drop(t(g) %*% V %*% g)), tolerance = 1e-4)
})

test_that("regCompTable(coefs = FALSE) leaves out the coefficients", {
  data("washington_roads", package = "flexCountReg")
  d <- washington_roads
  a <- countreg(Total_crashes ~ lnaadt + lnlength, data = d, family = "NB2",
                method = "BFGS")
  b <- countreg(Total_crashes ~ lnaadt + lnlength, data = d,
                family = "Poisson", method = "BFGS")
  tb <- regCompTable(list(a = a, b = b), coefs = FALSE, RSquare = FALSE)
  expect_false(any(c("(Intercept)", "lnaadt", "ln(alpha)") %in% tb$Parameter))
  expect_true(all(c("N Obs.", "LL", "AIC", "BIC") %in% tb$Parameter))
})

test_that("keep_data = FALSE stores neither the data nor closures over it", {
  data("washington_roads", package = "flexCountReg")
  d <- washington_roads
  # A formula keeps its environment (here it would be this test's frame,
  # with the data in it), so it is created in the global environment, as at
  # the console
  f <- Total_crashes ~ lnaadt + lnlength
  environment(f) <- globalenv()
  m1 <- countreg(f, data = d, family = "NB2", method = "BFGS")
  m0 <- countreg(f, data = d, family = "NB2", method = "BFGS",
                 keep_data = FALSE)
  expect_null(m0$data)
  expect_false(any(vapply(m0$model, is.function, logical(1))))
  expect_false(any(vapply(m1$model, is.function, logical(1))))
  expect_true(isTRUE(m0$model$lr_sample$fingerprinted))
  expect_false(any(vapply(m0$model$lr_sample$vars, length, integer(1)) ==
                     nrow(d)))
  expect_lt(length(serialize(m0, NULL)), length(serialize(d, NULL)))
  # The fingerprinted sample still recognises the same and other samples
  expect_identical(.fcr_same_sample(m0$model$lr_sample,
                                    m1$model$lr_sample), "")
  d2 <- d
  d2$lnaadt[1] <- d2$lnaadt[1] + 0.1
  m2 <- countreg(Total_crashes ~ lnaadt + lnlength, data = d2,
                 family = "NB2", method = "BFGS", keep_data = FALSE)
  expect_match(.fcr_same_sample(m0$model$lr_sample, m2$model$lr_sample),
               "lnaadt")
})

test_that("robust SEs with sampling weights do not depend on their scale", {
  data("washington_roads", package = "flexCountReg")
  d <- washington_roads
  set.seed(5)
  d$w <- runif(nrow(d), 0.5, 2)
  d$w100 <- 100 * d$w
  fit <- function(w, type) {
    countreg(Total_crashes ~ lnaadt + lnlength, data = d, family = "NB2",
             method = "BFGS", weights = w, weight_type = type)
  }
  s1 <- fit("w", "sampling")
  s100 <- fit("w100", "sampling")
  expect_identical(s1$model$vcov_type, "robust")      # the default
  expect_equal(coef(s1), coef(s100), tolerance = 1e-4)
  expect_equal(s1$model$se, s100$model$se, tolerance = 1e-3)
  expect_equal(s1$model$n_bic, nrow(d))
  # Frequency weights keep their meaning: 100 copies shrink SEs 10-fold
  f1 <- countreg(Total_crashes ~ lnaadt + lnlength, data = d, family = "NB2",
                 method = "BFGS", weights = "w", stderr = "robust")
  f100 <- countreg(Total_crashes ~ lnaadt + lnlength, data = d,
                   family = "NB2", method = "BFGS", weights = "w100",
                   stderr = "robust")
  expect_equal(f1$model$se / f100$model$se, rep(10, 4), tolerance = 1e-3,
               ignore_attr = TRUE)
  # With unit sampling weights the sandwich is the usual one times n / (n - 1)
  d$one <- 1
  u <- countreg(Total_crashes ~ lnaadt + lnlength, data = d, family = "NB2",
                method = "BFGS", weights = "one", weight_type = "sampling")
  r <- countreg(Total_crashes ~ lnaadt + lnlength, data = d, family = "NB2",
                method = "BFGS", stderr = "robust")
  n <- nrow(d)
  expect_equal(u$model$vcov, r$model$vcov * n / (n - 1), tolerance = 1e-6)
  expect_error(countreg(Total_crashes ~ lnaadt, data = d, family = "NB2",
                        weights = "w", weight_type = "sampling",
                        stderr = "boot", bootstraps = 5),
               "sampling")
  expect_output(print(summary(s1)), "pseudo-likelihood")
  expect_error(margEffTable(s1, se = "bootstrap", bootstraps = 2),
               "sampling")
})

test_that("the inverse-Gaussian sampler is positive for small shapes", {
  set.seed(1)
  x <- rinvgaus(1e4, mean = 1, shape = 1e-9)
  expect_true(all(x > 0 & is.finite(x)))
  x <- rinvgaus(5000, mean = 1, shape = 2)
  cdf <- function(q) stats::pnorm(sqrt(2 / q) * (q - 1)) +
    exp(4) * stats::pnorm(-sqrt(2 / q) * (q + 1))
  expect_gt(suppressWarnings(stats::ks.test(x, cdf)$p.value), 0.001)
})
