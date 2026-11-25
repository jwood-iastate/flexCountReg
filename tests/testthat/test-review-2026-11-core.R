# Regression tests for the November 2026 review: Poisson-Weibull quadrature
# for extreme shapes, DPO and GAMMACOUNT over their documented shape domain,
# countreg(bootstraps =) without stderr, and moments of user families.

pw_fam <- function() flexCountReg:::.fcr_get_family("PW")

# Reference P(Y = y) for the mean-parameterized Poisson-Weibull by
# integrate() over w = log(E), E ~ Exp(1), lambda = c E^(1 / alpha)
pw_ref_log <- function(y, mu, alpha) {
  logc <- log(mu) - lgamma(1 + 1 / alpha)
  lg <- function(w) {
    v <- stats::dpois(y, exp(logc + w / alpha), log = TRUE) + w - exp(w)
    v[is.na(v)] <- -Inf
    v
  }
  grid <- seq(-60, 6, length.out = 20001)
  hv <- lg(grid)
  j <- which.max(hv)
  hm <- hv[j]
  br <- sort(unique(c(-60, 6, pmin(pmax(grid[j] + c(-30, -10, -3, -1, -0.3,
    -0.1, -0.03, 0, 0.03, 0.1, 0.3, 1, 3, 10, 30), -60), 6))))
  tot <- 0
  for (i in seq_len(length(br) - 1)) {
    tot <- tot + stats::integrate(function(w) exp(lg(w) - hm), br[i],
                                  br[i + 1], rel.tol = 1e-12, abs.tol = 0,
                                  subdivisions = 2000)$value
  }
  hm + log(tot)
}

test_that("PW quadrature is finite and normalized for extreme shapes", {
  f <- pw_fam()
  expect_equal(f$dens(0, 1e-4, 40, NULL, engine = "quadrature"),
               0.9999, tolerance = 1e-6)
  # alpha = 100: nearly Poisson; the tail probability is close to dpois
  p20 <- f$dens(20, 0.5, 100, NULL, engine = "quadrature")
  expect_true(is.finite(p20))
  expect_equal(p20 / dpois(20, 0.5), 1, tolerance = 0.05)
  for (a in c(1e-4, 1e-3, 35, 50, 100, 200)) {
    for (mu in c(1e-4, 0.5, 5)) {
      y <- 0:400
      lp <- f$dens(y, rep(mu, length(y)), rep(a, length(y)), NULL,
                   engine = "quadrature", log = TRUE)
      expect_false(anyNA(lp))
      expect_true(all(lp <= 0))
      expect_equal(sum(exp(lp)), 1, tolerance = 1e-9,
                   label = sprintf("PMF sum at alpha %g, mu %g", a, mu))
    }
  }
})

test_that("PW quadrature matches integrate() over a grid of shapes", {
  f <- pw_fam()
  for (a in c(1e-3, 0.1, 1, 20, 100)) {
    for (mu in c(1e-3, 0.5, 50)) {
      y <- c(0, 1, 5, 30)
      lp <- f$dens(y, rep(mu, 4), rep(a, 4), NULL, engine = "quadrature",
                   log = TRUE)
      ref <- vapply(y, pw_ref_log, 1, mu = mu, alpha = a)
      use <- ref > -500
      expect_equal(lp[use], ref[use], tolerance = 1e-8,
                   label = sprintf("alpha %g, mu %g", a, mu))
    }
  }
  # Small alpha, y = 0: the integrand has a steep edge inside the bulk of
  # the mixing density (the fixed rule had a relative error of ~1e-4 in
  # 1 - P(Y = 0)); reference by integrating over the Weibull quantile
  c0 <- 0.5 / gamma(11)
  g <- function(u) exp(-c0 * (-log1p(-u))^10)
  br <- c(0, 10^seq(-12, 0, by = 0.05))
  p0 <- sum(vapply(seq_len(length(br) - 1), function(i)
    integrate(g, br[i], br[i + 1], rel.tol = 1e-13, abs.tol = 0)$value, 1))
  expect_equal(f$dens(0, 0.5, 0.1, NULL, engine = "quadrature", log = TRUE),
               log(p0), tolerance = 1e-9)
})

test_that("dpoisweibull and ppoisweibull handle extreme shapes", {
  expect_equal(dpoisweibull(0, mean_value = 1e-4, alpha = 40, sigma = 1),
               0.9999, tolerance = 1e-6)
  lp <- dpoisweibull(0:3, lambda = 1e-4, alpha = 50, sigma = 1, log = TRUE)
  expect_true(all(is.finite(lp) & lp <= 0))
  expect_equal(sum(dpoisweibull(0:200, lambda = 1e-4, alpha = 50,
                                sigma = 1)), 1, tolerance = 1e-10)
  expect_equal(ppoisweibull(5, lambda = 1e-4, alpha = 40, sigma = 1), 1,
               tolerance = 1e-10)
  # Shapes whose Gamma(1 + 1 / alpha) overflows (lambda parameterization)
  p <- expect_silent(dpoisweibull(0:2, lambda = 1, alpha = 1e-3, sigma = 1))
  expect_true(all(is.finite(p) & p > 0 & p < 1))
  expect_equal(ppoisweibull(2, lambda = 1, alpha = 1e-3, sigma = 1),
               sum(p), tolerance = 1e-12)
  # lambda and mean parameterizations agree where both are finite
  a <- 0.02
  m <- 2 * 1.5 * gamma(1 + 1 / a)
  expect_equal(dpoisweibull(0:20, lambda = 2, alpha = a, sigma = 1.5,
                            log = TRUE),
               pw_fam()$dens(0:20, rep(m, 21), rep(a, 21), NULL,
                             engine = "quadrature", log = TRUE),
               tolerance = 1e-10)
  # ordinary shapes are unchanged
  expect_equal(ppoisweibull(c(0, 2, 10), mean_value = 3, alpha = 2,
                            sigma = 1.5),
               c(0.1219522, 0.4863052, 0.9929785), tolerance = 1e-6)
  m2 <- 3 * 1.5 * gamma(1.5)
  expect_equal(ppoisweibull(c(0, 2, 10), lambda = 3, alpha = 2, sigma = 1.5),
               cumsum(pw_fam()$dens(0:10, rep(m2, 11), rep(2, 11), NULL,
                                    engine = "quadrature"))[c(1, 3, 11)],
               tolerance = 1e-12)
  expect_equal(qpoisweibull(c(0.1, 0.5, 0.9), lambda = 3, alpha = 2,
                            sigma = 1.5), c(1, 3, 8))
  expect_warning(r <- dpoisweibull(1, lambda = 1, alpha = -1, sigma = 1),
                 "NaN")
  expect_true(is.nan(r))
})

test_that("DPO solves the mean equation over the documented shape domain", {
  f <- flexCountReg:::.fcr_get_family("DPO")
  # was -Inf below phi = 0.003 (mean 0.5)
  expect_equal(f$dens(0, 0.5, 0.003, NULL, log = TRUE), -0.3127,
               tolerance = 1e-3)
  # (each row solves the mean equation; phi = 1e-6 takes ~0.1 s a row)
  cases <- list(c(1e-6, 0.5), c(1e-4, 0.5), c(1e-4, 5), c(0.002, 0.5),
                c(0.002, 5))
  for (cs in cases) {
    phi <- cs[1]; mu <- cs[2]
    ymax <- if (mu < 1) 80 else 300
    y <- 0:ymax
    p <- f$dens(y, rep(mu, length(y)), rep(phi, length(y)), NULL)
    up <- exp(f$cdf(ymax, mu, phi, NULL, lower.tail = FALSE))
    expect_equal(sum(p) + up, 1, tolerance = 1e-12)
    expect_lt(up, 1e-12)
    expect_equal(sum(y * p), mu, tolerance = 1e-10)
    expect_equal(f$variance(mu, phi, NULL), sum((y - mu)^2 * p),
                 tolerance = 1e-8)
  }
  # beyond the window budget (mu / phi far above 1e10) the density is 0
  expect_equal(f$dens(0, 1e6, 1e-6, NULL), 0)
})

test_that("GAMMACOUNT solves the mean equation for small shapes", {
  f <- flexCountReg:::.fcr_get_family("GAMMACOUNT")
  for (a in c(1e-6, 1e-3, 0.01)) {
    for (mu in c(0.5, 50)) {
      y <- 0:3000
      lp <- f$dens(y, rep(mu, length(y)), rep(a, length(y)), NULL,
                   log = TRUE)
      expect_true(all(is.finite(lp)))
      p <- exp(lp)
      expect_equal(sum(p), 1, tolerance = 1e-12)
      expect_equal(sum(y * p), mu, tolerance = 1e-10)
      expect_equal(f$variance(mu, a, NULL), sum((y - mu)^2 * p),
                   tolerance = 1e-8)
      expect_equal(exp(f$cdf(3, mu, a, NULL)), sum(p[1:4]),
                   tolerance = 1e-12)
    }
  }
  # Large means (the direct sum previously needed more than 5000 terms)
  expect_equal(exp(flexCountReg:::gc_mp_cpp(1e5, 1, 0L)), 1e5,
               tolerance = 1e-10)
  v <- flexCountReg:::gc_mp_cpp(c(1e5, 1e7), c(0.5, 0.01), 0L)
  expect_false(anyNA(v))
  m <- flexCountReg:::.gcount_moments(c(0.5, 0.01), exp(v))
  expect_equal(m$mean, c(1e5, 1e7), tolerance = 1e-9)
  expect_equal(flexCountReg:::gc_mp_cpp(c(1e5, 1e7), c(0.5, 0.01), 1L),
               m$var, tolerance = 1e-7)
  # draws use the log rate (the rate underflows at small shapes)
  set.seed(1)
  r <- f$rand(4000, 2, 1e-3, NULL)
  expect_equal(mean(r), 2, tolerance = 0.1)
})

test_that("countreg() uses bootstrap SEs when only `bootstraps` is given", {
  data("washington_roads", package = "flexCountReg")
  set.seed(3)
  m <- countreg(Total_crashes ~ lnaadt, washington_roads, family = "POISSON",
                method = "BFGS", bootstraps = 3)
  expect_equal(m$model$vcov_type, "bootstrap")
  expect_equal(m$model$n_boot_ok, 3)
})

test_that("moments of user families stop for unnormalized PMFs", {
  nb <- new_count_family("nbu",
    logpmf = function(y, mu, par)
      dnbinom(y, size = 1 / par$alpha, mu = mu, log = TRUE) + log1p(-1e-9),
    params = list(alpha = count_param("log", start = -1)))
  mu <- c(0.5, 2, 10)
  t <- system.time(m <- flexCountReg:::.fcr_user_moments(
    nb, mu, list(alpha = rep(0.5, 3))))[3]
  expect_lt(t, 2)
  expect_equal(m$mean, mu, tolerance = 1e-8)
  expect_equal(m$var, mu + 0.5 * mu^2, tolerance = 1e-8)
  # rows with no probability (invalid parameters) end quickly with NaN
  nbv <- new_count_family("nbv",
    logpmf = function(y, mu, par)
      dnbinom(y, size = 1 / par$alpha, mu = mu, log = TRUE),
    params = list(alpha = count_param("log", start = -1)),
    valid = function(mu, par) par$alpha < 5)
  t <- system.time(m <- flexCountReg:::.fcr_user_moments(
    nbv, c(1, 1), list(alpha = c(0.5, 10))))[3]
  expect_lt(t, 2)
  expect_equal(m$var, c(1.5, NaN))
  # a very heavy tail (infinite variance) warns at max_y
  heavy <- new_count_family("heavy",
    logpmf = function(y, mu, par) {
      # P(Y = y) proportional to (y + 1)^-2.5, ignoring mu
      -2.5 * log(y + 1) - log(1.341487257250917)
    },
    params = list(alpha = count_param("log", start = 0)),
    mean_is_mu = FALSE)
  expect_warning(flexCountReg:::.fcr_user_moments(heavy, 1, list(alpha = 1),
                                                  max_y = 1e4),
                 "without converging")
})
