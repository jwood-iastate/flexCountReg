# Tests for the count distributions added in 0.1.2: discrete Weibull,
# double Poisson, Delaporte, generalized Hermite, hyper-Poisson, gamma count,
# Polya-Aeppli and the negative binomial-generalized exponential mixture. Reference values are computed independently of the
# package code (brute-force convolutions, closed forms, base R).

fams <- list(
  dweib = list(d = ddweib, p = pdweib, q = qdweib, r = rdweib,
               args = list(theta = 0.8, shape = 1.3),
               bad = list(theta = 1.2, shape = 1)),
  dpo = list(d = ddpo, p = pdpo, q = qdpo, r = rdpo,
             args = list(lambda = 4, phi = 0.4),
             bad = list(lambda = 4, phi = -1)),
  delaporte = list(d = ddelaporte, p = pdelaporte, q = qdelaporte,
                   r = rdelaporte, args = list(mu = 5, rho = 0.6, k = 1.2),
                   bad = list(mu = 5, rho = 1.5, k = 1)),
  ghermite = list(d = dghermite, p = pghermite, q = qghermite, r = rghermite,
                  args = list(mu = 3, rho = 0.5, m = 3),
                  bad = list(mu = 3, rho = 0.5, m = 2.5)),
  hpois = list(d = dhpois, p = phpois, q = qhpois, r = rhpois,
               args = list(lambda = 4, gamma = 0.3),
               bad = list(lambda = -4, gamma = 1)),
  gcount = list(d = dgcount, p = pgcount, q = qgcount, r = rgcount,
                args = list(shape = 0.6, rate = 2, time = 2),
                bad = list(shape = 0, rate = 2, time = 1)),
  paeppli = list(d = dpaeppli, p = ppaeppli, q = qpaeppli, r = rpaeppli,
                 args = list(lambda = 2, prob = 0.4),
                 bad = list(lambda = 2, prob = 1)),
  nbge = list(d = dnbge, p = pnbge, q = qnbge, r = rnbge,
              args = list(mu = 3, shape = 1.5, alpha = 0.3),
              bad = list(mu = 3, shape = 1.5, alpha = -0.3))
)

call_f <- function(f, x, args, ...) do.call(f, c(list(x), args, list(...)))

# ---- special cases ----------------------------------------------------------

test_that("special cases reduce to base R distributions", {
  y <- 0:60
  expect_equal(ddweib(y, theta = 0.6, shape = 1), dgeom(y, prob = 0.4),
               tolerance = 1e-12)
  for (lambda in c(0.05, 0.7, 4, 25)) {
    expect_equal(ddpo(0:200, lambda = lambda, phi = 1, log = TRUE),
                 dpois(0:200, lambda = lambda, log = TRUE), tolerance = 1e-9)
    expect_equal(dhpois(y, lambda = lambda, gamma = 1), dpois(y, lambda),
                 tolerance = 1e-12)
    expect_equal(dpaeppli(y, lambda = lambda, prob = 0), dpois(y, lambda),
                 tolerance = 1e-12)
    expect_equal(dghermite(y, mu = lambda, rho = 0), dpois(y, lambda),
                 tolerance = 1e-12)
    expect_equal(dghermite(y, mu = lambda, rho = 0.7, m = 1),
                 dpois(y, lambda), tolerance = 1e-12)
    expect_equal(ddelaporte(y, mu = lambda, rho = 0, k = 2),
                 dpois(y, lambda), tolerance = 1e-12)
    expect_equal(ddelaporte(y, mu = lambda, rho = 1, k = 2),
                 dnbinom(y, size = 2, mu = lambda), tolerance = 1e-12)
  }
  expect_equal(dgcount(y, shape = 1, rate = 1.5, time = 2), dpois(y, 3),
               tolerance = 1e-12)
  expect_equal(pgcount(y, shape = 1, rate = 1.5, time = 2,
                       lower.tail = FALSE),
               ppois(y, 3, lower.tail = FALSE), tolerance = 1e-12)
  # Hermite with rho = 1: all mass on multiples of m
  h <- dghermite(0:30, mu = 2, rho = 1, m = 3)
  expect_equal(h[(0:30) %% 3 != 0], rep(0, 20))
  expect_equal(h[(0:30) %% 3 == 0], dpois(0:10, 2 / 3), tolerance = 1e-12)
  expect_equal(dghermite(5, 2, 1, 3, log = TRUE), -Inf)
  # NBGE with alpha = 0 is the Poisson-generalized exponential
  expect_equal(dnbge(y, mu = 2, shape = 1.5, alpha = 0),
               dpge(y, mean = 2, shape = 1.5, scale = 1), tolerance = 1e-12)
  # ... and is continuous as alpha -> 0
  expect_equal(dnbge(0:5, 2, 1.5, 1e-9), dnbge(0:5, 2, 1.5, 0),
               tolerance = 1e-8)
  # Delaporte approaches Poisson as k grows
  expect_equal(ddelaporte(y, 3, 0.5, 1e9), dpois(y, 3), tolerance = 1e-6)
})

# ---- independent references -------------------------------------------------

test_that("probabilities match independent calculations", {
  y <- 0:80
  # Discrete Weibull: theta^(y^b) - theta^((y+1)^b)
  th <- 0.85; b <- 0.7
  expect_equal(ddweib(y, th, b), th^(y^b) - th^((y + 1)^b),
               tolerance = 1e-10)
  expect_equal(pdweib(y, th, b), 1 - th^((y + 1)^b), tolerance = 1e-12)

  # Delaporte: brute-force convolution
  mu <- 6; rho <- 0.4; k <- 0.8
  pa <- dpois(0:80, (1 - rho) * mu)
  pb <- dnbinom(0:80, size = k, mu = rho * mu)
  conv <- vapply(y, function(v) sum(pa[1:(v + 1)] * pb[(v + 1):1]), 0)
  expect_equal(ddelaporte(y, mu, rho, k), conv, tolerance = 1e-12)

  # Hermite of order 3: brute-force convolution
  mu <- 4; rho <- 0.6; m <- 3
  l1 <- (1 - rho) * mu; lm <- rho * mu / m
  herm <- vapply(y, function(v) {
    j <- 0:(v %/% m)
    sum(dpois(j, lm) * dpois(v - m * j, l1))
  }, 0)
  expect_equal(dghermite(y, mu, rho, m), herm, tolerance = 1e-12)

  # Polya-Aeppli closed-form sum
  lam <- 2.5; q <- 0.35
  pae <- vapply(y, function(v) {
    if (v == 0) return(exp(-lam))
    j <- 1:v
    sum(exp(-lam + j * log(lam) - lgamma(j + 1) + lchoose(v - 1, j - 1) +
              (v - j) * log(q) + j * log1p(-q)))
  }, 0)
  expect_equal(dpaeppli(y, lam, q), pae, tolerance = 1e-12)

  # Double Poisson: weights normalized over a long grid
  lam <- 3.5; phi <- 0.45
  g <- 0:3000
  lw <- phi * g * log(lam) + (1 - phi) * ifelse(g > 0, g * log(g), 0) +
    (phi - 1) * g - lgamma(g + 1)
  w <- exp(lw - max(lw))
  expect_equal(ddpo(y, lam, phi), (w / sum(w))[y + 1], tolerance = 1e-12)

  # Hyper-Poisson: normalizing constant from gsl's 1F1
  lam <- 3; gam <- 2.5
  Z <- gsl::hyperg_1F1(1, gam, lam)
  ref <- exp(y * log(lam) - (lgamma(gam + y) - lgamma(gam))) / Z
  expect_equal(dhpois(y, lam, gam), ref, tolerance = 1e-12)

  # Gamma count: gamma distribution function differences
  a <- 1.7; s <- 5
  G <- function(c) ifelse(c == 0, 1, pgamma(s, c))
  expect_equal(dgcount(y, a, rate = 2.5, time = 2),
               G(y * a) - G((y + 1) * a), tolerance = 1e-9)

  # NBGE: integral over t = log(x) of the NB2 probability times the
  # density of t, by stats::integrate on a partition of the line
  nbge_ref <- function(v, mu, k, a) {
    A <- digamma(k + 1) - digamma(1)
    f <- function(t) {
      x <- exp(t)
      exp(dnbinom(v, size = 1 / a, mu = mu * x / A, log = TRUE) + log(k) +
            t - x + (k - 1) * log(-expm1(-x)))
    }
    br <- c(-60, -20, -8, -4, -2, -1, 0, 1, 2, 3, 4, 6)
    sum(vapply(seq_len(length(br) - 1), function(i)
      stats::integrate(f, br[i], br[i + 1], rel.tol = 1e-13, abs.tol = 0,
                       subdivisions = 2000)$value, 0))
  }
  yy <- c(0:10, 25, 40)
  for (pr in list(c(2, 1.5, 0.5), c(10, 0.4, 1.2), c(30, 2, 2))) {
    ref <- vapply(yy, nbge_ref, 0, mu = pr[1], k = pr[2], a = pr[3])
    expect_equal(dnbge(yy, pr[1], pr[2], pr[3]), ref, tolerance = 1e-9)
  }
})

test_that("documented moment identities hold", {
  # Hyper-Poisson: E[Y] = lambda - (gamma - 1)(1 - P(Y = 0))
  for (pr in list(c(3, 2), c(3, 0.4), c(15, 7.5))) {
    m <- flexCountReg:::.hpois_moments(pr[1], pr[2])$mean
    expect_equal(m, pr[1] - (pr[2] - 1) * (1 - dhpois(0, pr[1], pr[2])),
                 tolerance = 1e-10)
  }
  # Moment helpers agree with sums over the probabilities
  y <- 0:2000
  mom <- function(pm) {
    m <- sum(y * pm)
    c(m, sum(y^2 * pm) - m^2)
  }
  chk <- function(mo, pm) {
    expect_equal(c(mo$mean, mo$var), mom(pm), tolerance = 1e-9)
  }
  chk(flexCountReg:::.dweib_moments(0.8, 1.3), ddweib(y, 0.8, 1.3))
  chk(flexCountReg:::.dpo_moments(4, 0.4), ddpo(y, 4, 0.4))
  chk(flexCountReg:::.dpo_moments(4, 3), ddpo(y, 4, 3))
  chk(flexCountReg:::.hpois_moments(4, 0.3), dhpois(y, 4, 0.3))
  chk(flexCountReg:::.gcount_moments(0.6, 2, 2), dgcount(y, 0.6, 2, 2))
  chk(flexCountReg:::.delap_moments(5, 0.6, 1.2),
      ddelaporte(y, 5, 0.6, 1.2))
  chk(flexCountReg:::.ghermite_moments(3, 0.5, 3), dghermite(y, 3, 0.5, 3))
  chk(flexCountReg:::.pa_moments(2, 0.4), dpaeppli(y, 2, 0.4))
  y <- 0:6000
  chk(flexCountReg:::.nbge_moments(2, 1.5, 0.5), dnbge(y, 2, 1.5, 0.5))
})

test_that("discrete Weibull moments are accurate for heavy tails", {
  # shape = 0.35: the direct sum needs ~2e5 terms, beyond the 1e5 cap, so
  # the Euler-Maclaurin remainder is used; compare with a direct sum.
  th <- 0.5; b <- 0.35
  k <- 1:3e6
  f <- exp(log(th) * k^b)
  m1 <- sum(f)
  m2 <- sum((2 * k - 1) * f)
  mo <- flexCountReg:::.dweib_moments(th, b)
  expect_equal(mo$mean, m1, tolerance = 1e-9)
  expect_equal(mo$var, m2 - m1^2, tolerance = 1e-8)
})

test_that("the Polya-Aeppli distribution is Poisson-Tweedie with power 1.5", {
  # PGF of Y | L ~ Poisson(L), L compound Poisson with rate 2 sqrt(mu) / phi
  # and exponential jumps of mean phi sqrt(mu) / 2
  mu <- 3; phi <- 0.8
  rate <- 2 * sqrt(mu) / phi
  theta <- phi * sqrt(mu) / 2
  pgf_pt <- function(z) exp(-rate * theta * (1 - z) / (1 + theta * (1 - z)))
  pa <- flexCountReg:::.pt15_to_paeppli(mu, phi)
  pm <- dpaeppli(0:500, pa$lambda, pa$prob)
  for (z in c(0, 0.3, 0.7, 0.95)) {
    expect_equal(sum(z^(0:500) * pm), pgf_pt(z), tolerance = 1e-12)
  }
  mo <- flexCountReg:::.pa_moments(pa$lambda, pa$prob)
  expect_equal(c(mo$mean, mo$var), c(mu, mu + phi * mu^1.5))
})

# ---- shared contract ---------------------------------------------------------

test_that("pmf, cdf, survival and log scales are consistent", {
  y <- 0:400
  for (nm in names(fams)) {
    f <- fams[[nm]]
    pm <- call_f(f$d, y, f$args)
    cdf <- call_f(f$p, y, f$args)
    sf <- call_f(f$p, y, f$args, lower.tail = FALSE)
    expect_equal(sum(pm), 1, tolerance = 1e-12, label = nm)
    expect_equal(cumsum(pm), cdf, tolerance = 1e-12, label = nm)
    expect_equal(cdf + sf, rep(1, length(y)), tolerance = 1e-12, label = nm)
    big <- pm > 1e-290   # log(pm) loses precision for subnormal pm
    expect_equal(call_f(f$d, y, f$args, log = TRUE)[big], log(pm[big]),
                 label = nm)
    expect_equal(call_f(f$p, y[1:30], f$args, log.p = TRUE), log(cdf[1:30]),
                 label = nm)
    # non-integer q is floored
    expect_equal(call_f(f$p, 2.7, f$args), cdf[3], label = nm)
  }
})

test_that("upper tails keep their relative precision", {
  for (nm in names(fams)) {
    f <- fams[[nm]]
    q <- if (nm == "nbge") 400 else 150   # NBGE has a heavier tail
    ls <- call_f(f$p, q, f$args, lower.tail = FALSE, log.p = TRUE)
    # independent: log-sum of the log pmf beyond q
    lp <- call_f(f$d, (q + 1):(q + 4000), f$args, log = TRUE)
    lp <- lp[is.finite(lp)]
    ref <- max(lp) + log(sum(exp(lp - max(lp))))
    expect_true(is.finite(ls) && ls < log(1e-15), label = nm)
    expect_equal(ls, ref, tolerance = 1e-9, label = nm)
  }
})

test_that("quantiles are the smallest count reaching the probability", {
  pr <- c(1e-12, 0.01, 0.3, 0.5, 0.77, 0.99, 1 - 1e-12)
  for (nm in names(fams)) {
    f <- fams[[nm]]
    qq <- call_f(f$q, pr, f$args)
    F <- call_f(f$p, qq, f$args)
    Fm <- call_f(f$p, qq - 1, f$args)
    expect_true(all(F >= pr * (1 - 1e-10) & Fm < pr), label = nm)
    # equivalent specifications through the upper tail and log scale
    expect_equal(call_f(f$q, 1 - pr[2:6], f$args, lower.tail = FALSE),
                 qq[2:6], label = nm)
    expect_equal(call_f(f$q, log(pr), f$args, log.p = TRUE), qq,
                 label = nm)
    expect_equal(call_f(f$q, c(0, 1), f$args), c(0, Inf), label = nm)
    # very small upper tail probabilities
    qt <- call_f(f$q, -200, f$args, lower.tail = FALSE, log.p = TRUE)
    expect_true(call_f(f$p, qt, f$args, lower.tail = FALSE,
                       log.p = TRUE) <= -200 + 1e-9, label = nm)
    expect_true(call_f(f$p, qt - 1, f$args, lower.tail = FALSE,
                       log.p = TRUE) > -200, label = nm)
  }
})

test_that("random generation matches the distribution", {
  for (nm in names(fams)) {
    f <- fams[[nm]]
    set.seed(123)
    x1 <- call_f(f$r, 20000, f$args)
    set.seed(123)
    x2 <- call_f(f$r, 20000, f$args)
    expect_identical(x1, x2, label = nm)
    expect_true(all(x1 >= 0 & x1 == round(x1)), label = nm)
    y <- 0:2000
    pm <- call_f(f$d, y, f$args)
    m <- sum(y * pm)
    v <- sum(y^2 * pm) - m^2
    expect_lt(abs(mean(x1) - m), 4.5 * sqrt(v / 20000), label = nm)
    # observed frequencies of small counts (chi-square style bound)
    tab <- tabulate(x1 + 1, 6) / 20000
    se <- sqrt(pm[1:6] * (1 - pm[1:6]) / 20000)
    expect_true(all(abs(tab - pm[1:6]) < 4.5 * se + 1e-12), label = nm)
    expect_length(call_f(f$r, c(5, 1, 1), f$args), 3)
    expect_length(call_f(f$r, 0, f$args), 0)
  }
})

test_that("argument handling follows base R conventions", {
  for (nm in names(fams)) {
    f <- fams[[nm]]
    # recycling and zero-length arguments
    expect_length(call_f(f$d, 0:4, f$args), 5)
    expect_length(call_f(f$d, numeric(0), f$args), 0)
    expect_length(call_f(f$p, numeric(0), f$args), 0)
    expect_length(call_f(f$q, numeric(0), f$args), 0)
    # missing values
    expect_true(is.na(call_f(f$d, NA, f$args)), label = nm)
    expect_true(is.na(call_f(f$p, NA, f$args)), label = nm)
    expect_true(is.na(call_f(f$q, NA, f$args)), label = nm)
    na_args <- f$args
    na_args[[1]] <- NA
    expect_true(is.na(call_f(f$d, 1, na_args)), label = nm)
    # outside the support
    expect_equal(call_f(f$d, c(-1, Inf), f$args), c(0, 0), label = nm)
    expect_equal(call_f(f$p, c(-1, Inf), f$args), c(0, 1), label = nm)
    expect_warning(z <- call_f(f$d, 0.5, f$args), "non-integer",
                   label = nm)
    expect_equal(z, 0)
    # invalid parameters: NaN and a warning
    expect_warning(z <- call_f(f$d, 1, f$bad), "NaN", label = nm)
    expect_true(is.nan(z), label = nm)
    expect_warning(z <- call_f(f$p, 1, f$bad), "NaN", label = nm)
    expect_true(is.nan(z), label = nm)
    expect_warning(z <- call_f(f$q, 0.5, f$bad), "NaN", label = nm)
    expect_true(is.nan(z), label = nm)
    expect_warning(z <- call_f(f$q, 1.5, f$args), "NaN", label = nm)
    expect_true(is.nan(z), label = nm)
    expect_warning(z <- call_f(f$r, 2, f$bad), "NAs produced", label = nm)
    expect_true(all(is.na(z)), label = nm)
    expect_error(call_f(f$p, 1, f$args, lower.tail = NA), "lower.tail")
  }
})

test_that("vectorized parameters are evaluated elementwise", {
  x <- c(0, 3, 7)
  expect_equal(ddelaporte(x, mu = c(1, 4, 9), rho = c(0.2, 0.5, 0.9), k = 2),
               c(ddelaporte(0, 1, 0.2, 2), ddelaporte(3, 4, 0.5, 2),
                 ddelaporte(7, 9, 0.9, 2)))
  expect_equal(ddpo(x, lambda = c(1, 4, 9), phi = c(0.5, 1, 2)),
               c(ddpo(0, 1, 0.5), ddpo(3, 4, 1), ddpo(7, 9, 2)))
  expect_equal(qhpois(0.5, lambda = c(1, 4, 9), gamma = c(0.5, 1, 2)),
               c(qhpois(0.5, 1, 0.5), qhpois(0.5, 4, 1), qhpois(0.5, 9, 2)))
  expect_equal(dghermite(4, mu = 2, rho = 0.5, m = 2:4),
               c(dghermite(4, 2, 0.5, 2), dghermite(4, 2, 0.5, 3),
                 dghermite(4, 2, 0.5, 4)))
})

test_that("large means and extreme parameters stay finite", {
  expect_equal(sum(ddpo(0:300000, 1e5, 0.7)), 1, tolerance = 1e-10)
  lp <- dhpois(c(0, 1e5), 1e5, 50, log = TRUE)
  expect_true(all(is.finite(lp)))
  expect_true(is.finite(pdweib(1e6, 0.8, 1.3, lower.tail = FALSE,
                               log.p = TRUE)))
  expect_equal(pdweib(1e6, 0.8, 1.3, lower.tail = FALSE, log.p = TRUE),
               log(0.8) * (1e6 + 1)^1.3)
  expect_equal(sum(ddelaporte(0:5000, 1000, 0.5, 3)) +
                 pdelaporte(5000, 1000, 0.5, 3, lower.tail = FALSE), 1,
               tolerance = 1e-12)
})
