# Property tests for the count distributions' d/p/q/r functions: the PMF sums
# to one, the CDF is the cumulative PMF (also in the upper tail and on the
# log scale), the quantile function inverts the CDF, log densities agree,
# and random draws reproduce the mean.

K <- 0:150

check_count_family <- function(name, d, p, q, r, n_r = 20000, tol_r = 0.05,
                               p_tail = TRUE, p_log = TRUE, d_log = TRUE,
                               q_tail = FALSE) {
  pmf <- d(K)
  expect_equal(sum(pmf), 1, tolerance = 1e-4, label = paste(name, "PMF sum"))
  expect_true(all(pmf >= 0), label = paste(name, "PMF >= 0"))
  ks <- c(0, 1, 2, 5)
  cdf <- p(ks)
  expect_equal(cdf, cumsum(pmf)[ks + 1], tolerance = 1e-6,
               label = paste(name, "CDF"))
  if (p_tail) {
    expect_equal(p(ks, lower.tail = FALSE), 1 - cdf, tolerance = 1e-6,
                 label = paste(name, "upper tail"))
  }
  if (p_log) {
    expect_equal(p(ks, log.p = TRUE), log(cdf), tolerance = 1e-6,
                 label = paste(name, "log CDF"))
  }
  if (d_log) {
    expect_equal(d(ks, log = TRUE), log(pmf[ks + 1]), tolerance = 1e-8,
                 label = paste(name, "log PMF"))
  }
  # The quantile at each CDF value is that count (largest-to-smallest checks)
  probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  qs <- q(probs)
  cdf_all <- cumsum(pmf)
  expect_true(all(cdf_all[qs + 1] >= probs - 1e-8),
              label = paste(name, "q: F(q) >= p"))
  below <- qs > 0
  expect_true(all(cdf_all[qs[below]] < probs[below] + 1e-8),
              label = paste(name, "q: F(q - 1) < p"))
  if (q_tail) {
    expect_equal(q(1 - probs, lower.tail = FALSE), qs,
                 label = paste(name, "q upper tail"))
    expect_equal(q(log(probs), log.p = TRUE), qs,
                 label = paste(name, "q log.p"))
  }
  set.seed(1)
  draws <- r(n_r)
  expect_length(draws, n_r)
  expect_equal(mean(draws), sum(K * pmf), tolerance = tol_r,
               label = paste(name, "r mean"))
}

test_that("Generalized Poisson 1 and 2 (over- and underdispersed)", {
  for (phi in c(0.4, -0.2)) {
    check_count_family(
      paste("GP1", phi),
      d = function(x, ...) dgp1(x, mu = 2, phi = phi, ...),
      p = function(q, ...) pgp1(q, mu = 2, phi = phi, ...),
      q = function(pr, ...) qgp1(pr, mu = 2, phi = phi, ...),
      r = function(n) rgp1(n, mu = 2, phi = phi), q_tail = TRUE)
  }
  for (alpha in c(0.3, -0.1)) {
    check_count_family(
      paste("GP2", alpha),
      d = function(x, ...) dgp2(x, mu = 2, alpha = alpha, ...),
      p = function(q, ...) pgp2(q, mu = 2, alpha = alpha, ...),
      q = function(pr, ...) qgp2(pr, mu = 2, alpha = alpha, ...),
      r = function(n) rgp2(n, mu = 2, alpha = alpha), q_tail = TRUE)
  }
  # phi = 0 and alpha = 0 are the Poisson
  expect_equal(dgp1(0:10, 2, 0), dpois(0:10, 2))
  expect_equal(pgp1(0:10, 2, 0), ppois(0:10, 2))
  expect_equal(qgp1(c(0.2, 0.7), 2, 0), qpois(c(0.2, 0.7), 2))
  expect_equal(dgp2(0:10, 2, 0), dpois(0:10, 2))
  expect_equal(pgp2(0:10, 2, 0), ppois(0:10, 2))
  # Edge probabilities
  expect_equal(qgp1(0, 2, 0.4), 0)
  expect_equal(qgp1(1, 2, 0.4), Inf)
  expect_true(is.nan(qgp1(1.5, 2, 0.4)))
  expect_equal(pgp1(-1, 2, 0.4), 0)
  expect_equal(pgp2(-1, 2, 0.3), 0)
  expect_warning(v <- dgp2(1, mu = 2, alpha = -0.6), "Invalid")  # 1 + a mu < 0
  expect_true(is.na(v))
  expect_warning(dgp1(1, mu = -1, phi = 0.2), "mu")
  expect_warning(dgp1(1, mu = 1, phi = -2), "phi")
})

test_that("Conway-Maxwell-Poisson by mean and by rate", {
  check_count_family(
    "COM mu",
    d = function(x, ...) dcom(x, mu = 2, nu = 0.7, ...),
    p = function(q, ...) pcom(q, mu = 2, nu = 0.7, ...),
    q = function(pr) qcom(pr, mu = 2, nu = 0.7),
    r = function(n) rcom(n, mu = 2, nu = 0.7))
  check_count_family(
    "COM lambda",
    d = function(x, ...) dcom(x, lambda = 1.5, nu = 1.3, ...),
    p = function(q, ...) pcom(q, lambda = 1.5, nu = 1.3, ...),
    q = function(pr) qcom(pr, lambda = 1.5, nu = 1.3),
    r = function(n) rcom(n, lambda = 1.5, nu = 1.3))
  # nu = 1 is the Poisson with rate lambda
  expect_equal(dcom(0:8, lambda = 2, nu = 1), dpois(0:8, 2), tolerance = 1e-8)
  expect_equal(sum((0:200) * dcom(0:200, mu = 3, nu = 0.5)), 3,
               tolerance = 1e-6)
})

test_that("Poisson-Weibull by lambda and by the Weibull mean and SD", {
  check_count_family(
    "PW",
    d = function(x, ...) dpoisweibull(x, lambda = 1.5, alpha = 2, sigma = 1,
                                      ...),
    p = function(q, ...) ppoisweibull(q, lambda = 1.5, alpha = 2, sigma = 1,
                                      ...),
    q = function(pr) qpoisweibull(pr, lambda = 1.5, alpha = 2, sigma = 1),
    r = function(n) rpoisweibull(n, lambda = 1.5, alpha = 2, sigma = 1),
    tol_r = 0.02)
  # Mean of the mixture is lambda * sigma * Gamma(1 + 1/alpha)
  pmf <- dpoisweibull(K, lambda = 1.5, alpha = 2, sigma = 1, ndraws = 5000)
  expect_equal(sum(K * pmf), 1.5 * gamma(1.5), tolerance = 5e-3)
  # The Weibull mean and SD determine alpha and sigma
  pars <- calculate_params(mean_value = 2, sd_value = 1)
  expect_equal(pars[2] * gamma(1 + 1 / pars[1]), 2, tolerance = 1e-6)
  expect_equal(pars[2] * sqrt(gamma(1 + 2 / pars[1]) -
                                gamma(1 + 1 / pars[1])^2), 1, tolerance = 1e-3)
  a <- dpoisweibull(0:5, mean_value = 2, sd_value = 1)
  b <- dpoisweibull(0:5, lambda = 1, alpha = pars[1], sigma = pars[2])
  expect_equal(a, b)
  expect_equal(ppoisweibull(3, mean_value = 2, sd_value = 1), sum(a[1:4]))
  expect_equal(qpoisweibull(0.5, mean_value = 2, sd_value = 1),
               qpoisweibull(0.5, lambda = 1, alpha = pars[1], sigma = pars[2]))
  set.seed(2)
  expect_equal(mean(rpoisweibull(20000, mean_value = 2, sd_value = 1)), 2,
               tolerance = 0.03)
  expect_equal(ppoisweibull(-1, lambda = 1, alpha = 2, sigma = 1), 0)
  expect_equal(ppoisweibull(-1, lambda = 1, alpha = 2, sigma = 1,
                            log.p = TRUE), -Inf)
  expect_warning(calculate_params(), "Insufficient parameters")
  expect_error(qpoisweibull(0.5, lambda = 1), "alpha and sigma")
  expect_error(qpoisweibull(0.5, alpha = 2, sigma = 1), "lambda or mean")
})

test_that("Poisson-generalized-exponential", {
  h <- randtoolbox::halton(2000)
  check_count_family(
    "PGE",
    d = function(x, ...) dpge(x, mean = 2, shape = 2, scale = 1,
                              haltons = h, ...),
    p = function(q, ...) ppge(q, mean = 2, shape = 2, scale = 1,
                              haltons = h, ...),
    q = function(pr) qpge(pr, mean = 2, shape = 2, scale = 1, ndraws = 2000),
    r = function(n) rpge(n, mean = 2, shape = 2, scale = 1))
  expect_error(dpge(1, mean = -1), "mean")
  expect_error(dpge(1, shape = 0), "shape")
  expect_error(dpge(1, scale = -1), "scale")
  expect_error(dpge(1, haltons = c(0, 0.5)), "haltons")
  expect_error(dpge(1, ndraws = 0), "ndraws")
  expect_equal(dpge(numeric(0)), numeric(0))
  expect_warning(v <- ppge(1, mean = -1), "mean")
  expect_true(is.nan(v))
  expect_equal(ppge(-1, haltons = h), 0)
})

test_that("Poisson-lognormal", {
  h <- randtoolbox::halton(2000)
  check_count_family(
    "PLN",
    d = function(x, ...) dpLnorm(x, mean = 2, sigma = 0.6, hdraws = h, ...),
    p = function(q, ...) ppLnorm(q, mean = 2, sigma = 0.6, ndraws = 2000,
                                 ...),
    q = function(pr) qpLnorm(pr, mean = 2, sigma = 0.6, ndraws = 2000),
    r = function(n) rpLnorm(n, mean = 2, sigma = 0.6), tol_r = 0.03)
  expect_error(dpLnorm(1, mean = -1), "mean")
  expect_error(dpLnorm(1, sigma = 0), "sigma")
  expect_equal(dpLnorm(numeric(0)), numeric(0))
  expect_warning(ppLnorm(1, mean = -1), "greater than 0")
  expect_equal(ppLnorm(-1), 0)
  skip_if_not_installed("poilog")
  expect_equal(dpLnorm(0:3, mean = 2, sigma = 0.6, engine = "poilog"),
               dpLnorm(0:3, mean = 2, sigma = 0.6, ndraws = 20000),
               tolerance = 1e-3)
})

test_that("Poisson-inverse-gamma and Poisson-inverse-Gaussian", {
  check_count_family(
    "PInvGamma",
    d = function(x, ...) dpinvgamma(x, mu = 2, eta = 0.4, ...),
    p = function(q, ...) ppinvgamma(q, mu = 2, eta = 0.4, ...),
    q = function(pr) qpinvgamma(pr, mu = 2, eta = 0.4),
    r = function(n) rpinvgamma(n, mu = 2, eta = 0.4), tol_r = 0.05)
  for (form in c("Type 1", "Type 2")) {
    check_count_family(
      paste("PIG", form),
      d = function(x, ...) dpinvgaus(x, mu = 2, eta = 0.5, form = form, ...),
      p = function(q, ...) ppinvgaus(q, mu = 2, eta = 0.5, form = form, ...),
      q = function(pr) qpinvgaus(pr, mu = 2, eta = 0.5, form = form),
      r = function(n) rpinvgaus(n, mu = 2, eta = 0.5, form = form))
  }
  expect_equal(ppinvgamma(-1, 2, 0.4), 0)
  expect_warning(v <- ppinvgamma(c(1, 2), mu = c(-1, 2), eta = 1), "mu")
  expect_true(is.nan(v[1]))
  expect_equal(v[2], ppinvgamma(2, 2, 1))
  expect_warning(ppinvgamma(1, mu = 1, eta = -1), "eta")
  expect_error(rpinvgamma(-1, 1, 1), "non-negative")
  expect_equal(rpinvgamma(0, 1, 1), numeric(0))
  expect_error(rpinvgamma(2, mu = -1), "mu")
  expect_error(rpinvgamma(2, eta = -1), "eta")
  expect_error(rpinvgaus(2, mu = -1), "mu")
  expect_error(rpinvgaus(2, eta = -1), "eta")
  expect_error(rpinvgaus(c(-1)), "non-negative")
  expect_warning(v <- ppinvgaus(1, mu = -1), "mu")
  expect_true(is.nan(v))
  expect_warning(ppinvgaus(1, eta = -1), "eta")
})

test_that("Poisson-Lindley family", {
  check_count_family(
    "PL",
    d = function(x, ...) dplind(x, mean = 2, theta = 1.5, ...),
    p = function(q, ...) pplind(q, mean = 2, theta = 1.5, ...),
    q = function(pr) qplind(pr, mean = 2, theta = 1.5),
    r = function(n) rplind(n, mean = 2, theta = 1.5))
  h <- randtoolbox::halton(2000)
  check_count_family(
    "PLL",
    d = function(x, ...) dplindLnorm(x, mean = 2, theta = 1.5, sigma = 0.4,
                                     hdraws = h, ...),
    p = function(q, ...) pplindLnorm(q, mean = 2, theta = 1.5, sigma = 0.4,
                                     ndraws = 2000, ...),
    q = function(pr) qplindLnorm(pr, mean = 2, theta = 1.5, sigma = 0.4,
                                 ndraws = 2000),
    r = function(n) rplindLnorm(n, mean = 2, theta = 1.5, sigma = 0.4),
    tol_r = 0.05)
  check_count_family(
    "PLG",
    d = function(x, ...) dplindGamma(x, mean = 2, theta = 1.5, alpha = 0.4,
                                     ...),
    p = function(q, ...) pplindGamma(q, mean = 2, theta = 1.5, alpha = 0.4,
                                     ...),
    q = function(pr) qplindGamma(pr, mean = 2, theta = 1.5, alpha = 0.4),
    r = function(n) rplindGamma(n, mean = 2, theta = 1.5, alpha = 0.4))
  # lambda parameterizations agree with the mean parameterization
  lam <- 2 * 1.5 * 2.5 / 3.5
  expect_equal(pplind(0:4, lambda = lam, theta = 1.5),
               pplind(0:4, mean = 2, theta = 1.5))
  expect_equal(qplindLnorm(0.5, lambda = lam, theta = 1.5, sigma = 0.4),
               qplindLnorm(0.5, mean = 2, theta = 1.5, sigma = 0.4))
  expect_equal(rplind(3, lambda = c(1, 2, 3), theta = 1)[0], numeric(0))
  # alpha = 0 is the Poisson-Lindley limit; mean = 0 is degenerate at zero
  expect_equal(dplindGamma(0:5, mean = 2, theta = 1.5, alpha = 0),
               dplind(0:5, mean = 2, theta = 1.5), tolerance = 1e-6)
  expect_equal(dplindGamma(0:2, mean = 0, theta = 1.5, alpha = 0.4),
               c(1, 0, 0))
  expect_equal(qplindGamma(c(0, 1), 2, 1.5, 0.4), c(0, Inf))
  expect_warning(expect_true(is.nan(qplindGamma(2, 2, 1.5, 0.4))))
  expect_true(is.na(qplindGamma(NA_real_, 2, 1.5, 0.4)))
  expect_error(qplindGamma(NA, 2, 1.5, 0.4), "numeric")
  expect_error(rplindGamma(-1), "non-negative")
  expect_warning(v <- pplind(1, mean = -1), "positive")
  expect_true(is.nan(v))
  expect_warning(pplind(1, lambda = -1), "positive")
  expect_equal(pplind(-1, 2, 1.5), 0)
  expect_warning(v <- qplind(0.5, mean = -1), "greater than 0")
  expect_true(is.nan(v))
  expect_warning(v <- qplindLnorm(0.5, mean = -1), "greater than 0")
  expect_true(is.nan(v))
})

test_that("Sichel and generalized Waring", {
  check_count_family(
    "Sichel",
    d = function(x, ...) dsichel(x, mu = 2, sigma = 1, gamma = -0.5, ...),
    p = function(q, ...) psichel(q, mu = 2, sigma = 1, gamma = -0.5, ...),
    q = function(pr, ...) qsichel(pr, mu = 2, sigma = 1, gamma = -0.5, ...),
    r = function(n) rsichel(n, mu = 2, sigma = 1, gamma = -0.5),
    q_tail = TRUE)
  check_count_family(
    "GW",
    d = function(x, ...) dgwar(x, mu = 2, k = 2, rho = 4, ...),
    p = function(q, ...) pgwar(q, mu = 2, k = 2, rho = 4, ...),
    q = function(pr) qgwar(pr, mu = 2, k = 2, rho = 4),
    r = function(n) rgwar(n, mu = 2, k = 2, rho = 4), tol_r = 0.06)
  expect_error(dgwar(1, mu = -1, k = 1, rho = 3), "mu")
  expect_error(dgwar(1, mu = 1, k = -1, rho = 3), "k")
  expect_error(dgwar(1, mu = 1, k = 1, rho = 0.5), "rho")
  expect_error(pgwar(1, mu = -1, k = 1, rho = 3), "mu")
  expect_error(qgwar(2, mu = 1, k = 1, rho = 3), "p must")
  expect_equal(dgwar(numeric(0), 1, 1, 3), numeric(0))
})

test_that("continuous helper distributions", {
  # Lindley
  x <- seq(0.1, 5, by = 0.5)
  expect_equal(integrate(dlindley, 0, Inf, theta = 1.5)$value, 1,
               tolerance = 1e-6)
  expect_equal(plindley(x, 1.5),
               vapply(x, function(u) integrate(dlindley, 0, u,
                                               theta = 1.5)$value, 0),
               tolerance = 1e-6)
  expect_equal(plindley(x, 1.5, lower.tail = FALSE), 1 - plindley(x, 1.5))
  expect_equal(dlindley(x, 1.5, log = TRUE), log(dlindley(x, 1.5)))
  skip_if_not_installed("lamW")
  pr <- c(0.1, 0.5, 0.9)
  expect_equal(plindley(qlindley(pr, 1.5), 1.5), pr, tolerance = 1e-8)
  expect_equal(qlindley(pr, 1.5, lower.tail = FALSE), qlindley(1 - pr, 1.5))
  expect_equal(qlindley(log(pr), 1.5, log.p = TRUE), qlindley(pr, 1.5))
})

test_that("inverse gamma and inverse Gaussian helpers", {
  x <- c(0.5, 1, 2)
  expect_equal(integrate(dinvgamma, 0, Inf, shape = 3, scale = 2)$value, 1,
               tolerance = 1e-6)
  expect_equal(pinvgamma(x, 3, 2), 1 - pgamma(1 / x, shape = 3, rate = 2),
               tolerance = 1e-10)
  expect_equal(pinvgamma(x, 3, 2, lower.tail = FALSE), 1 - pinvgamma(x, 3, 2))
  expect_equal(pinvgamma(x, 3, 2, log.p = TRUE), log(pinvgamma(x, 3, 2)))
  expect_equal(dinvgamma(x, 3, 2, log = TRUE), log(dinvgamma(x, 3, 2)))
  pr <- c(0.2, 0.5, 0.8)
  expect_equal(pinvgamma(qinvgamma(pr, 3, 2), 3, 2), pr, tolerance = 1e-8)
  expect_equal(qinvgamma(pr, 3, 2, lower.tail = FALSE), qinvgamma(1 - pr, 3, 2))
  expect_equal(qinvgamma(log(pr), 3, 2, log.p = TRUE), qinvgamma(pr, 3, 2))
  set.seed(3)
  expect_equal(mean(rinvgamma(50000, 3, 2)), 2 / (3 - 1), tolerance = 0.02)
  set.seed(4)
  w <- rinvgaus(50000, mean = 2, shape = 3)
  expect_equal(mean(w), 2, tolerance = 0.02)
  expect_equal(var(w), 2^3 / 3, tolerance = 0.05)
})
