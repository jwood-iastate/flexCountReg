# Regression tests for the corrections from the development review of
# flexCountReg 0.1.3 (findings F01-F28).

test_that("F01: PGE and PW estimate only the shape parameter", {
  h <- randtoolbox::halton(200)
  # The scale of the mixing distribution cancels out of the mean-parameterized
  # PMF, so it cannot be estimated
  expect_equal(sapply(c(0.1, 1, 10, 100), function(s)
    dpge(3, mean = 2, shape = 1.5, scale = s, haltons = h)),
    rep(dpge(3, mean = 2, shape = 1.5, scale = 1, haltons = h), 4))
  expect_equal(sapply(c(0.1, 1, 10, 100), function(s)
    dpWeib_cpp(3, mean = 2, alpha = 1.5, sigma = s, h = h)),
    rep(dpWeib_cpp(3, mean = 2, alpha = 1.5, sigma = 1, h = h), 4))
  d <- washington_roads
  for (fam in c("PGE", "PW")) {
    m <- countreg(Total_crashes ~ lnaadt, d, family = fam, offset = "lnlength",
                  ndraws = 100, method = "BFGS")
    expect_length(coef(m), 3)
    expect_equal(attr(logLik(m), "df"), 3)
    expect_true(all(is.finite(sqrt(diag(vcov(m))))))
    expect_error(countreg(Total_crashes ~ lnaadt, d, family = fam,
                          dis_param_formula_2 = ~ speed50),
                 "cancels out of the likelihood")
    expect_error(countreg.rp(Total_crashes ~ lnaadt,
                             rpar_formula = ~ -1 + speed50, data = d,
                             family = fam, dis_param_formula_2 = ~ speed50),
                 "cancels out of the likelihood")
  }
})

test_that("F02: underdispersed GP1/GP2 use the moments of the implemented PMF", {
  y <- 0:300
  # The renormalized PMF has mean 0.803 and 0.948, not mu = 0.75
  m1 <- .fcr_family_mean("GP1", 0.75, 1 - 0.5)
  expect_equal(m1, sum(y * dgp1(y, 0.75, -0.5)), tolerance = 1e-12)
  expect_equal(m1, 0.80305, tolerance = 1e-5)
  expect_equal(.fcr_family_mean("GP2", 0.75, 1 - 0.9),
               sum(y * dgp2(y, 0.75, -0.9)), tolerance = 1e-12)
  v1 <- .fcr_family_variance("GP1", 0.75, 0.5)
  expect_equal(v1, sum(y^2 * dgp1(y, 0.75, -0.5)) - m1^2, tolerance = 1e-10)
  # All mass at zero: mean and variance are zero
  expect_equal(.fcr_family_mean("GP1", 0.2, 0.2), 0)
  # Nonnegative dispersion keeps the textbook moments
  expect_equal(.fcr_family_mean("GP1", 2, 1.3), 2)
  expect_equal(.fcr_family_variance("GP1", 2, 1.3), 2 * 1.3^2)

  # A fitted underdispersed model: predictions are the PMF means and the
  # Pearson residuals have variance close to 1
  set.seed(11)
  n <- 3000
  x <- runif(n)
  sim <- data.frame(x = x, y = rgp1(n, mu = exp(0.6 + 0.5 * x), phi = -0.3))
  m <- countreg(y ~ x, sim, family = "GP1", method = "BFGS")
  phi <- exp(coef(m)[["ln(phi+1)"]]) - 1
  expect_lt(phi, 0)
  mu <- exp(coef(m)[1] + coef(m)[2] * sim$x[1:5])
  pmf_mean <- vapply(mu, function(u) sum(y * dgp1(y, u, phi)), 0)
  expect_equal(predict(m, newdata = sim[1:5, ]), pmf_mean, tolerance = 1e-8)
  expect_equal(mean(predict(m)), mean(sim$y), tolerance = 0.01)
  expect_equal(var(residuals(m, type = "pearson")), 1, tolerance = 0.1)
})

test_that("F02: random-parameter GP predictions average the exact mean", {
  set.seed(12)
  n <- 300
  x <- runif(n)
  sim <- data.frame(x = x, y = rgp1(n, mu = exp(0.6 + 0.5 * x), phi = -0.3))
  # F09: the default start (dispersion 0) gives a finite likelihood
  rp <- countreg.rp(y ~ 1, rpar_formula = ~ -1 + x, data = sim,
                    family = "GP1", ndraws = 30, method = "BFGS",
                    max.iters = 100)
  expect_lt(exp(coef(rp)[["ln(phi+1)"]]) - 1, 0)
  sim_p <- predict(rp, newdata = sim[1:10, ], method = "Simulated")
  expect_message(ex <- predict(rp, newdata = sim[1:10, ], method = "Exact"),
                 "using Simulated")
  expect_equal(ex, sim_p)
  expect_equal(mean(predict(rp, newdata = sim)), mean(sim$y), tolerance = 0.05)
})
