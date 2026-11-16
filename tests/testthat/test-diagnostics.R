# Probability predictions, quantile residuals, rootograms, and a check of
# every registered family's density, mean and variance.

sim_nb2 <- function(n = 3000, seed = 11) {
  set.seed(seed)
  x <- rnorm(n)
  data.frame(x = x, y = rnbinom(n, size = 1 / 0.5, mu = exp(0.3 + 0.4 * x)))
}

test_that("predict(type = 'prob') gives P(Y = k) per row", {
  m <- countreg(Total_crashes ~ lnaadt + lnlength, washington_roads,
                method = "BFGS")
  P <- predict(m, newdata = washington_roads[1:10, ], type = "prob",
               max_count = 60)
  expect_equal(dim(P), c(10, 61))
  expect_equal(rowSums(P), rep(1, 10), tolerance = 1e-6)
  mu <- predict(m, newdata = washington_roads[1:10, ])
  alpha <- exp(coef(m)["ln(alpha)"])
  expect_equal(P[, "2"], dnbinom(2, size = 1 / alpha, mu = mu))
  nd <- washington_roads[1:3, ]
  nd$lnaadt[2] <- NA
  P2 <- predict(m, newdata = nd, type = "prob", max_count = 5)
  expect_true(all(is.na(P2[2, ])))
})

test_that("quantile residuals are standard normal for a correct model", {
  d <- sim_nb2()
  m <- countreg(y ~ x, d, family = "NB2", method = "BFGS")
  set.seed(2)
  r <- residuals(m, type = "quantile")
  expect_length(r, nrow(d))
  expect_equal(mean(r), 0, tolerance = 0.06)
  expect_equal(sd(r), 1, tolerance = 0.06)
  expect_gt(stats::ks.test(r, "pnorm")$p.value, 0.01)
  # A Poisson fit to the same overdispersed data has too-wide residuals
  mp <- countreg(y ~ x, d, family = "Poisson", method = "BFGS")
  expect_gt(sd(residuals(mp, type = "quantile")), 1.1)
})

test_that("rootogram compares observed and expected frequencies", {
  d <- sim_nb2()
  m <- countreg(y ~ x, d, family = "NB2", method = "BFGS")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  r <- rootogram(m)
  expect_s3_class(r, "flexCountReg_rootogram")
  expect_equal(sum(r$observed), nrow(d))
  expect_equal(r$observed, tabulate(d$y + 1, nbins = max(d$y) + 1))
  # Expected frequencies are close to the observed for the true model
  expect_lt(max(abs(r$expected - r$observed)[1:5] / r$observed[1:5]), 0.1)
  expect_no_error(rootogram(m, style = "standing", max_count = 8))
})

test_that("every family's density sums to 1 with the stated mean and variance", {
  h <- randtoolbox::halton(5000)
  pars <- list(POISSON = c(NA, NA), GP1 = c(1.3, NA), GP2 = c(1.2, NA),
               NB1 = c(0.5, NA), NB2 = c(0.5, NA), NBP = c(0.5, 1.5),
               PLN = c(0.5, NA), PGE = c(2, NA), PIG1 = c(0.5, NA),
               PIG2 = c(0.5, NA), PIG = c(0.3, NA), PL = c(1.5, NA),
               PLG = c(1.5, 0.4), PLL = c(1.5, 0.5), PW = c(2, NA),
               SI = c(exp(0.5), 1), GW = c(2, 3), COM = c(0.7, NA))
  expect_setequal(names(pars), .fcr_families)
  y <- 0:2000
  mu <- 2
  for (fam in names(pars)) {
    a <- if (is.na(pars[[fam]][1])) NULL else pars[[fam]][1]
    s <- if (is.na(pars[[fam]][2])) NULL else pars[[fam]][2]
    p <- get_probFunc(fam)(y = y, predicted = rep(mu, length(y)), alpha = a,
                           sigma = s, haltons = h,
                           normed_haltons = stats::qnorm(h),
                           engine = "halton")
    # Simulation-based densities (Halton integration) are less exact
    tol <- if (fam %in% c("PLN", "PGE", "PLL", "PW")) 0.03 else 1e-4
    expect_equal(sum(p), 1, tolerance = tol, label = paste(fam, "sum"))
    m1 <- sum(y * p) / sum(p)
    expect_equal(m1, mu * .fcr_mean_adjust(fam, a, s), tolerance = tol,
                 label = paste(fam, "mean"))
    v <- sum(y^2 * p) / sum(p) - m1^2
    expect_equal(v, .fcr_family_variance(fam, mu, a, s), tolerance = tol,
                 label = paste(fam, "variance"))
  }
})

test_that("simulate() matches the fitted mean with underreporting", {
  m <- countreg(Total_crashes ~ lnaadt + lnlength, washington_roads,
                underreport_formula = ~ ShouldWidth04, method = "BFGS")
  sims <- simulate(m, nsim = 200, seed = 3)
  expect_equal(mean(as.matrix(sims)), mean(fitted(m)), tolerance = 0.03)
})
