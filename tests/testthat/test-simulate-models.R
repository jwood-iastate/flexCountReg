# simulate() for the random-parameter, random-effects and fixed-effects
# estimators (R/simulate_models.R)

sim_fits <- function() {
  data("washington_roads", package = "flexCountReg", envir = environment())
  w <- washington_roads
  f <- Total_crashes ~ lnaadt + offset(lnlength)
  list(
    data = w,
    repois = repois(f, "ID", w, method = "BFGS"),
    renb = renb(f, "ID", w, method = "BFGS"),
    poislind_re = poislind_re(f, "ID", w, method = "BFGS"),
    negmultinom = negmultinom(f, group_var = "ID", data = w, family = "NM1"),
    countreg_rp = countreg_rp(f, rpar_formula = ~ -1 + speed50, data = w,
                              family = "NB2", ndraws = 100, method = "BFGS",
                              group_var = "ID"),
    countreg_fe = suppressMessages(countreg_fe(Total_crashes ~ lnaadt,
                                               fixef = ~ ID, data = w)))
}

test_that("every estimator simulates its training rows reproducibly", {
  skip_on_cran()
  fits <- sim_fits()
  for (nm in setdiff(names(fits), "data")) {
    m <- fits[[nm]]
    s <- simulate(m, nsim = 3, seed = 1)
    expect_s3_class(s, "data.frame")
    expect_equal(dim(s), c(nobs(m), 3), info = nm)
    expect_true(all(s >= 0 & s == round(s)), info = nm)
    expect_identical(s, simulate(m, nsim = 3, seed = 1), info = nm)
    # Simulated means are the model's expected counts
    S <- as.matrix(simulate(m, nsim = 300, seed = 2))
    pred <- if (nm == "countreg_rp") predict(m, method = "Exact") else
      suppressMessages(predict(m, newdata = flexCountReg:::.fcr_training_rows(m)))
    expect_equal(mean(S), mean(pred), tolerance = 0.03, info = nm)
  }
})

test_that("panel models draw one random effect per panel", {
  skip_on_cran()
  fits <- sim_fits()
  m <- fits$repois
  d <- flexCountReg:::.fcr_training_rows(m)
  S <- as.matrix(simulate(m, nsim = 400, seed = 3))
  # Rows of a panel share their effect: the panel totals vary more than
  # the sum of the rows' variances
  tot <- rowsum(S, d$ID)
  expect_gt(mean(apply(tot, 1, stats::var)),
            1.1 * mean(rowsum(apply(S, 1, stats::var), d$ID)))
  # Variance of a row: mu + alpha mu^2
  a <- exp(coef(m)[["lnalpha"]])
  mu <- predict(m)
  i <- which.max(mu)
  expect_equal(stats::var(S[i, ]), mu[i] + a * mu[i]^2, tolerance = 0.25)
})

test_that("a refit on simulated data recovers the parameters", {
  skip_on_cran()
  fits <- sim_fits()
  m <- fits$repois
  d <- flexCountReg:::.fcr_training_rows(m)
  d$Total_crashes <- simulate(m, nsim = 1, seed = 4)$sim_1
  r <- repois(Total_crashes ~ lnaadt + offset(lnlength), "ID", d,
              method = "BFGS")
  z <- (coef(r) - coef(m)) / sqrt(diag(vcov(m)))
  expect_true(all(abs(z) < 4))
})
