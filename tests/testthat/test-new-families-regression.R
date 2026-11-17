# The count families added to countreg() and countreg.rp() in 0.1.2. All are
# parameterized by the mean mu = exp(X beta); the registry passes the
# distribution parameters as exp(W delta).

new_fams <- c("DW", "DPO", "DELAPORTE", "HERMITE", "HERMITE3", "HERMITE4",
              "HYPERPOIS", "GAMMACOUNT", "PA", "NBGE")

# Distribution parameters (natural scale, as the registry receives them)
fam_pars <- list(
  DW = list(1.3, NULL), DPO = list(2.5, NULL), DELAPORTE = list(1.5, 2),
  HERMITE = list(0.4, NULL), HERMITE3 = list(0.4, NULL),
  HERMITE4 = list(0.4, NULL), HYPERPOIS = list(0.4, NULL),
  GAMMACOUNT = list(2, NULL), PA = list(0.6, NULL), NBGE = list(1.5, 0.4))

reg_dens <- function(fam, y, mu, log = FALSE) {
  pr <- fam_pars[[fam]]
  n <- length(y)
  flexCountReg:::get_probFunc(fam)(
    y = y, predicted = rep(mu, n), alpha = rep(pr[[1]], n),
    sigma = if (!is.null(pr[[2]])) rep(pr[[2]], n), log = log)
}

test_that("every new family is registered with its parameters", {
  fams <- flexCountReg:::.fcr_families
  expect_true(all(new_fams %in% fams))
  params <- lapply(new_fams, flexCountReg:::get_params)
  names(params) <- new_fams
  expect_identical(params$DELAPORTE, list("logit(rho)", "ln(k)"))
  expect_identical(params$NBGE, list("ln(shape)", "ln(alpha)"))
  expect_identical(params$DPO, list("ln(phi)", NULL))
})

test_that("the regression families are parameterized by the mean", {
  y <- 0:3000
  for (fam in new_fams) {
    for (mu in c(0.3, 4, 25)) {
      pm <- reg_dens(fam, y, mu)
      m <- sum(y * pm)
      v <- sum(y^2 * pm) - m^2
      pr <- fam_pars[[fam]]
      expect_equal(sum(pm), 1, tolerance = 1e-10, label = paste(fam, mu))
      expect_equal(m, mu, tolerance = 1e-9, label = paste(fam, mu))
      expect_equal(
        flexCountReg:::.fcr_family_variance(fam, mu, pr[[1]], pr[[2]]),
        v, tolerance = 1e-8, label = paste(fam, mu))
      # log scale is consistent
      expect_equal(reg_dens(fam, 0:20, mu, log = TRUE), log(pm[1:21]),
                   tolerance = 1e-12, label = paste(fam, mu))
    }
  }
})

test_that("regression densities match the public distribution functions", {
  y <- 0:40
  mu <- 4
  rho <- stats::plogis(log(1.5))
  expect_equal(reg_dens("DELAPORTE", y, mu),
               ddelaporte(y, mu = mu, rho = rho, k = 2))
  expect_equal(reg_dens("HERMITE3", y, mu),
               dghermite(y, mu = mu, rho = stats::plogis(log(0.4)), m = 3))
  pa <- flexCountReg:::.pt15_to_paeppli(mu, 0.6)
  expect_equal(reg_dens("PA", y, mu), dpaeppli(y, pa$lambda, pa$prob))
  expect_equal(reg_dens("NBGE", y, mu), dnbge(y, mu, 1.5, 0.4))
  # Mean-solved families: the natural parameter found from mu
  lam <- exp(flexCountReg:::tab_mp_cpp(0, mu, 2.5, 0L, 4L))
  expect_equal(reg_dens("DPO", y, mu), ddpo(y, lam, 2.5), tolerance = 1e-12)
  expect_equal(flexCountReg:::.dpo_moments(lam, 2.5)$mean, mu,
               tolerance = 1e-12)
  lam <- exp(flexCountReg:::tab_mp_cpp(0, mu, 0.4, 1L, 4L))
  expect_equal(reg_dens("HYPERPOIS", y, mu), dhpois(y, lam, 0.4),
               tolerance = 1e-12)
  cc <- exp(flexCountReg:::dw_mp_cpp(mu, 1.3, 0L))
  expect_equal(reg_dens("DW", y, mu), ddweib(y, exp(-cc), 1.3),
               tolerance = 1e-12)
  s <- exp(flexCountReg:::gc_mp_cpp(mu, 2, 0L))
  expect_equal(reg_dens("GAMMACOUNT", y, mu), dgcount(y, 2, s),
               tolerance = 1e-12)
})

test_that("the mean equations are solved across the parameter range", {
  mu <- c(1e-4, 0.05, 1, 30, 2000)
  # DW: sum exp(-c k^b) directly from c (theta = exp(-c) can be too close
  # to 1 to recover c from it); shape 0.3 has heavy tails, so the direct
  # check uses the moderate means only
  for (b in c(0.3, 1, 4)) {
    cc <- exp(flexCountReg:::dw_mp_cpp(mu, rep(b, 5), 0L))
    use <- if (b < 1) 1:3 else 1:5
    m <- vapply(use, function(i) sum(exp(-cc[i] * (1:2e6)^b)), 0)
    expect_equal(m, mu[use], tolerance = 1e-10, label = paste("DW", b))
  }
  for (par in c(0.2, 1, 5)) {
    lam <- exp(flexCountReg:::tab_mp_cpp(rep(0, 5), mu, rep(par, 5), 0L, 4L))
    expect_equal(flexCountReg:::.dpo_moments(lam, par)$mean, mu,
                 tolerance = 1e-10, label = paste("DPO", par))
    lam <- exp(flexCountReg:::tab_mp_cpp(rep(0, 5), mu, rep(par, 5), 1L, 4L))
    expect_equal(flexCountReg:::.hpois_moments(lam, par)$mean, mu,
                 tolerance = 1e-10, label = paste("HP", par))
    s <- exp(flexCountReg:::gc_mp_cpp(mu, rep(par, 5), 0L))
    expect_equal(flexCountReg:::.gcount_moments(par, s, 1)$mean, mu,
                 tolerance = 1e-10, label = paste("GC", par))
  }
  # Outside the numerical domain: density 0 (not an error or a hang)
  f <- flexCountReg:::get_probFunc("DPO")
  expect_equal(f(y = 1, predicted = 3, alpha = 1e7, sigma = NULL,
                 log = TRUE), -Inf)
  expect_equal(f(y = 1, predicted = 1e12, alpha = 2, sigma = NULL,
                 log = TRUE), -Inf)
})

test_that("cdf hooks match sums of the regression densities", {
  for (fam in new_fams) {
    pr <- fam_pars[[fam]]
    cdf <- flexCountReg:::.fcr_family_registry[[fam]]$cdf
    q <- c(0, 2, 7, 30)
    pm <- reg_dens(fam, 0:3000, 4)
    lo <- cdf(q, rep(4, 4), rep(pr[[1]], 4),
              if (!is.null(pr[[2]])) rep(pr[[2]], 4), TRUE)
    up <- cdf(q, rep(4, 4), rep(pr[[1]], 4),
              if (!is.null(pr[[2]])) rep(pr[[2]], 4), FALSE)
    expect_equal(exp(lo), cumsum(pm)[q + 1], tolerance = 1e-10, label = fam)
    expect_equal(exp(up), 1 - cumsum(pm)[q + 1], tolerance = 1e-8,
                 label = fam)
  }
})

test_that("binomial thinning maps the parameters as documented", {
  p <- 0.6
  ytrue <- 0:400
  k <- 0:40
  thin_pmf <- function(pm) {
    vapply(k, function(v) sum(pm * stats::dbinom(v, ytrue, p)), 0)
  }
  for (fam in c("DELAPORTE", "NBGE", "HERMITE", "PA")) {
    pr <- fam_pars[[fam]]
    th <- flexCountReg:::.fcr_thin_params(fam, pr[[1]], pr[[2]], p)
    n <- length(k)
    reported <- flexCountReg:::get_probFunc(fam)(
      y = k, predicted = rep(p * 3, n), alpha = rep(th$alpha, n),
      sigma = if (!is.null(th$sigma)) rep(th$sigma, n))
    expect_equal(reported, thin_pmf(reg_dens(fam, ytrue, 3)),
                 tolerance = 1e-9, label = fam)
  }
  expect_error(flexCountReg:::.fcr_check_thinning("DW"), "not available")
  expect_error(flexCountReg:::.fcr_check_thinning("HERMITE3"),
               "not available")
})

test_that("countreg() recovers the parameters of simulated data", {
  skip_on_cran()
  set.seed(2026)
  n <- 400
  d <- data.frame(x = stats::rnorm(n))
  mu <- exp(0.8 + 0.5 * d$x)
  for (fam in new_fams) {
    pr <- fam_pars[[fam]]
    d$y <- flexCountReg:::.fcr_family_r(fam, n, mu, pr[[1]], pr[[2]])
    fit <- countreg(y ~ x, data = d, family = fam, method = "BFGS")
    expect_true(isTRUE(fit$model$converged), label = fam)
    b <- stats::coef(fit)[1:2]
    se <- sqrt(diag(stats::vcov(fit)))[1:2]
    expect_true(all(abs(b - c(0.8, 0.5)) < 4 * se), label = fam)
    # predictions are the mean
    expect_equal(as.numeric(stats::predict(fit, type = "response")),
                 as.numeric(exp(b[1] + b[2] * d$x)), tolerance = 1e-8,
                 label = fam)
  }
})

test_that("Poisson is nested in DPO, HYPERPOIS and GAMMACOUNT", {
  skip_on_cran()
  set.seed(5)
  d <- data.frame(x = stats::rnorm(300))
  d$y <- stats::rpois(300, exp(0.5 + 0.3 * d$x))
  for (fam in c("DPO", "HYPERPOIS", "GAMMACOUNT")) {
    fit <- countreg(y ~ x, data = d, family = fam, method = "BFGS")
    lr <- flexCountReg:::.fcr_lr_vs_poisson(fit)
    expect_false(is.null(lr), label = fam)
    expect_equal(lr$df, 1)
    expect_false(lr$boundary)
    expect_gte(lr$p.value, 0)
  }
  # Not registered as nested: Delaporte reaches Poisson only at rho = 0,
  # where k is not identified. The Hessian is then nearly singular; whether
  # it is singular in floating point (and so warns) depends on the platform,
  # so only the missing LR test is checked.
  fit <- suppressWarnings(
    countreg(y ~ x, data = d, family = "DELAPORTE", method = "BFGS"))
  expect_null(flexCountReg:::.fcr_lr_vs_poisson(fit))
})

test_that("underreporting models are available for thinning families", {
  skip_on_cran()
  set.seed(11)
  n <- 500
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  p <- stats::plogis(1 + d$z)
  ytrue <- rdelaporte(n, mu = exp(1 + 0.4 * d$x), rho = 0.5, k = 2)
  d$y <- stats::rbinom(n, ytrue, p)
  fit <- countreg(y ~ x, data = d, family = "DELAPORTE",
                  underreport_formula = ~ z, method = "BFGS")
  expect_true(isTRUE(fit$model$converged))
  expect_error(countreg(y ~ x, data = d, family = "DW",
                        underreport_formula = ~ z), "not available")
})

test_that("countreg.rp() fits the new families", {
  skip_on_cran()
  set.seed(3)
  n <- 300
  d <- data.frame(x = stats::rnorm(n))
  b <- 0.4 + 0.3 * stats::rnorm(n)
  d$y <- rdelaporte(n, mu = exp(0.6 + b * d$x), rho = 0.5, k = 2)
  for (fam in c("DELAPORTE", "PA", "DW")) {
    fit <- countreg.rp(y ~ 1, rpar_formula = ~ -1 + x, data = d,
                       family = fam, ndraws = 50, method = "BFGS")
    expect_true(isTRUE(fit$model$converged), label = fam)
    expect_true(all(is.finite(stats::coef(fit))), label = fam)
    pr <- stats::predict(fit, type = "response")
    expect_true(all(is.finite(pr) & pr > 0), label = fam)
  }
})
