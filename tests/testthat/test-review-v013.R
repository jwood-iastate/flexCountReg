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

test_that("F03: underreporting uses the distribution of thinned counts", {
  set.seed(21)
  n <- 4e5
  mu <- 4
  p <- 0.4
  h <- randtoolbox::halton(2000)
  check <- function(fam, y_true, alpha, sigma = NULL) {
    y_obs <- rbinom(n, y_true, p)
    emp <- tabulate(y_obs + 1, nbins = 13) / n
    thin <- .fcr_thin_params(fam, alpha, sigma, p)
    fit <- get_probFunc(fam)(y = 0:12, predicted = rep(p * mu, 13),
                             alpha = thin$alpha, sigma = thin$sigma,
                             haltons = h, normed_haltons = qnorm(h),
                             engine = "halton")
    expect_lt(max(abs(emp - fit)), 0.003, label = fam)
    # The old observation model (dispersion unchanged) is visibly different
    old <- get_probFunc(fam)(y = 0:12, predicted = rep(p * mu, 13),
                             alpha = alpha, sigma = sigma, haltons = h,
                             normed_haltons = qnorm(h), engine = "halton")
    expect_gt(max(abs(emp - old)), 0.01, label = paste(fam, "old"))
  }
  check("NB1", rnbinom(n, size = mu / 1.5, mu = mu), 1.5)
  check("NBP", rnbinom(n, size = mu^(2 - 1.5) / 0.8, mu = mu), 0.8, 1.5)
  check("PIG1", rpinvgaus(n, mu = mu, eta = 1.5), 1.5)
  # p = 1 recovers the original family
  expect_equal(.fcr_thin_params("NBP", 0.8, 1.5, 1)$alpha, 0.8)
})

test_that("F03: families not closed under thinning reject underreporting", {
  d <- washington_roads
  for (fam in c("GP1", "GP2", "COM", "GW")) {
    expect_error(countreg(Total_crashes ~ lnaadt, d, family = fam,
                          underreport_formula = ~ speed50),
                 "not of the same family")
  }
})

test_that("F03/F04: NB1 underreporting: thinned dispersion, identification", {
  set.seed(22)
  n <- 6000
  x <- rnorm(n)
  z <- rbinom(n, 1, 0.5)
  mu <- exp(1 + 0.4 * x)
  y_true <- rnbinom(n, size = mu / 1.5, mu = mu)
  y <- rbinom(n, y_true, plogis(-0.5 + 1.5 * z))
  d <- data.frame(y, x, z)
  # The count intercept and reporting intercept are nearly collinear (only
  # the logit curvature separates them): a warning, and the truth lies
  # within the (wide) uncertainty
  expect_warning(
    m <- countreg(y ~ x, d, family = "NB1", underreport_formula = ~ z,
                  method = "BFGS"),
    "weakly identified")
  # The well-identified slope is recovered
  expect_equal(unname(coef(m)["x"]), 0.4, tolerance = 0.1)
  expect_equal(unname(sqrt(vcov(m)["x", "x"])), 0.013, tolerance = 0.2)
  # Pearson residuals use the thinned variance
  expect_equal(var(residuals(m, type = "pearson")), 1, tolerance = 0.1)
  # A constant reporting probability is not identified at all
  expect_error(countreg(y ~ x, d, family = "NB1", underreport_formula = ~ 1),
               "intercept-only `underreport_formula` is not identified")
})

test_that("F05: no LR test against Poisson when Poisson is not nested", {
  d <- washington_roads
  pl <- countreg(Total_crashes ~ lnaadt, d, family = "PL", method = "BFGS")
  expect_null(summary(pl)$lr_test)
  expect_message(r <- regCompTest(pl), "not nested")
  expect_true(is.na(r$LR_pvalue))
  expect_true(is.finite(r$AICbase))
  # NB2 with a scalar dispersion: boundary test (chi-bar-square)
  nb <- countreg(Total_crashes ~ lnaadt, d, method = "BFGS")
  lr <- summary(nb)$lr_test
  expect_true(lr$boundary)
  expect_equal(lr$p.value, 0.5 * pchisq(lr$statistic, 1, lower.tail = FALSE))
  # With covariates in the base model too: 0.5 chi2(df - 1) + 0.5 chi2(df)
  r <- regCompTest(nb)
  expect_equal(r$LRdof, 2)
  expect_equal(r$LR_type, "boundary")
  expect_equal(r$LR_pvalue,
               0.5 * pchisq(r$LR, 1, lower.tail = FALSE) +
                 0.5 * pchisq(r$LR, 2, lower.tail = FALSE))
  # A dispersion formula makes the Poisson restriction nonregular
  nbf <- countreg(Total_crashes ~ lnaadt, d, method = "BFGS",
                  dis_param_formula_1 = ~ speed50)
  expect_null(summary(nbf)$lr_test)
  # GP1: interior restriction (phi = 0), ordinary chi-square
  gp <- countreg(Total_crashes ~ lnaadt, d, family = "GP1", method = "BFGS")
  lr <- summary(gp)$lr_test
  expect_false(lr$boundary)
  expect_equal(lr$p.value, pchisq(lr$statistic, 1, lower.tail = FALSE))
})

test_that("F05: two-model LR tests check nesting and the sample", {
  d <- washington_roads
  nb1 <- countreg(Total_crashes ~ lnaadt, d, family = "NB1", method = "BFGS")
  nbp <- countreg(Total_crashes ~ lnaadt, d, family = "NBP", method = "BFGS")
  pl <- countreg(Total_crashes ~ lnaadt, d, family = "PL", method = "BFGS")
  plg <- countreg(Total_crashes ~ lnaadt, d, family = "PLG", method = "BFGS")
  # NB1 is NBP with P = 1 (interior)
  r <- regCompTest(nb1, model2 = nbp)
  expect_equal(r$LR_type, "interior")
  expect_equal(r$LR_pvalue, pchisq(r$LR, 1, lower.tail = FALSE))
  # PL is PLG with alpha = 0 (boundary)
  expect_equal(.fcr_nesting(.fcr_model_spec(pl), .fcr_model_spec(plg))$type,
               "boundary")
  # Poisson within NB2, two-model form: boundary test
  po <- countreg(Total_crashes ~ lnaadt, d, family = "POISSON",
                 method = "BFGS")
  nb2 <- countreg(Total_crashes ~ lnaadt, d, method = "BFGS")
  r <- regCompTest(nb2, model2 = po)
  expect_equal(r$LR_type, "boundary")
  expect_equal(r$LR_pvalue, 0.5 * pchisq(r$LR, 1, lower.tail = FALSE))
  # PL and NB1 have the same number of parameters; PL vs NBP is not nested
  expect_warning(r <- regCompTest(pl, model2 = nbp), "not nested")
  expect_true(is.na(r$LR_pvalue))
  # Same number of rows but different observations
  d2 <- d
  d2$Total_crashes <- rev(d2$Total_crashes)
  other <- countreg(Total_crashes ~ lnaadt, d2, family = "NBP",
                    method = "BFGS")
  expect_error(regCompTest(nb1, model2 = other), "same observations")
})

test_that("F06: likelihoods use log-probabilities computed on the log scale", {
  # log P(Y = 1000 | mu = 1) is about -5913, but P itself underflows to 0
  for (fam in c("POISSON", "NB2", "NB1", "GP1", "COM", "PIG2", "PL")) {
    a <- if (fam == "POISSON") NULL else if (fam == "GP1") 1.2 else 0.5
    lp <- get_probFunc(fam)(y = 1000, predicted = 1, alpha = a, sigma = NULL,
                            haltons = NULL, normed_haltons = NULL,
                            engine = "halton", log = TRUE)
    expect_true(is.finite(lp) && lp < -500, label = fam)
  }
  expect_equal(get_probFunc("POISSON")(1000, 1, log = TRUE),
               dpois(1000, 1, log = TRUE))
  # A fit with an extreme count has a finite log-likelihood equal to glm()'s
  d <- washington_roads
  d$Total_crashes[1] <- 1000
  m <- countreg(Total_crashes ~ lnaadt, d, family = "POISSON",
                method = "BFGS")
  g <- glm(Total_crashes ~ lnaadt, poisson, d)
  expect_equal(m$model$maximum, as.numeric(logLik(g)), tolerance = 1e-6)
  # Zero-model adjustments on the log scale match the ordinary scale and stay
  # finite where 1 - P(Y = 0) is tiny
  y <- c(0, 0, 3, 5)
  lp <- dpois(y, 2, log = TRUE)
  lf0 <- rep(dpois(0, 2, log = TRUE), 4)
  eta <- c(-1, 0.5, 2, -3)
  for (type in c("inflated", "hurdle")) {
    expect_equal(.fcr_zero_adjust_log(y, lp, lf0, eta, type),
                 log(.fcr_zero_adjust(y, exp(lp), exp(lf0), plogis(eta),
                                      type)), tolerance = 1e-12)
  }
  # Hurdle with mean 1e-10: P(Y > 0) = 1e-10, and log P(Y = 1) is finite
  mu <- 1e-10
  h <- .fcr_zero_adjust_log(1, dpois(1, mu, log = TRUE),
                            dpois(0, mu, log = TRUE), 0, "hurdle")
  expect_equal(h, log(0.5) + dpois(1, mu, log = TRUE) -
                 log(-expm1(-mu)), tolerance = 1e-10)
})

test_that("F07: COM uses log normalizers, consistent support and a checked solver", {
  # nu = 1 is Poisson, also where Z = exp(1000) overflows
  expect_equal(dcom(1000, lambda = 1000, nu = 1), dpois(1000, 1000),
               tolerance = 1e-10)
  expect_equal(dcom(c(0, 500, 1000, 1500), lambda = 1000, nu = 1, log = TRUE),
               dpois(c(0, 500, 1000, 1500), 1000, log = TRUE),
               tolerance = 1e-10)
  expect_equal(pcom(c(900, 1100), lambda = 1000, nu = 1),
               ppois(c(900, 1100), 1000), tolerance = 1e-10)
  expect_equal(qcom(c(0.1, 0.5, 0.99), lambda = 1000, nu = 1),
               qpois(c(0.1, 0.5, 0.99), 1000))
  expect_equal(com_moments_vec_cpp(c(3, 1000), c(1, 1)),
               cbind(c(3, 1000), c(3, 1000)), tolerance = 1e-8,
               ignore_attr = TRUE)
  # The rate solver reaches the requested mean
  for (nu in c(0.05, 0.5, 1.7, 5)) {
    for (mu in c(1e-3, 0.5, 5, 100, 1000)) {
      lam <- find_lambda_cpp(mu, nu)
      expect_equal(com_expect_cpp(lam, nu), mu, tolerance = 1e-7,
                   label = paste(mu, nu))
    }
  }
  expect_true(is.nan(find_lambda_cpp(-1, 1)))
  # Moments agree with the PMF
  p <- dcom(0:200, mu = 4, nu = 0.6)
  expect_equal(sum(p), 1, tolerance = 1e-12)
  expect_equal(sum((0:200) * p), 4, tolerance = 1e-7)
  expect_equal(.fcr_family_variance("COM", 4, 0.6),
               sum((0:200)^2 * p) - 16, tolerance = 1e-6)
})

test_that("F08: Cholesky labels, starts and decoding use one packing", {
  for (K in 1:4) {
    set.seed(K)
    L <- matrix(0, K, K)
    L[lower.tri(L, diag = TRUE)] <- rnorm(K * (K + 1) / 2)
    v <- .fcr_chol_pack(L)
    expect_equal(.fcr_chol_unpack(v, K), L)
    rp <- letters[seq_len(K)]
    nm <- .fcr_chol_names(rp)
    # Each named element is the matching element of L
    for (k in seq_along(nm)) {
      ij <- strsplit(sub("^Chol:", "", nm[k]), "_")[[1]]
      expect_equal(v[k], L[match(ij[1], rp), match(ij[2], rp)])
    }
    expect_equal(.fcr_chol_unpack(.fcr_chol_pack(diag(0.1, K)), K),
                 diag(0.1, K))
  }
  expect_error(.fcr_chol_unpack(1:4, 3), "needs 6 values")
  # Three correlated random parameters: the default start is 0.1 * I, and
  # named start values are placed by name
  d <- washington_roads
  d$x3 <- d$ShouldWidth04
  rp3 <- ~ -1 + lnaadt + speed50 + x3
  m <- suppressWarnings(countreg.rp(Total_crashes ~ 1, rpar_formula = rp3,
                                    data = d, correlated = TRUE, ndraws = 20,
                                    method = "BFGS", max.iters = 0))
  ch <- coef(m)[grepl("^Chol:", names(coef(m)))]
  L0 <- .fcr_chol_unpack(unname(ch), 3)
  expect_equal(L0, diag(0.1, 3))
  expect_equal(names(ch)[c(1, 4, 6)],
               c("Chol:lnaadt_lnaadt", "Chol:speed50_speed50", "Chol:x3_x3"))
  # The simulated covariance of the draws is L L'
  L <- matrix(c(0.3, 0.1, -0.2, 0, 0.4, 0.05, 0, 0, 0.25), 3)
  h <- randtoolbox::halton(20000, 3)
  dr <- generate_random_draws(
    hdraws = h, random_coefs_means = c(0, 0, 0),
    rand_var_params = .fcr_chol_pack(L), rpardists = rep("n", 3),
    rpar = c("a", "b", "c"), X_rand = diag(3), correlated = TRUE)
  expect_equal(cov(t(dr$xb_rand_mat)), L %*% t(L), tolerance = 0.02)
})

test_that("F10: a free-location random slope cannot also be a fixed slope", {
  d <- washington_roads
  expect_error(countreg.rp(Total_crashes ~ lnaadt + speed50,
                           rpar_formula = ~ -1 + speed50, data = d,
                           rpardists = c(speed50 = "n")),
               "Only the sum of the fixed coefficient")
  expect_error(countreg.rp(Total_crashes ~ lnaadt + speed50,
                           rpar_formula = ~ -1 + speed50, data = d,
                           rpardists = c(speed50 = "u")),
               "both `formula` and `rpar_formula`")
  # Gamma random parameters have a fixed support, so a fixed shift is allowed
  expect_s3_class(
    suppressWarnings(countreg.rp(Total_crashes ~ lnaadt + speed50,
                                 rpar_formula = ~ -1 + speed50, data = d,
                                 rpardists = c(speed50 = "g"), ndraws = 20,
                                 method = "BFGS", max.iters = 20)),
    "countreg_rp")
})

test_that("F12: negmultinom predictions use the training blueprint", {
  d <- washington_roads
  d$fac <- factor(ifelse(d$speed50 == 1, "fast", "slow"))
  m <- negmultinom(Total_crashes ~ lnaadt + fac, "ID", d,
                   offset = "lnlength")
  # A single row, a factor subset and reordered rows predict the same values
  all_p <- predict(m, newdata = d)
  expect_equal(predict(m, newdata = d[5, ]), all_p[5])
  slow <- d[d$fac == "slow", ][1:3, ]
  expect_equal(predict(m, newdata = slow), all_p[as.numeric(rownames(slow))])
  o <- c(10, 3, 7)
  expect_equal(predict(m, newdata = d[o, ]), all_p[o])
  # A missing offset column is an error, not a silently dropped offset
  expect_error(predict(m, newdata = d[1:3, setdiff(names(d), "lnlength")]))
  # Formula and argument offsets agree; both together are an error
  m2 <- negmultinom(Total_crashes ~ lnaadt + fac + offset(lnlength), "ID", d)
  expect_equal(coef(m2), coef(m), tolerance = 1e-6)
  expect_error(negmultinom(Total_crashes ~ lnaadt + offset(lnlength), "ID", d,
                           offset = "lnlength"), "not both")
  # EB predictions on the training data match the fitted EB values
  expect_equal(predict(m, newdata = d, method = "EB"), m$model$eb.fitted,
               tolerance = 1e-8)
  # Named start values are matched by name
  sv <- rev(coef(m))
  m3 <- negmultinom(Total_crashes ~ lnaadt + fac, "ID", d,
                    offset = "lnlength", start.vals = sv)
  expect_equal(coef(m3), coef(m), tolerance = 1e-5)
})

test_that("F13: panel-model marginal-effect SEs follow the delta method", {
  d <- washington_roads
  for (fit in list(
    renb(Total_crashes ~ lnaadt + speed50, d, group_var = "ID",
         method = "BFGS"),
    poisLind.re(Total_crashes ~ lnaadt + speed50, d, group_var = "ID",
                method = "BFGS"),
    repois(Total_crashes ~ lnaadt + speed50, d, group_var = "ID",
           method = "BFGS"))) {
    me <- margEffTable(fit, digits = 8)
    b <- coef(fit)[["speed50"]]
    se_b <- sqrt(vcov(fit)["speed50", "speed50"])
    expect_equal(me$std_error[me$term == "speed50"],
                 100 * exp(b) * se_b, tolerance = 1e-4,
                 label = class(fit)[1])
    expect_gt(me$std_error[me$term == "lnaadt"], 0)
  }
})

test_that("F24: negmultinom convergence uses the optimizer's codes", {
  m <- suppressWarnings(negmultinom(Total_crashes ~ lnaadt, "ID",
                                    washington_roads, method = "BFGS",
                                    max.iters = 2))
  expect_false(m$model$converged)
  m <- negmultinom(Total_crashes ~ lnaadt, "ID", washington_roads)
  expect_true(m$model$converged)
  expect_equal(m$model$bic, BIC(m))
})

test_that("F11: auxiliary formulas keep their trained data-dependent bases", {
  d <- washington_roads
  m <- countreg(Total_crashes ~ lnaadt, d, offset = "lnlength",
                dis_param_formula_1 = ~ poly(AADT, 2), method = "BFGS")
  p_all <- predict(m, newdata = d, type = "prob", max_count = 3)
  # One row alone, or in a different batch, gives the same prediction
  for (i in c(1, 50, 200)) {
    expect_equal(predict(m, newdata = d[i, ], type = "prob", max_count = 3),
                 p_all[i, , drop = FALSE], ignore_attr = TRUE)
  }
  o <- c(200, 3, 77, 15)
  expect_equal(predict(m, newdata = d[o, ], type = "prob", max_count = 3),
               p_all[o, ], ignore_attr = TRUE)
  # Pearson variances use the same dispersion as in estimation
  expect_equal(residuals(m, type = "pearson"),
               as.numeric(d$Total_crashes - fitted(m)) /
                 sqrt(.fcr_fitted_variance(m)))
  # scale() in a zero model
  z <- countreg(Animal ~ lnaadt, d, family = "POISSON",
                zero_formula = ~ scale(AADT), method = "NR")
  zp <- predict(z, newdata = d, type = "zero")
  expect_equal(predict(z, newdata = d[c(9, 2), ], type = "zero"), zp[c(9, 2)])
})

test_that("F14: the marginal-effect bootstrap resamples independent units", {
  d <- washington_roads
  # A cluster drawn twice becomes two panels with fresh ids
  set.seed(1)
  rs <- .margEff_resample_data(d[1:40, ], cluster_var = "ID")
  ids <- split(rs$data$ID, rs$data$.fcr_boot_id)
  expect_true(all(vapply(ids, function(v) length(unique(v)) == 1, NA)))
  expect_equal(sum(lengths(ids)), length(rs$idx))
  expect_true(any(duplicated(vapply(ids, function(v) as.character(v[1]), ""))))
  # Panel models resample panels by default, and each refit sees one panel
  # per sampled cluster
  re <- renb(Total_crashes ~ lnaadt + speed50, d, group_var = "ID",
             method = "BFGS")
  calls <- list()
  trace_refit <- function(object, boot_data, panel = NULL, weights = NULL) {
    calls[[length(calls) + 1]] <<- list(panel = panel,
                                        n = length(unique(boot_data[[panel]])))
    stop("skip")
  }
  local_mocked_bindings(.margEff_refit = trace_refit)
  suppressWarnings(margEffTable(re, vars = "lnaadt", se = "bootstrap",
                                bootstraps = 3))
  expect_length(calls, 3)
  expect_equal(calls[[1]]$panel, ".fcr_boot_id")
  expect_equal(calls[[1]]$n, length(unique(d$ID)))
  expect_error(margEffTable(re, se = "bootstrap", cluster_var = "nope"),
               "not found")
})

test_that("F14: refits resample numeric weights and do not bootstrap again", {
  d <- washington_roads
  w <- rep(1:2, length.out = nrow(d))
  m <- countreg(Total_crashes ~ lnaadt, d, weights = w, method = "BFGS",
                stderr = "boot", bootstraps = 5)
  set.seed(2)
  rs <- .margEff_resample_data(d)
  r <- .margEff_refit(m, rs$data, weights = w[rs$idx])
  ref <- countreg(Total_crashes ~ lnaadt, rs$data, weights = w[rs$idx],
                  method = "BFGS")
  expect_equal(coef(r), coef(ref), tolerance = 1e-8)
  expect_equal(r$model$vcov_type, "hessian")
})

test_that("F15: bootstrap replicates must converge and be finite", {
  set.seed(1)
  k <- 0
  refit <- function(idx) {
    k <<- k + 1
    est <- c(a = mean(idx), b = 1)
    if (k %% 4 == 1) attr(est, "converged") <- FALSE
    if (k %% 4 == 2) est["b"] <- Inf
    est
  }
  expect_warning(dr <- .fcr_bootstrap(10, 12, c("a", "b"), refit),
                 "6 of 12 bootstrap replicates failed \\(not converged: 3, non-finite estimates: 3\\)")
  expect_equal(sum(stats::complete.cases(dr)), 6)
  expect_equal(as.vector(attr(dr, "failures")[c("not converged",
                                                 "non-finite estimates")]),
               c(3, 3))
  # A refit stopped at the iteration limit is not accepted
  d <- washington_roads
  expect_error(suppressWarnings(
    countreg(Total_crashes ~ lnaadt, d, method = "BFGS", stderr = "boot",
             bootstraps = 4, max.iters = 1)),
    "Fewer than 2 bootstrap replicates succeeded")
})

test_that("F16: frequency weights match the expanded data", {
  d <- washington_roads[1:300, ]
  w <- rep(1:3, length.out = nrow(d))
  ex <- d[rep(seq_len(nrow(d)), w), ]
  a <- countreg(Total_crashes ~ lnaadt, d, weights = w, method = "BFGS",
                stderr = "robust")
  b <- countreg(Total_crashes ~ lnaadt, ex, method = "BFGS",
                stderr = "robust")
  expect_equal(coef(a), coef(b), tolerance = 1e-6)
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)),
               tolerance = 1e-8)
  expect_equal(vcov(a), vcov(b), tolerance = 1e-4)
  expect_equal(BIC(a), BIC(b), tolerance = 1e-8)
  # Without weights the robust covariance is sandwich::sandwich()
  u <- countreg(Total_crashes ~ lnaadt, d, method = "BFGS", stderr = "robust")
  expect_equal(vcov(u), sandwich::sandwich(u$model), tolerance = 1e-6,
               ignore_attr = TRUE)
  # Panel frequency weights (constant within panels): panels are counted
  # with their weights in BIC
  wp <- as.numeric(factor(d$ID)) %% 2 + 1
  rp <- countreg.rp(Total_crashes ~ lnaadt, rpar_formula = ~ -1 + speed50,
                    data = d, weights = wp, group_var = "ID", ndraws = 20,
                    method = "BFGS", max.iters = 50)
  expect_equal(attr(logLik(rp), "nobs"), sum(wp[!duplicated(d$ID)]))
})

test_that("F17: hurdle simulation never puts zeros in the positive part", {
  set.seed(5)
  n <- 4000
  d <- data.frame(x = rnorm(n))
  # Tiny count means: P(Y > 0) of the count part is about 1e-3
  d$y <- ifelse(runif(n) < 0.4, 1 + rpois(n, 0.001), 0)
  m <- countreg(y ~ 1, d, family = "POISSON", zero_formula = ~ 1,
                zero_type = "hurdle", method = "BFGS")
  p0 <- predict(m, type = "zero")   # P(Y = 0) of the hurdle model
  sims <- simulate(m, nsim = 20, seed = 1)
  # The zero share matches P(Y = 0)
  expect_equal(mean(unlist(sims) == 0), mean(p0), tolerance = 0.02)
  # Positive draws from the truncated Poisson with a tiny mean are almost
  # all 1, and their mean matches the truncated mean
  pos <- unlist(sims)[unlist(sims) > 0]
  mu <- exp(coef(m)[["(Intercept)"]])
  expect_equal(mean(pos), mu / -expm1(-mu), tolerance = 0.01)
  # Mean predictions use the precise 1 - P(Y = 0) of the count part
  expect_equal(predict(m)[1], (1 - p0[1]) * mu / -expm1(-mu),
               tolerance = 1e-8)
})

test_that("F18: every count p-function follows the same contract", {
  fs <- list(
    pplind = function(q, ...) pplind(q, mean = 2, theta = 1.5, ...),
    ppinvgaus = function(q, ...) ppinvgaus(q, mu = 2, eta = 0.5, ...),
    ppinvgamma = function(q, ...) ppinvgamma(q, mu = 2, eta = 0.5, ...),
    ppLnorm = function(q, ...) ppLnorm(q, mean = 2, sigma = 0.5, ...),
    pplindLnorm = function(q, ...)
      pplindLnorm(q, mean = 2, theta = 1.5, sigma = 0.5, ndraws = 200, ...),
    pplindGamma = function(q, ...)
      pplindGamma(q, mean = 2, theta = 1.5, alpha = 0.4, ...),
    ppge = function(q, ...) ppge(q, mean = 2, shape = 1.5, scale = 1, ...),
    ppoisweibull = function(q, ...)
      ppoisweibull(q, lambda = 2, alpha = 1.5, sigma = 1, ...),
    psichel = function(q, ...) psichel(q, mu = 2, sigma = 1, gamma = -0.5,
                                       ...),
    pgwar = function(q, ...) pgwar(q, mu = 2, k = 2, rho = 3, ...),
    pcom = function(q, ...) pcom(q, mu = 2, nu = 0.7, ...),
    pgp1 = function(q, ...) pgp1(q, mu = 2, phi = 0.3, ...),
    pgp2 = function(q, ...) pgp2(q, mu = 2, alpha = 0.3, ...))
  for (nm in names(fs)) {
    f <- fs[[nm]]
    p <- f(c(NA, -Inf, -1, Inf))
    expect_true(is.na(p[1]), label = nm)
    expect_equal(p[2:4], c(0, 0, 1), label = nm)
    expect_equal(f(c(-1, Inf), lower.tail = FALSE), c(1, 0), label = nm)
    expect_equal(f(-1, log.p = TRUE), -Inf, label = nm)
    lo <- f(0:6)
    expect_true(all(diff(lo) >= 0), label = paste(nm, "monotone"))
    expect_equal(lo + f(0:6, lower.tail = FALSE), rep(1, 7),
                 tolerance = 1e-6, label = paste(nm, "tails"))
    expect_identical(f(numeric(0)), numeric(0), label = nm)
  }
  # Small upper tails are summed directly (1 - CDF would be 0 or noise)
  expect_equal(pplind(60, mean = 2, theta = 1.5, lower.tail = FALSE,
                      log.p = TRUE),
               log(sum(dplind(61:3000, mean = 2, theta = 1.5))),
               tolerance = 1e-8)
  expect_equal(pgp1(40, mu = 2, phi = 0.3, lower.tail = FALSE),
               sum(dgp1(41:2000, 2, 0.3)), tolerance = 1e-8)
  expect_equal(ppinvgaus(80, mu = 2, eta = 0.5, lower.tail = FALSE,
                         log.p = TRUE),
               log(sum(dpinvgaus(81:3000, mu = 2, eta = 0.5))),
               tolerance = 1e-8)
})

test_that("F19: CURE and augment(data = ) keep rows aligned", {
  d <- washington_roads
  d$lnaadt[3] <- NA
  m <- suppressMessages(countreg(Total_crashes ~ lnaadt, d, method = "BFGS"))
  # augment() with reordered data predicts in that order
  o <- sample(nrow(d))
  a <- augment(m, data = d[o, ])
  keep <- !is.na(d$lnaadt[o])
  expect_equal(nrow(a), sum(keep))
  expect_equal(a$.fitted, as.vector(predict(m, newdata = d[o, ][keep, ])))
  expect_equal(a$.resid, d$Total_crashes[o][keep] - a$.fitted)
  ap <- augment(m, data = d[o, ], type.residuals = "pearson")
  expect_equal(ap$.resid[order(o[keep])],
               as.numeric(residuals(m, type = "pearson")))
  # CURE uses one mask for outcome, predictions and covariate
  skip_if_not_installed("cureplots")
  d$speed50[10] <- NA
  expect_s3_class(cureplot(m, data = d, indvar = "speed50"), "ggplot")
  expect_error(cureplot(m, indvar = "nope"), "not found")
  expect_error(cureplot(m, method = "Exact"), "applies only")
})
