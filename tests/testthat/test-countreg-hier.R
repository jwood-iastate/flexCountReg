# Hierarchical count models: countreg_hier().

sal <- function() {
  data("stream_salamanders", package = "flexCountReg", envir = environment())
  s <- stream_salamanders
  s$sample <- factor(s$sample)
  s
}

test_that("the Gauss-Hermite rule integrates normal moments exactly", {
  r <- .fcr_gh_rule(7)
  expect_equal(sum(r$w), 1)
  expect_equal(sum(r$w * r$z^2), 1)
  expect_equal(sum(r$w * r$z^4), 3)
  expect_equal(sum(r$w * r$z^12), 10395)   # 11!!, exact up to degree 13
  g <- .fcr_gh_grid(3, 2)
  expect_equal(nrow(g$Z), 9)
  expect_equal(sum(g$w * g$Z[, 1]^2 * g$Z[, 2]^2), 1)
})

test_that("analytic eta derivatives match differences of the log-density", {
  set.seed(1)
  y <- rpois(100, 3)
  eta <- rnorm(100, 1, 0.5)
  a <- rep(0.7, 100)
  h <- 1e-5
  for (fam in c("POISSON", "NB1", "NB2")) {
    f <- .fcr_hier_logf(fam, y, rep(1, 100))
    lf <- function(e) {
      .fcr_get_family(fam)$dens(y, exp(e), a, NULL, log = TRUE)
    }
    d <- f(eta, a, NULL)
    expect_equal(d$l, lf(eta), tolerance = 1e-10)
    expect_equal(d$d1, (lf(eta + h) - lf(eta - h)) / (2 * h),
                 tolerance = 1e-6)
    expect_equal(d$d2, (f(eta + h, a, NULL)$d1 - f(eta - h, a, NULL)$d1) /
                   (2 * h), tolerance = 1e-6)
  }
  # The generic (difference) branch agrees with the analytic NB2
  nb <- .fcr_hier_logf("NB2", y, rep(2, 100))(eta, a, NULL)
  lf <- function(e) 2 * dnbinom(y, size = 1 / a, mu = exp(e), log = TRUE)
  h <- 1e-3
  d2 <- (-lf(eta + 2 * h) + 16 * lf(eta + h) - 30 * lf(eta) +
           16 * lf(eta - h) - lf(eta - 2 * h)) / (12 * h^2)
  expect_equal(nb$d2, d2, tolerance = 1e-7)
})

test_that("theta and (sd, cor) convert both ways", {
  terms <- list(list(q = 2L, ntheta = 3L, toff = 0L, cnms = c("a", "b"),
                     group = "g"),
                list(q = 1L, ntheta = 1L, toff = 3L, cnms = "(Intercept)",
                     group = "h"))
  th <- c(0.6, -0.2, 0.3, 0.4)
  vc <- .fcr_hier_vc_vector(terms, th)
  expect_named(vc, c("sd_a|g", "sd_b|g", "cor_b.a|g", "sd_(Intercept)|h"))
  expect_equal(vc[["sd_a|g"]], 0.6)
  expect_equal(vc[["cor_b.a|g"]], -0.2 / sqrt(0.2^2 + 0.3^2))
  expect_equal(.fcr_hier_theta_from_vc(terms, vc), th)
  # A negative diagonal gives the same covariance
  expect_equal(.fcr_hier_vc_vector(terms, c(-0.6, 0.2, 0.3, -0.4)), vc)
})

test_that("Poisson random intercepts match lme4::glmer (AGHQ and Laplace)", {
  skip_if_not_installed("lme4")
  data("washington_roads", package = "flexCountReg", envir = environment())
  f <- Total_crashes ~ lnaadt + speed50 + offset(lnlength) + (1 | ID)
  h9 <- countreg_hier(f, data = washington_roads, family = "POISSON")
  g9 <- suppressWarnings(lme4::glmer(f, data = washington_roads,
                                     family = poisson, nAGQ = 9))
  expect_equal(h9$model$engine, "aghq")
  expect_equal(h9$model$maximum, as.numeric(logLik(g9)), tolerance = 1e-6)
  expect_equal(unname(fixef(h9)), unname(lme4::fixef(g9)), tolerance = 1e-4)
  expect_equal(unname(coef(h9)["sd_(Intercept)|ID"]),
               sqrt(lme4::VarCorr(g9)$ID[1]), tolerance = 1e-4)
  h1 <- countreg_hier(f, data = washington_roads, family = "POISSON",
                      engine = "laplace")
  g1 <- suppressWarnings(lme4::glmer(f, data = washington_roads,
                                     family = poisson))
  expect_equal(h1$model$maximum, as.numeric(logLik(g1)), tolerance = 1e-6)
  expect_equal(unname(fixef(h1)), unname(lme4::fixef(g1)), tolerance = 1e-4)
  expect_equal(unname(h1$model$se[1:3]), unname(sqrt(diag(as.matrix(vcov(g1))))),
               tolerance = 5e-3)
  # One quadrature node is the Laplace approximation
  ha <- countreg_hier(f, data = washington_roads, family = "POISSON",
                      engine = "aghq", nAGQ = 1)
  expect_equal(ha$model$maximum, h1$model$maximum, tolerance = 1e-8)
})

test_that("NB1, NB2, crossed, nested and slope models match glmmTMB", {
  skip_if_not_installed("glmmTMB")
  skip_on_cran()
  s <- sal()
  tmb <- function(f, fam) {
    suppressWarnings(glmmTMB::glmmTMB(f, data = s, family = fam))
  }
  sds <- function(m) {
    unlist(lapply(glmmTMB::VarCorr(m)$cond, attr, "stddev"))
  }
  # NB2 random intercept
  f <- count ~ spp + mined + (1 | site)
  h <- countreg_hier(f, data = s, family = "NB2", engine = "laplace")
  m <- tmb(f, glmmTMB::nbinom2)
  expect_equal(h$model$maximum, as.numeric(logLik(m)), tolerance = 1e-7)
  expect_equal(unname(coef(h)),
               unname(c(glmmTMB::fixef(m)$cond, -log(sigma(m)), sds(m))),
               tolerance = 1e-4)
  expect_equal(unname(h$model$se[1:8]),
               unname(sqrt(diag(vcov(m)$cond))), tolerance = 5e-3)
  # NB1 (variance mu (1 + alpha), as glmmTMB's nbinom1)
  h <- countreg_hier(f, data = s, family = "NB1", engine = "laplace")
  m <- tmb(f, glmmTMB::nbinom1)
  expect_equal(h$model$maximum, as.numeric(logLik(m)), tolerance = 1e-7)
  expect_equal(unname(coef(h)),
               unname(c(glmmTMB::fixef(m)$cond, log(sigma(m)), sds(m))),
               tolerance = 1e-3)
  # Crossed site and sample effects
  f <- count ~ spp + mined + (1 | site) + (1 | sample)
  h <- countreg_hier(f, data = s, family = "POISSON")
  expect_equal(h$model$engine, "laplace")
  m <- tmb(f, poisson)
  expect_equal(h$model$maximum, as.numeric(logLik(m)), tolerance = 1e-7)
  expect_equal(unname(coef(h)),
               unname(c(glmmTMB::fixef(m)$cond, sds(m))), tolerance = 1e-3)
  # Uncorrelated intercept and slope (double bar)
  f <- count ~ cover + mined + (cover || site)
  h <- suppressMessages(countreg_hier(f, data = s, family = "NB2"))
  m <- tmb(count ~ cover + mined + (1 | site) + (0 + cover | site),
           glmmTMB::nbinom2)
  expect_equal(h$model$maximum, as.numeric(logLik(m)), tolerance = 1e-6)
  expect_length(h$model$re_terms, 2)
  # Correlated intercept and slope (Laplace and AGHQ with one node agree)
  f <- count ~ cover + mined + (1 + cover | site)
  h <- suppressMessages(suppressWarnings(
    countreg_hier(f, data = s, family = "POISSON", engine = "laplace")))
  ha <- suppressMessages(suppressWarnings(
    countreg_hier(f, data = s, family = "POISSON", nAGQ = 1)))
  m <- tmb(f, poisson)
  expect_equal(h$model$maximum, as.numeric(logLik(m)), tolerance = 1e-6)
  expect_equal(ha$model$maximum, h$model$maximum, tolerance = 1e-6)
})

test_that("nested terms expand and agree with the explicit form", {
  s <- sal()
  h1 <- suppressMessages(countreg_hier(count ~ mined + (1 | site/sample),
                                       data = s, family = "POISSON"))
  h2 <- suppressMessages(countreg_hier(
    count ~ mined + (1 | site) + (1 | site:sample), data = s,
    family = "POISSON"))
  expect_equal(h1$model$maximum, h2$model$maximum, tolerance = 1e-8)
  expect_equal(sort(unname(h1$model$n_groups)), c(23, 92))
})

test_that("estimates and predictions do not depend on row order", {
  s <- sal()
  f <- count ~ spp + mined + (1 | site)
  h <- countreg_hier(f, data = s, family = "NB2")
  set.seed(3)
  o <- sample(nrow(s))
  h2 <- countreg_hier(f, data = s[o, ], family = "NB2")
  expect_equal(coef(h2), coef(h), tolerance = 1e-5)
  expect_equal(predict(h2, newdata = s), predict(h, newdata = s),
               tolerance = 1e-5)
})

test_that("predictions condition on the chosen levels", {
  s <- sal()
  h <- countreg_hier(count ~ spp + mined + (1 | site) + (1 | sample),
                     data = s, family = "POISSON")
  expect_equal(predict(h), fitted(h))
  re <- ranef(h)
  expect_named(re, c("site", "sample"))
  eta0 <- predict(h, re.form = NA, type = "link")
  eta_site <- predict(h, re.form = ~ (1 | site), type = "link")
  expect_equal(eta_site - eta0,
               re$site[as.character(s$site), "(Intercept)"])
  eta_all <- predict(h, type = "link")
  expect_equal(eta_all - eta_site,
               re$sample[as.character(s$sample), "(Intercept)"])
  # Marginal mean over new groups: exp(eta + (var_site + var_sample) / 2)
  v <- unname(coef(h)[c("sd_(Intercept)|site", "sd_(Intercept)|sample")]^2)
  expect_equal(predict(h, re.form = NA, marginal = TRUE),
               exp(eta0 + sum(v) / 2))
  # New groups
  nd <- s[1:2, ]
  nd$site <- factor(c("R-1", "zzz"), levels = c(levels(s$site), "zzz"))
  p <- predict(h, newdata = nd, re.form = ~ (1 | site), type = "link")
  expect_equal(p[2], predict(h, newdata = nd[2, ], re.form = NA,
                             type = "link"))
  expect_error(predict(h, newdata = nd, allow.new.levels = FALSE),
               "not seen in estimation")
  expect_error(predict(h, re.form = ~ (1 | year)), "not in the model")
})

test_that("methods: summary, VarCorr, simulate, tidy and the LR test", {
  s <- sal()
  h <- countreg_hier(count ~ spp + mined + (1 | site), data = s,
                     family = "NB2")
  expect_s3_class(h, c("countreg_hier", "flexCountReg"))
  expect_output(print(summary(h)), "Random effects")
  vc <- VarCorr(h)
  expect_equal(attr(vc$site, "stddev")[[1]],
               unname(coef(h)["sd_(Intercept)|site"]))
  expect_output(print(vc), "site")
  expect_length(fixef(h), 8)
  expect_equal(dim(attr(ranef(h)$site, "postVar")), c(1, 1, 23))
  sim <- simulate(h, nsim = 3, seed = 1)
  expect_equal(dim(sim), c(nrow(s), 3))
  expect_true(all(sim >= 0))
  sim2 <- simulate(h, nsim = 3, seed = 1, use_re = "new")
  expect_false(isTRUE(all.equal(sim, sim2)))
  td <- tidy(h)
  expect_equal(td$component[td$term == "sd_(Intercept)|site"],
               "random_effects")
  expect_equal(stats::nobs(h), nrow(s))
  expect_equal(attr(logLik(h), "df"), 10)
  lr <- re_lrtest(h)
  expect_equal(lr$df, 1)
  expect_gt(lr$statistic, 10)
  m0 <- countreg(count ~ spp + mined, data = s, family = "NB2",
                 method = "BFGS")
  expect_equal(lr$logLik_without, m0$model$maximum, tolerance = 1e-5)
})

test_that("a variance at zero is reported as a singular fit", {
  s <- sal()
  expect_message(
    h <- countreg_hier(count ~ spp + mined + (1 | site) + (1 | sample),
                       data = s, family = "NB2"),
    "Singular fit")
  expect_true(h$model$singular)
  expect_lt(coef(h)[["sd_(Intercept)|sample"]], 1e-3)
  expect_true(is.na(h$model$se[["sd_(Intercept)|sample"]]))
  expect_output(print(summary(h)), "Singular fit")
})

test_that("other families, distribution formulas and weights", {
  s <- sal()
  hp <- countreg_hier(count ~ mined + (1 | site), data = s, family = "PIG2")
  expect_true(hp$model$converged)
  expect_true("ln(eta)" %in% names(coef(hp)))
  hd <- countreg_hier(count ~ mined + (1 | site), data = s, family = "NB2",
                      dis_param_formula_1 = ~ mined)
  expect_true("ln(alpha):minedyes" %in% names(coef(hd)))
  # Frequency weights equal duplicated rows
  w <- rep(1:2, length.out = nrow(s))
  hw <- countreg_hier(count ~ mined + (1 | site), data = s,
                      family = "POISSON", weights = w)
  hdup <- countreg_hier(count ~ mined + (1 | site),
                        data = s[rep(seq_len(nrow(s)), w), ],
                        family = "POISSON")
  expect_equal(hw$model$maximum, hdup$model$maximum, tolerance = 1e-7)
})

test_that("input errors are clear", {
  s <- sal()
  expect_error(countreg_hier(count ~ mined, data = s), "no random-effect")
  expect_error(countreg_hier(count ~ mined + (1 | site) + (1 | sample),
                             data = s, engine = "aghq"),
               "single grouping term")
  expect_error(countreg_hier(count ~ mined + (1 | nosuch), data = s),
               "not found")
  expect_error(countreg_hier(count ~ te(cover, DOP) + (1 | site), data = s),
               "use s\\(\\) or t2")
  expect_error(countreg_hier(count ~ cover + (1 | site), data = s,
                             dis_param_formula_1 = ~ s(cover)),
               "distribution-parameter formulas")
  expect_error(countreg_hier(count ~ cover + (1 + cover | site), data = s,
                             re_dist = "gamma"), "one effect")
  expect_error(countreg_hier(count ~ cover + (1 | site) + (1 | sample),
                             data = s, re_dist = "uniform"),
               "single grouping term")
  expect_error(countreg_hier(count ~ cover + (1 | site), data = s,
                             re_dist = "lognormal"), "must be one of")
})

# ---- Phase 2 ---------------------------------------------------------------

test_that("smooth terms match mgcv's marginal-likelihood fit", {
  skip_if_not_installed("mgcv")
  set.seed(2)
  d <- data.frame(g = factor(sample(40, 800, TRUE)), x = runif(800))
  d$y <- rpois(800, exp(0.3 + sin(2 * pi * d$x) + rnorm(40, 0, 0.5)[d$g]))
  h <- countreg_hier(y ~ s(x) + (1 | g), data = d, family = "POISSON")
  gm <- mgcv::gam(y ~ s(x) + s(g, bs = "re"), data = d, family = poisson,
                  method = "ML")
  st <- h$model$smooth_table
  expect_equal(st$edf, unname(summary(gm)$edf[1]), tolerance = 2e-3)
  expect_equal(h$model$smooth_lambda, unname(gm$sp[1]), tolerance = 1e-2)
  expect_equal(unname(coef(h)["sd_(Intercept)|g"]),
               unname(sqrt(gm$sig2 / gm$sp[2])), tolerance = 2e-3)
  se <- smooth_estimates(h, n = 6, seed = 1)
  pg <- predict(gm, newdata = data.frame(x = se$x, g = "1"), type = "terms",
                se.fit = TRUE)
  expect_equal(se$fit, unname(pg$fit[, "s(x)"]), tolerance = 5e-3)
  expect_equal(se$se, unname(pg$se.fit[, "s(x)"]), tolerance = 5e-3)
  # Predictions use the smooth on its original basis
  nd <- d[1:5, ]
  expect_equal(predict(h, newdata = nd), fitted(h)[1:5], tolerance = 1e-8)
  expect_output(print(summary(h)), "Smooth terms")
  pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_silent(plot_smooths(h))
  # A smooth alone (no grouping term) is a GAM fitted by marginal likelihood
  h0 <- countreg_hier(y ~ s(x), data = d, family = "POISSON")
  expect_length(h0$model$re_terms, 0)
  expect_gt(h0$model$smooth_table$edf, 3)
  lr <- re_lrtest(h)
  expect_equal(lr$logLik_without, h0$model$maximum, tolerance = 1e-5)
})

test_that("the simulation engine matches quadrature and closed forms", {
  s <- sal()
  f <- count ~ spp + mined + (1 | site)
  ha <- countreg_hier(f, data = s, family = "POISSON")
  hn <- countreg_hier(f, data = s, family = "POISSON", engine = "msl",
                      re_draws = 1000)
  expect_equal(hn$model$engine, "msl")
  expect_equal(hn$model$maximum, ha$model$maximum, tolerance = 1e-4)
  expect_equal(coef(hn), coef(ha), tolerance = 0.02)
  # Poisson with gamma group effects is repois()
  hg <- countreg_hier(f, data = s, family = "POISSON", re_dist = "gamma",
                      re_draws = 1000)
  rp <- repois(count ~ spp + mined, group_var = "site", data = s,
               method = "BFGS")
  expect_equal(hg$model$maximum, rp$model$maximum, tolerance = 1e-4)
  expect_equal(unname(coef(hg)[1:8]), unname(coef(rp)[1:8]),
               tolerance = 0.02)
  expect_equal(unname(coef(hg)["sd_(Intercept)|site"]),
               unname(sqrt(exp(coef(rp)["lnalpha"]))), tolerance = 0.02)
  # Gamma effects have mean 1 on the count scale
  expect_equal(predict(hg, re.form = NA, marginal = TRUE),
               predict(hg, re.form = NA))
  # Uniform effects: marginal mean from the uniform's mgf
  hu <- countreg_hier(f, data = s, family = "POISSON", re_dist = "uniform",
                      re_draws = 300)
  sdu <- coef(hu)[["sd_(Intercept)|site"]]
  a <- sqrt(3) * sdu
  expect_equal(predict(hu, re.form = NA, marginal = TRUE),
               predict(hu, re.form = NA) * sinh(a) / a)
  expect_equal(dim(attr(ranef(hu)$site, "postVar")), c(1, 1, 23))
  expect_output(print(hu), "uniform")
  sim <- simulate(hu, 2, seed = 1, use_re = "new")
  expect_equal(dim(sim), c(nrow(s), 2))
})

test_that("shifted Halton draws are evenly spread and reproducible", {
  U <- .fcr_hier_draws(50, 200, 2, seed = 7)
  expect_length(U, 2)
  expect_equal(dim(U[[1]]), c(50, 200))
  expect_true(all(U[[1]] >= 0 & U[[1]] < 1))
  # Each group's draws are a rotated low-discrepancy set
  expect_lt(max(abs(rowMeans(U[[1]]) - 0.5)), 0.01)
  expect_identical(.fcr_hier_draws(50, 200, 2, seed = 7), U)
  set.seed(1)
  a <- runif(1)
  set.seed(1)
  .fcr_hier_draws(5, 10, 1, seed = 99)
  expect_identical(runif(1), a)
})

test_that("predictions from a separate history", {
  data("washington_roads", package = "flexCountReg", envir = environment())
  w <- washington_roads
  past <- w[w$Year < 2018, ]
  fut <- w[w$Year == 2018, ]
  h <- countreg_hier(Total_crashes ~ lnaadt + speed50 + offset(lnlength) +
                       (1 | ID), data = past, family = "POISSON")
  # The training data as history gives the fitted effects
  expect_equal(predict(h, newdata = fut, history = past),
               predict(h, newdata = fut), tolerance = 1e-6)
  one <- past[past$Year == min(past$Year), ]
  p1 <- predict(h, newdata = fut, history = one)
  expect_false(isTRUE(all.equal(p1, predict(h, newdata = fut))))
  # Groups without history are new groups
  none <- one[0, ]
  expect_equal(predict(h, newdata = fut, history = none),
               predict(h, newdata = fut, re.form = NA))
  expect_error(predict(h, newdata = fut,
                       history = one[, setdiff(names(one), "Total_crashes")]),
               "outcome")
})

test_that("profile intervals agree with glmmTMB and are asymmetric", {
  skip_if_not_installed("glmmTMB")
  skip_on_cran()
  s <- sal()
  f <- count ~ spp + mined + (1 | site)
  h <- countreg_hier(f, data = s, family = "NB2", engine = "laplace")
  ci <- confint(h, method = "profile")
  m <- suppressWarnings(glmmTMB::glmmTMB(f, data = s,
                                         family = glmmTMB::nbinom2))
  ref <- exp(confint(m, parm = "theta_", method = "profile")[, 1:2])
  expect_equal(unname(ci[1, ]), unname(ref), tolerance = 2e-3)
  wald <- confint(h)["sd_(Intercept)|site", ]
  expect_gt(ci[1, 1], wald[1])
  expect_error(confint(h, "sppPR:nope", method = "profile"),
               "not available")
})

test_that("correlated random effects (cre)", {
  data("washington_roads", package = "flexCountReg", envir = environment())
  w <- washington_roads
  h <- countreg_hier(Total_crashes ~ lnaadt + speed50 + offset(lnlength) +
                       (1 | ID), data = w, family = "POISSON", cre = ~ lnaadt)
  expect_true("mean_lnaadt" %in% names(coef(h)))
  tst <- cre_test(h)
  expect_s3_class(tst, "htest")
  nd <- w[1:4, c("ID", "lnaadt", "speed50", "lnlength")]
  expect_equal(predict(h, newdata = nd), fitted(h)[1:4], tolerance = 1e-8)
  expect_error(countreg_hier(Total_crashes ~ lnaadt + (1 | ID) + (1 | Year),
                             data = w, cre = ~ lnaadt), "cre_group")
})

test_that("marginal effects of hierarchical and fixed-effects fits", {
  data("washington_roads", package = "flexCountReg", envir = environment())
  w <- washington_roads
  h <- countreg_hier(Total_crashes ~ lnaadt + speed50 + offset(lnlength) +
                       (1 | ID), data = w, family = "POISSON")
  me <- margEffTable(h)
  expect_setequal(me$term, c("lnaadt", "speed50"))
  expect_true(all(me$std_error > 0))
  fe <- suppressWarnings(suppressMessages(countreg_fe(
    Total_crashes ~ lnaadt + speed50 + offset(lnlength), fixef = ~ ID,
    data = w)))
  expect_identical(margEffTable(fe)$term, "lnaadt")
})
