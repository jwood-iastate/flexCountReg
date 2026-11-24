# Penalized smooth terms (generalized additive count models).

gam_data <- function() {
  data("washington_roads", package = "flexCountReg", envir = environment())
  washington_roads
}
f_s <- Total_crashes ~ s(lnaadt, k = 10) + speed50 + offset(lnlength)

test_that("smooth bases are built once and rebuilt for new data", {
  skip_if_not_installed("mgcv")
  w <- gam_data()
  p <- prepare_model_data(f_s, w, aux_formulas = list(
    dis1 = ~ s(lnlength, k = 5)))
  expect_identical(colnames(p$X)[3:11], paste0("s(lnaadt).", 1:9))
  expect_identical(p$smooth$main[[1]]$cols, 3:11)
  expect_identical(p$smooth$aux$dis1[[1]]$cols, 2:5)
  rows <- c(5, 100, 1000)
  f <- forge_new_data(p$blueprint, w[rows, ])
  expect_equal(unname(f$X), unname(p$X[rows, ]))
  expect_equal(unname(f$aux$dis1), unname(p$aux$dis1[rows, ]))
  nd <- w[rows, ]
  nd$lnaadt[2] <- NA
  expect_true(all(is.na(forge_new_data(p$blueprint, nd)$X[2, 3:11])))
})

test_that("fixed smoothing parameters reproduce mgcv (Poisson and NB2)", {
  skip_if_not_installed("mgcv")
  w <- gam_data()
  gp <- mgcv::gam(f_s, data = w, family = stats::poisson(), method = "REML")
  a <- countreg(f_s, data = w, family = "POISSON", sp = gp$sp)
  expect_equal(coef(a), coef(gp), tolerance = 1e-5, ignore_attr = TRUE)
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(gp)),
               tolerance = 1e-8)
  expect_equal(sum(a$model$edf[3:11]), sum(gp$edf[3:11]), tolerance = 1e-4)
  expect_equal(unname(sqrt(diag(vcov(a)))), unname(sqrt(diag(vcov(gp)))),
               tolerance = 1e-4)
  # NB2 with the dispersion fixed at our estimate
  b <- countreg(f_s, data = w, family = "NB2", sp = gp$sp)
  gt <- mgcv::gam(f_s, data = w, sp = gp$sp,
                  family = mgcv::negbin(exp(-coef(b)[["ln(alpha)"]])))
  expect_equal(coef(b)[names(coef(gt))], coef(gt), tolerance = 1e-4)
  # Smooth functions and their standard errors equal mgcv's term values
  e <- smooth_estimates(a, n = 40, nsim = 500, seed = 1)
  pg <- stats::predict(gp, newdata = data.frame(lnaadt = e$lnaadt,
                                                speed50 = 0, lnlength = 0),
                       type = "terms", se.fit = TRUE)
  expect_equal(e$fit, unname(pg$fit[, "s(lnaadt)"]), tolerance = 1e-5)
  expect_equal(e$se, unname(pg$se.fit[, "s(lnaadt)"]), tolerance = 1e-5)
  expect_true(all(e$supr - e$fit > e$upr - e$fit))
})

test_that("selected smoothing parameters are close to mgcv's REML", {
  skip_if_not_installed("mgcv")
  w <- gam_data()
  gp <- mgcv::gam(f_s, data = w, family = stats::poisson(), method = "REML")
  a <- countreg(f_s, data = w, family = "POISSON")
  expect_true(a$model$smooth_converged)
  expect_lt(abs(log(a$model$smooth_lambda / gp$sp)), 0.1)
  expect_equal(sum(a$model$edf[3:11]), sum(gp$edf[3:11]), tolerance = 0.02)
  # A second, unrelated smooth (dispersion) does not disturb the first
  b <- countreg(f_s, data = w, family = "NB2",
                dis_param_formula_1 = ~ s(lnlength, k = 6))
  b1 <- countreg(f_s, data = w, family = "NB2")
  st <- b$model$smooth_table
  expect_identical(st$parameter, c("mean", "ln(alpha)"))
  expect_equal(st$edf[1], b1$model$smooth_table$edf[1], tolerance = 0.1)
  # Predictions use the stored bases
  rows <- c(3, 70, 900)
  expect_equal(predict(b, newdata = w[rows, ]), predict(b)[rows])
})

test_that("summary, logLik and covariance of smooth fits", {
  skip_if_not_installed("mgcv")
  w <- gam_data()
  a <- countreg(f_s, data = w, family = "NB2", stderr = "robust")
  s <- summary(a)
  expect_false(any(grepl("^s\\(lnaadt\\)", rownames(s$coefficients))))
  expect_identical(s$smooth_table$term, "s(lnaadt)")
  expect_output(print(s), "Smooth terms")
  expect_equal(attr(logLik(a), "df"), sum(a$model$edf))
  expect_lt(attr(logLik(a), "df"), length(coef(a)))
  expect_identical(a$model$vcov_type, "robust")
  expect_error(countreg(f_s, data = w, family = "NB2", stderr = "boot",
                        bootstraps = 5), "not yet available")
  z <- countreg(Total_crashes ~ s(lnaadt, k = 6) + offset(lnlength),
                data = w, family = "NB2", zero_formula = ~ s(lnlength, k = 5))
  expect_identical(z$model$smooth_table$parameter, c("mean", "Zero"))
})

test_that("smooths in random-parameter models", {
  skip_on_cran()
  skip_if_not_installed("mgcv")
  w <- gam_data()
  set.seed(2)
  w$z <- stats::runif(nrow(w))
  rp <- suppressWarnings(countreg_rp(
    Total_crashes ~ s(lnaadt, k = 8) + offset(lnlength),
    rpar_formula = ~ -1 + speed50, data = w, family = "NB2",
    group_var = "ID", ndraws = 30))
  expect_identical(rp$model$smooth_table$term, "s(lnaadt)")
  rows <- c(2, 50, 700)
  expect_equal(predict(rp, newdata = w[rows, ]), predict(rp)[rows])
  rh <- suppressWarnings(countreg_rp(
    Total_crashes ~ lnaadt + offset(lnlength), rpar_formula = ~ -1 + speed50,
    het_mean_formula = ~ s(z, k = 5), data = w, family = "NB2",
    group_var = "ID", ndraws = 20))
  expect_identical(rh$model$smooth_table$parameter, "HetMean:speed50")
  expect_error(countreg_rp(Total_crashes ~ lnaadt,
                           rpar_formula = ~ -1 + s(speed50), data = w,
                           ndraws = 10), "rpar_formula")
})

test_that("plot_smooths draws every smooth", {
  skip_if_not_installed("mgcv")
  w <- gam_data()
  m <- countreg(Total_crashes ~ te(lnaadt, lnlength, k = c(4, 4)),
                data = w, family = "POISSON",
                sp = c(1, 1))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  e <- plot_smooths(m, nsim = 200)
  expect_identical(unique(e$term), "te(lnaadt,lnlength)")
  expect_error(smooth_estimates(m, term = "s(nope)"), "No smooth term")
  m0 <- countreg(Total_crashes ~ lnaadt, data = w, family = "POISSON",
                 method = "BFGS")
  expect_error(smooth_estimates(m0), "no smooth terms")
})

test_that("smooths in the random-effects panel models", {
  skip_on_cran()
  skip_if_not_installed("mgcv")
  w <- gam_data()
  f <- Total_crashes ~ s(lnaadt, k = 8) + speed50 + offset(lnlength)
  rp <- repois(f, group_var = "ID", data = w)
  nm <- suppressWarnings(negmultinom(f, group_var = "ID", data = w,
                                     family = "NM2"))
  # NM2 and the Poisson-gamma random-effects model are the same model
  expect_equal(nm$model$smooth_table$edf, rp$model$smooth_table$edf,
               tolerance = 1e-3)
  expect_equal(unname(coef(nm)[1:9]), unname(coef(rp)[1:9]),
               tolerance = 1e-3)
  rn <- suppressWarnings(renb(f, group_var = "ID", data = w))
  for (m in list(rp, nm, rn)) {
    expect_identical(m$model$smooth_table$term, "s(lnaadt)")
    expect_identical(m$model$vcov_type, "bayesian")
    rows <- c(4, 80, 1200)
    expect_equal(predict(m, newdata = w[rows, ]),
                 predict(m, newdata = w)[rows])
  }
  expect_error(repois(f, group_var = "ID", data = w, stderr = "boot",
                      bootstraps = 3), "not yet available")
  w$site_len <- stats::ave(w$lnlength, w$ID)
  nd <- suppressWarnings(negmultinom(
    Total_crashes ~ lnaadt + offset(lnlength), group_var = "ID", data = w,
    dis_param_formula_1 = ~ s(site_len, k = 5), stderr = "robust"))
  expect_identical(nd$model$smooth_table$parameter, "ln(alpha)")
})

test_that("the analytic repois score matches the numeric one", {
  w <- gam_data()
  p <- prepare_model_data(Total_crashes ~ lnaadt + speed50 + offset(lnlength),
                          w, group_var = "ID")
  cl <- split(seq_along(p$y), p$group, drop = TRUE)
  par <- c(-8, 1, -0.3, -1.2)
  set.seed(1)
  wt <- sample(1:3, length(cl), TRUE)
  a <- .repois_grad(par, p$y, p$X, p$offset, cl, wt)
  n <- maxLik::numericGradient(function(b)
    .repois_loglik(b, p$y, p$X, p$offset, cl, wt), par)
  expect_lt(max(abs(a - n)) / max(abs(n)), 1e-7)
})

test_that("smooths in fixed-effects models match mgcv with dummies", {
  skip_if_not_installed("mgcv")
  set.seed(7)
  G <- 150
  Tn <- 5
  d <- data.frame(id = factor(rep(seq_len(G), each = Tn)),
                  t = factor(rep(seq_len(Tn), G)))
  a_i <- stats::rnorm(G, 0, 0.6)
  d$x <- stats::runif(nrow(d), -2, 2) + 0.3 * a_i[d$id]
  d$z <- stats::rnorm(nrow(d))
  d$y <- stats::rpois(nrow(d), exp(0.3 + a_i[d$id] + sin(1.5 * d$x) +
                                     0.3 * d$z))
  f <- y ~ s(x, k = 8) + z
  a <- suppressMessages(countreg_fe(f, fixef = ~ id, data = d))
  expect_gt(a$model$smooth_table$edf, 3)       # a curve, not a line
  keep <- d$id %in% names(a$model$fixef$fe1)
  dk <- droplevels(d[keep, ])
  g <- mgcv::gam(y ~ s(x, k = 8) + z + id, data = dk,
                 family = stats::poisson(), sp = a$model$smooth_lambda)
  b <- names(coef(a))
  expect_equal(coef(a), coef(g)[b], tolerance = 1e-6)
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(g)),
               tolerance = 1e-8)
  # Demeaning gives the same fit; two ways match mgcv with both dummies
  dm <- suppressMessages(countreg_fe(f, fixef = ~ id, data = d,
                                     method = "demean",
                                     sp = a$model$smooth_lambda))
  expect_equal(coef(dm), coef(a), tolerance = 1e-6)
  t2 <- suppressMessages(countreg_fe(f, fixef = ~ id + t, data = d,
                                     sp = 1))
  k2 <- d$id %in% names(t2$model$fixef$fe1)
  g2 <- mgcv::gam(y ~ s(x, k = 8) + z + id + t, data = droplevels(d[k2, ]),
                  family = stats::poisson(), sp = 1)
  expect_equal(coef(t2), coef(g2)[b], tolerance = 1e-6)
  # Predictions use the bases built on the estimation rows
  expect_equal(predict(a, newdata = dk[c(3, 30), ]),
               unname(fitted(g)[c(3, 30)]), tolerance = 1e-6)
  e <- smooth_estimates(a, n = 20, nsim = 200, seed = 1)
  expect_identical(unique(e$term), "s(x)")
  # A smooth of a group-constant variable is absorbed by the fixed effects
  d$xbar <- stats::ave(d$x, d$id)
  expect_error(suppressMessages(countreg_fe(y ~ s(xbar) + z, fixef = ~ id,
                                            data = d)),
               "absorbed by the fixed effects")
  expect_error(suppressMessages(countreg_fe(f, fixef = ~ id, data = d,
                                            jackknife = TRUE)),
               "jackknife")
})

test_that("smooths in finite mixtures and latent classes", {
  skip_on_cran()
  skip_if_not_installed("mgcv")
  w <- gam_data()
  f <- Total_crashes ~ s(lnaadt, k = 8) + offset(lnlength)
  # One component is the countreg() model
  m1 <- suppressMessages(countreg_mix(f, data = w, family = "NB2",
                                      components = 1, n_starts = 1,
                                      seed = 1))
  cr <- countreg(f, data = w, family = "NB2")
  expect_equal(as.numeric(logLik(m1)), as.numeric(logLik(cr)),
               tolerance = 1e-6)
  expect_equal(m1$model$smooth_table$edf, cr$model$smooth_table$edf,
               tolerance = 1e-4)
  m2 <- suppressMessages(suppressWarnings(countreg_mix(
    f, data = w, family = "POISSON", components = 2, n_starts = 2,
    seed = 1)))
  st <- m2$model$smooth_table
  expect_identical(st$parameter, c("Component1:mean", "Component2:mean"))
  expect_identical(m2$model$vcov_type, "bayesian")
  expect_equal(attr(logLik(m2), "df"), m2$model$df_total)
  rows <- c(4, 90, 1300)
  expect_equal(predict(m2, newdata = w[rows, ]), predict(m2)[rows])
  expect_output(print(summary(m2)), "Smooth terms")
  lc <- suppressMessages(suppressWarnings(countreg_lc(
    Total_crashes ~ lnaadt + offset(lnlength), data = w, family = "POISSON",
    classes = 2, membership_formula = ~ s(lnlength, k = 5), n_starts = 2,
    seed = 1)))
  expect_match(lc$model$smooth_table$parameter, "^Mixing:")
  expect_error(suppressMessages(countreg_mix(
    f, data = w, family = "POISSON", components = 2, n_starts = 1, seed = 1,
    stderr = "boot", bootstraps = 3)), "not yet available")
})

test_that("only calls are smooth terms, not variables named s", {
  expect_false(.fcr_has_smooth(~ -1 + s))
  expect_false(.fcr_has_smooth(y ~ x + te))
  expect_true(.fcr_has_smooth(~ s(x)))
  expect_true(.fcr_has_smooth(y ~ log(x) + te(a, b)))
  expect_false(.fcr_has_smooth(NULL))
  d <- data.frame(y = rpois(50, 2), x = rnorm(50), s = runif(50))
  p <- prepare_model_data(y ~ x + s, d, aux_formulas = list(dis1 = ~ s))
  expect_false(.fcr_prep_has_smooth(p))
})
