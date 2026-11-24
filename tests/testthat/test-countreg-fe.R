# Fixed-effects count regression (countreg_fe()).

fe_wr <- function() {
  data("washington_roads", package = "flexCountReg", envir = environment())
  w <- washington_roads
  set.seed(1)
  w$x2 <- stats::rnorm(nrow(w))
  w$yr <- factor(sample(2016:2018, nrow(w), TRUE))
  w$ID <- as.character(w$ID)
  w
}

nb_panel <- function(G = 150, Tn = 5, seed = 11) {
  set.seed(seed)
  d <- data.frame(id = rep(seq_len(G), each = Tn), t = rep(seq_len(Tn), G))
  a_i <- stats::rnorm(G, 0, 0.7)
  d$x1 <- stats::rnorm(nrow(d)) + 0.5 * a_i[d$id]
  d$x2 <- stats::rnorm(nrow(d))
  d$y <- stats::rnbinom(nrow(d), size = 2,
                        mu = exp(0.5 + a_i[d$id] + 0.6 * d$x1 - 0.4 * d$x2))
  d
}

test_that("Poisson FE equals the dummy-variable Poisson (one and two ways)", {
  w <- fe_wr()
  f <- Total_crashes ~ lnaadt + x2 + offset(lnlength)
  a <- suppressMessages(countreg_fe(f, fixef = ~ ID, data = w))
  expect_s3_class(a, "countreg_fe")
  expect_identical(a$model$method, "conditional")
  keep <- w$ID %in% names(a$model$fixef$fe1)
  g <- stats::glm(Total_crashes ~ lnaadt + x2 + factor(ID) +
                    offset(lnlength) - 1, data = w[keep, ],
                  family = stats::poisson())
  b <- c("lnaadt", "x2")
  expect_equal(coef(a), coef(g)[b], tolerance = 1e-10)
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(g)),
               tolerance = 1e-10)
  expect_identical(attr(logLik(a), "df"), attr(logLik(g), "df"))
  expect_equal(vcov(a, type = "hessian"), vcov(g)[b, b], tolerance = 1e-8)
  # Cluster-robust by site: the sandwich of the dummy model, with G the
  # number of sites in the fit
  G <- length(unique(w$ID[keep]))
  v0 <- sandwich::vcovCL(g, cluster = w$ID[keep], type = "HC0",
                         cadjust = FALSE)[b, b]
  expect_equal(vcov(a), v0 * G / (G - 1), tolerance = 1e-6,
               ignore_attr = TRUE)
  expect_equal(predict(a, newdata = w[keep, ]), unname(fitted(g)),
               tolerance = 1e-8)
  # Dropped all-zero and single-observation sites are recorded
  expect_gt(length(a$model$dropped$ID), 0)
  expect_false(any(a$model$dropped$ID %in% names(a$model$fixef$fe1)))
  # Demeaning gives the same fit
  d <- suppressMessages(countreg_fe(f, fixef = ~ ID, data = w,
                                    method = "demean"))
  expect_equal(coef(d), coef(a), tolerance = 1e-10)
  expect_equal(vcov(d), vcov(a), tolerance = 1e-8)
  # Two ways
  t2 <- suppressMessages(countreg_fe(f, fixef = ~ ID + yr, data = w))
  expect_identical(t2$model$method, "demean")
  k2 <- w$ID %in% names(t2$model$fixef$fe1)
  g2 <- stats::glm(Total_crashes ~ lnaadt + x2 + factor(ID) + yr +
                     offset(lnlength) - 1, data = w[k2, ],
                   family = stats::poisson())
  expect_equal(coef(t2), coef(g2)[b], tolerance = 1e-8)
  expect_equal(as.numeric(logLik(t2)), as.numeric(logLik(g2)),
               tolerance = 1e-8)
  expect_equal(predict(t2, newdata = w[k2, ]), unname(fitted(g2)),
               tolerance = 1e-7)
})

test_that("NB2 FE equals the dummy-variable NB2 (glm.nb)", {
  d <- nb_panel()
  n <- suppressMessages(countreg_fe(y ~ x1 + x2, fixef = ~ id, data = d,
                                    family = "NB2"))
  expect_identical(n$model$method, "concentrated")
  keep <- d$id %in% names(n$model$fixef$fe1)
  gn <- MASS::glm.nb(y ~ x1 + x2 + factor(id) - 1, data = d[keep, ])
  expect_equal(coef(n)[1:2], coef(gn)[c("x1", "x2")], tolerance = 1e-4)
  expect_equal(coef(n)[["ln(alpha)"]], -log(gn$theta), tolerance = 1e-3)
  expect_equal(as.numeric(logLik(n)), as.numeric(logLik(gn)),
               tolerance = 1e-6)
  expect_equal(predict(n, newdata = d[keep, ]), unname(fitted(gn)),
               tolerance = 1e-3)
  # glm.nb() reports the expected information; with the same weights
  # mu / (1 + alpha mu) our estimates give the same standard errors
  mu <- predict(n, newdata = d[keep, ])
  W <- mu / (1 + exp(coef(n)[["ln(alpha)"]]) * mu)
  Xt <- .fcr_fe_demean_cols(as.matrix(d[keep, c("x1", "x2")]),
                            list(factor(d$id[keep])), W)
  expect_equal(sqrt(diag(solve(crossprod(Xt * sqrt(W))))),
               sqrt(diag(vcov(gn)))[1:2], tolerance = 1e-3,
               ignore_attr = TRUE)
  # A user-defined NB2 gives the same fit
  nbu <- new_count_family(
    "feNB2", params = list(alpha = count_param("log", -1)),
    logpmf = function(y, mu, par)
      stats::dnbinom(y, size = 1 / par$alpha, mu = mu, log = TRUE))
  u <- suppressMessages(countreg_fe(y ~ x1 + x2, fixef = ~ id, data = d,
                                    family = nbu))
  expect_equal(unname(coef(u)), unname(coef(n)), tolerance = 1e-4)
})

test_that("Poisson by concentration equals the conditional likelihood", {
  d <- nb_panel()
  a <- suppressMessages(countreg_fe(y ~ x1 + x2, fixef = ~ id, data = d))
  b <- suppressMessages(countreg_fe(y ~ x1 + x2, fixef = ~ id, data = d,
                                    method = "concentrated"))
  expect_equal(coef(b), coef(a), tolerance = 1e-6)
  expect_equal(as.numeric(logLik(b)), as.numeric(logLik(a)),
               tolerance = 1e-8)
  expect_equal(vcov(b), vcov(a), tolerance = 1e-4)
})

test_that("predictions, new groups and the jackknife", {
  d <- nb_panel()
  n <- suppressMessages(countreg_fe(y ~ x1 + x2, fixef = ~ id, data = d,
                                    family = "NB2"))
  nd <- d[1:2, ]
  nd$id <- 99999
  expect_message(p <- predict(n, newdata = nd), "not in the estimation")
  expect_true(all(is.na(p)))
  pm <- predict(n, newdata = nd, newgroup = "mean_fe")
  expect_true(all(is.finite(pm)))
  expect_error(predict(n, newdata = nd, type = "prob"), "response")
  jk <- suppressMessages(countreg_fe(y ~ x1 + x2, fixef = ~ id, data = d,
                                     family = "NB2", jackknife = TRUE,
                                     time = "t"))
  j <- jk$model$jackknife
  expect_equal(coef(jk), 2 * j$uncorrected - (j$half1 + j$half2) / 2)
  expect_equal(j$uncorrected, coef(n))
  d3 <- d[d$t <= 3, ]
  expect_warning(suppressMessages(countreg_fe(
    y ~ x1 + x2, fixef = ~ id, data = d3, family = "NB2",
    jackknife = TRUE, time = "t")), "at least 4")
})

test_that("inputs are checked and unidentified regressors removed", {
  d <- nb_panel()
  d$z <- ave(d$x1, d$id)          # constant within groups
  expect_warning(m <- suppressMessages(countreg_fe(
    y ~ x1 + z, fixef = ~ id, data = d)), "not identified")
  expect_identical(m$model$not_identified, "z")
  expect_identical(names(coef(m)), "x1")
  expect_error(countreg_fe(y ~ x1, fixef = "id", data = d), "one-sided")
  expect_error(countreg_fe(y ~ x1, fixef = ~ nope, data = d), "not found")
  d$t2 <- factor(d$t)
  expect_error(countreg_fe(y ~ x1, fixef = ~ id + t2, data = d,
                           family = "NB2"), "one fixed-effect dimension")
  expect_error(countreg_fe(y ~ x1, fixef = ~ id, data = d,
                           family = "NB2", method = "conditional"),
               "Poisson")
  out <- suppressMessages(countreg_fe(y ~ x1, fixef = ~ id, data = d,
                                      stderr = "hessian"))
  expect_identical(out$model$vcov_type, "hessian")
  expect_output(print(summary(out)), "Poisson fixed effects")
})

test_that("estimates and clustered SEs equal fixest", {
  skip_if_not_installed("fixest")
  w <- fe_wr()
  f <- Total_crashes ~ lnaadt + x2 + offset(lnlength)
  for (fx in c("ID", "ID + yr")) {
    a <- suppressMessages(countreg_fe(f, fixef = stats::as.formula(
      paste("~", fx)), data = w))
    fe <- suppressMessages(fixest::fepois(
      stats::as.formula(paste("Total_crashes ~ lnaadt + x2 +",
                              "offset(lnlength) |", fx)), data = w))
    v <- stats::vcov(fe, cluster = ~ ID,
                      ssc = fixest::ssc(adj = FALSE, cluster.adj = TRUE))
    expect_equal(coef(a), coef(fe), tolerance = 1e-10)
    expect_equal(as.numeric(logLik(a)), as.numeric(logLik(fe)),
                 tolerance = 1e-10)
    expect_equal(vcov(a), v, tolerance = 1e-6, ignore_attr = TRUE)
  }
  d <- nb_panel()
  n <- suppressMessages(countreg_fe(y ~ x1 + x2, fixef = ~ id, data = d,
                                    family = "NB2"))
  fn <- suppressMessages(fixest::fenegbin(y ~ x1 + x2 | id, data = d))
  expect_equal(coef(n)[1:2], coef(fn)[c("x1", "x2")], tolerance = 1e-4)
  expect_equal(exp(-coef(n)[["ln(alpha)"]]), fn$theta, tolerance = 1e-4,
               ignore_attr = TRUE)
  expect_equal(as.numeric(logLik(n)), as.numeric(logLik(fn)),
               tolerance = 1e-8)
})
