# User-defined count families (new_count_family(), check_count_family()).

user_nb2 <- function(name = "myNB2") {
  new_count_family(
    name, params = list(alpha = count_param("log", start = -1)),
    logpmf = function(y, mu, par)
      stats::dnbinom(y, size = 1 / par$alpha, mu = mu, log = TRUE))
}

wr <- function() {
  data("washington_roads", package = "flexCountReg", envir = environment())
  washington_roads
}

f_main <- Total_crashes ~ lnaadt + lnlength

test_that("the constructor validates its inputs", {
  expect_error(new_count_family("", function(y, mu, par) 0), "name")
  expect_error(new_count_family("a", "not a function"), "logpmf")
  expect_error(new_count_family("a", function(y, mu, par) 0,
                                params = list(count_param(), count_param())),
               "named list")
  expect_error(new_count_family(
    "a", function(y, mu, par) 0,
    params = list(a = count_param(), b = count_param(), c = count_param())),
    "at most two")
  expect_error(new_count_family("a", function(y, mu, par) 0,
                                params = list(a = 1)), "count_param")
  expect_error(count_param(start = NA), "start")
  expect_error(count_param("probit"), "should be one of")
  fam <- user_nb2()
  expect_s3_class(fam, "fcr_count_family")
  expect_match(fam$key, "^user:myNB2:[0-9a-f]{12}$")
  expect_output(print(fam), "myNB2.*alpha \\(log link\\)")
})

test_that("a user NB2 reproduces the built-in NB2", {
  w <- wr()
  nb <- user_nb2()
  a <- countreg(f_main, data = w, family = nb, method = "BFGS")
  b <- countreg(f_main, data = w, family = "NB2", method = "BFGS")
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)),
               tolerance = 1e-8)
  expect_equal(coef(a), coef(b), tolerance = 1e-4)
  expect_identical(names(coef(a)), names(coef(b)))
  expect_identical(a$model$family_spec, nb)
  expect_equal(predict(a, newdata = w[1:5, ], type = "prob"),
               predict(b, newdata = w[1:5, ], type = "prob"),
               tolerance = 1e-5)
  set.seed(3)
  u <- stats::runif(nrow(w))
  expect_equal(.fcr_quantile_residuals(a, u), .fcr_quantile_residuals(b, u),
               tolerance = 1e-4)
  expect_true(all(is.finite(as.matrix(simulate(a, nsim = 2, seed = 1)))))
  expect_output(print(summary(a)), "family: myNB2")
  # Variance and mean by summing the PMF
  pv <- .fcr_family_variance(a$model$family, c(0.5, 4), c(0.4, 0.4))
  expect_equal(pv, c(0.5, 4) + 0.4 * c(0.5, 4)^2, tolerance = 1e-8)
  # Dispersion formula and hurdle model
  a <- countreg(f_main, data = w, family = nb, method = "BFGS",
                dis_param_formula_1 = ~ speed50)
  b <- countreg(f_main, data = w, family = "NB2", method = "BFGS",
                dis_param_formula_1 = ~ speed50)
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)),
               tolerance = 1e-8)
  a <- countreg(f_main, data = w, family = nb, method = "BFGS",
                zero_formula = ~ speed50, zero_type = "hurdle")
  b <- countreg(f_main, data = w, family = "NB2", method = "BFGS",
                zero_formula = ~ speed50, zero_type = "hurdle")
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)),
               tolerance = 1e-8)
  expect_error(countreg(f_main, data = w, family = nb,
                        underreport_formula = ~ speed50),
               "user-defined family \"myNB2\"")
})

test_that("user families work in random-parameter and mixture models", {
  skip_on_cran()
  w <- wr()
  nb <- user_nb2()
  a <- suppressWarnings(countreg_rp(
    Total_crashes ~ lnaadt, rpar_formula = ~ -1 + speed50, data = w,
    family = nb, group_var = "ID", ndraws = 30))
  b <- suppressWarnings(countreg_rp(
    Total_crashes ~ lnaadt, rpar_formula = ~ -1 + speed50, data = w,
    family = "NB2", group_var = "ID", ndraws = 30))
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)),
               tolerance = 1e-6)
  expect_equal(predict(a, newdata = w, method = "Individual"),
               predict(b, newdata = w, method = "Individual"),
               tolerance = 1e-3, ignore_attr = TRUE)
  # Mixtures (a subset and few starts: user families have no analytic
  # gradient in mixtures, so these fits are slower than the built-in NB2)
  ws <- w[1:500, ]
  a <- suppressMessages(suppressWarnings(
    countreg_mix(f_main, data = ws, family = nb, components = 2, seed = 1,
                 n_starts = 2)))
  b <- suppressMessages(suppressWarnings(
    countreg_mix(f_main, data = ws, family = "NB2", components = 2,
                 seed = 1, n_starts = 2)))
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)),
               tolerance = 1e-6)
  # A user family and a built-in one as components
  m <- suppressMessages(suppressWarnings(countreg_mix(
    f_main, data = ws, family = list(nb, "POISSON"), seed = 1,
    n_starts = 2)))
  expect_length(m$model$family_spec, 1L)
  expect_identical(m$model$family[1], nb$key)
})

test_that("two parameters, a logit link and a mean that is not mu", {
  w <- wr()
  zinb <- new_count_family(
    "my-ZINB",
    params = list(alpha = count_param("log", -1),
                  pi = count_param("logit", -2)),
    mean_is_mu = FALSE,
    logpmf = function(y, mu, par) {
      lnb <- stats::dnbinom(y, size = 1 / par$alpha, mu = mu, log = TRUE)
      ifelse(y == 0, log(par$pi + (1 - par$pi) * exp(lnb)),
             log1p(-par$pi) + lnb)
    })
  a <- suppressWarnings(countreg(f_main, data = w, family = zinb,
                                 method = "BFGS"))
  b <- suppressWarnings(countreg(f_main, data = w, family = "NB2",
                                 zero_formula = ~ 1, method = "BFGS"))
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)),
               tolerance = 1e-6)
  expect_identical(names(coef(a))[4:5], c("ln(alpha)", "logit(pi)"))
  expect_equal(coef(a)[["logit(pi)"]], coef(b)[["Zero:(Intercept)"]],
               tolerance = 1e-2)
  # The mean (1 - pi) mu is computed by summing the PMF
  expect_equal(predict(a, newdata = w[1:20, ]),
               predict(b, newdata = w[1:20, ]), tolerance = 1e-3)
})

test_that("families with the same name do not collide", {
  w <- wr()
  nb <- user_nb2("same")
  other <- new_count_family(
    "same", params = list(theta = count_param()),
    logpmf = function(y, mu, par)
      stats::dnbinom(y, size = par$theta, mu = mu, log = TRUE))
  expect_false(identical(nb$key, other$key))
  a <- countreg(f_main, data = w, family = nb, method = "BFGS")
  b <- countreg(f_main, data = w, family = other, method = "BFGS")
  # size = 1 / alpha against size = theta: ln(theta) = -ln(alpha)
  expect_equal(coef(b)[["ln(theta)"]], -coef(a)[["ln(alpha)"]],
               tolerance = 1e-3)
})

test_that("a saved fit brings its family into a new session", {
  w <- wr()
  m <- countreg(f_main, data = w, family = user_nb2("saved"),
                method = "BFGS")
  p <- predict(m, newdata = w[1:3, ])
  # Forget the family, as a new R session would
  rm(list = m$model$family, envir = .fcr_user_families)
  expect_null(.fcr_user_families[[m$model$family]])
  expect_equal(predict(m, newdata = w[1:3, ]), p)
  expect_false(is.null(.fcr_user_families[[m$model$family]]))
  skip_on_cran()
  skip_if_not_installed("callr")
  f <- tempfile(fileext = ".rds")
  saveRDS(m, f)
  p2 <- callr::r(function(f, nd) {
    loadNamespace("flexCountReg")
    m <- readRDS(f)
    stats::predict(m, newdata = nd)
  }, args = list(f = f, nd = w[1:3, ]))
  expect_equal(p2, p)
})

test_that("a user family key alone cannot be resolved", {
  fam <- user_nb2("orphan")
  key <- .fcr_register_family(fam)
  rm(fam)
  rm(list = key, envir = .fcr_user_families)
  expect_error(.fcr_get_family(key), "not available in this R session")
})

test_that("check_count_family() finds wrong supplied functions", {
  good <- new_count_family(
    "NBv", params = list(alpha = count_param()),
    logpmf = function(y, mu, par)
      stats::dnbinom(y, size = 1 / par$alpha, mu = mu, log = TRUE),
    variance = function(mu, par) mu + par$alpha * mu^2,
    cdf = function(q, mu, par) stats::pnbinom(q, size = 1 / par$alpha,
                                              mu = mu),
    rand = function(n, mu, par) stats::rnbinom(n, size = 1 / par$alpha,
                                               mu = mu))
  set.seed(1)
  out <- expect_output(check_count_family(
    good, mu = c(0.5, 10), par = list(alpha = c(0.5, 2))))
  expect_true(all(out$ok))
  bad <- good
  bad$variance <- function(mu, par) mu + par$alpha * mu
  set.seed(1)
  expect_warning(expect_output(check_count_family(
    bad, mu = c(0.5, 10), par = list(alpha = c(0.5, 2)))),
    "variance")
  expect_error(check_count_family(good, par = list()), "alpha")
  expect_error(check_count_family("NB2"), "new_count_family")
})
